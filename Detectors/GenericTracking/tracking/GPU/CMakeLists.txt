# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# The device tracker for CUDA and, hipified from the same sources, for HIP:
# relocatable device code, device-linked with O2's out-of-line device
# functions (MatLayerCylSet's LUT query and what it calls, from O2's GPU
# external provider) into one shared library per backend, which O2's GPU
# backend links; loaded with it, the library registers its runtime
# (RuntimeFactory.h).

if(NOT CUDA_ENABLED AND NOT HIP_ENABLED)
  return()
endif()

# gtrk_hipify(<output-variable> <source.cu>...): one .hip in the build tree
# per .cu source; only the .cu files are translated, so the headers they
# include stay vendor-neutral.
function(gtrk_hipify output)
  set(generated)
  foreach(source IN LISTS ARGN)
    get_filename_component(absolute "${source}" ABSOLUTE)
    get_filename_component(name "${source}" NAME_WE)
    set(hip "${CMAKE_CURRENT_BINARY_DIR}/hipified/${name}.hip")
    add_custom_command(OUTPUT "${hip}"
                       COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_BINARY_DIR}/hipified"
                       COMMAND ${hip_HIPIFY_PERL_EXECUTABLE} -quiet-warnings "-o=${hip}" "${absolute}"
                       DEPENDS "${absolute}"
                       COMMENT "Hipifying ${source}"
                       VERBATIM)
    list(APPEND generated "${hip}")
  endforeach()
  set(${output} ${generated} PARENT_SCOPE)
endfunction()

# gtrk_device_sources(<output-variable> <CUDA|HIP> <source.cu>...): the
# sources as written for CUDA, hipified for HIP.
function(gtrk_device_sources output backend)
  set(sources ${ARGN})
  if(backend STREQUAL "HIP")
    gtrk_hipify(sources ${ARGN})
    set_source_files_properties(${sources} PROPERTIES LANGUAGE HIP)
  endif()
  set(${output} ${sources} PARENT_SCOPE)
endfunction()

# gtrk_device_target(<target> <CUDA|HIP> <SHARED|OBJECT> [<exported symbol pattern>...]):
# a library of device code for one backend, device-linked with O2's external
# provider as ITS's is, but without device LTO, which moves results in their
# last bits. The CUDA runtime is the
# shared one, so every library of a process talks to the same runtime state;
# cross-execution-space calls are errors. A shared library exports only its
# entries: ITS's GPU library carries the provider too, and its relocatable
# device code's registration symbols would otherwise bind to one copy in a
# process that loads both.
set(GTRK_DEVICE_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/cuda" CACHE INTERNAL "")
function(gtrk_device_target target backend type)
  # The device headers next to the CUDA sources: a hipified .hip lives in the
  # build tree, so a quoted include does not find them beside it.
  target_include_directories(${target} PRIVATE "${CMAKE_SOURCE_DIR}" "${GTRK_DEVICE_HEADERS}")
  string(TOLOWER "${backend}" backend_namespace)
  target_compile_definitions(${target} PRIVATE GTRK_BACKEND_NAMESPACE=${backend_namespace} GTRK_BACKEND_${backend})
  set_target_gpu_arch("${backend}" ${target})
  # The architecture's parameters (GPUCA_WARP_SIZE): GPUDefParametersWrapper.h.
  target_link_libraries(${target} PRIVATE O2::GPUTracking)
  set_target_properties(${target} PROPERTIES POSITION_INDEPENDENT_CODE ON)
  if(type STREQUAL "SHARED")
    target_link_libraries(${target} PRIVATE O2::GPUTracking${backend}ExternalProvider)
    set(version_script "${CMAKE_CURRENT_BINARY_DIR}/${target}.exports")
    list(JOIN ARGN "; " exported)
    if(exported)
      file(GENERATE OUTPUT "${version_script}" CONTENT "{ global: ${exported}; local: *; };\n")
    else()
      file(GENERATE OUTPUT "${version_script}" CONTENT "{ local: *; };\n")
    endif()
    target_link_options(${target} PRIVATE "LINKER:--version-script=${version_script}")
    set_property(TARGET ${target} APPEND PROPERTY LINK_DEPENDS "${version_script}")
  endif()
  if(backend STREQUAL "CUDA")
    set_target_properties(${target} PROPERTIES CUDA_SEPARABLE_COMPILATION ON CUDA_RUNTIME_LIBRARY Shared)
    if(type STREQUAL "SHARED")
      set_target_properties(${target} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    endif()
    target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Werror=cross-execution-space-call>)
  else()
    # FindO2GPU compiles HIP with -ffast-math outside the deterministic mode,
    # and that includes -ffinite-math-only: comparisons with an infinite
    # start value or a NaN guard become undefined. Only the finite-math
    # assumption is taken out (CUDA's -use_fast_math makes none).
    target_compile_options(${target} PRIVATE $<$<COMPILE_LANGUAGE:HIP>:-fgpu-rdc -fno-finite-math-only>)
    # hipCUB, CUB's interface on rocPRIM, for the scans.
    target_link_libraries(${target} PRIVATE hip::hipcub)
    if(type STREQUAL "SHARED")
      target_link_options(${target} PRIVATE -fgpu-rdc --hip-link)
      target_link_libraries(${target} PRIVATE hip::host)
    endif()
  endif()
endfunction()

foreach(backend IN ITEMS CUDA HIP)
  if(NOT ${backend}_ENABLED)
    continue()
  endif()
  string(TOLOWER ${backend} lower)

  # One device tracker per layout and unit, compiled in parallel.
  gtrk_device_sources(sources ${backend} cuda/DeviceRuntime.cu cuda/DeviceTrackerITS.cu cuda/DeviceTrackerMFT.cu cuda/DeviceTrackerTRKFT3.cu)
  o2_add_library(GTrkTracking${backend}
                 SOURCES ${sources}
                 PUBLIC_LINK_LIBRARIES O2::GTrkTracking
                 TARGETVARNAME targetName)
  gtrk_device_target(${targetName} ${backend} SHARED)
endforeach()
