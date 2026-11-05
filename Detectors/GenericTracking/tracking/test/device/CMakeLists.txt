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

# The device tests on each backend, labelled gpu because they need one.

if(NOT CUDA_ENABLED AND NOT HIP_ENABLED)
  return()
endif()

set(fixtures "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/.." "${CMAKE_CURRENT_SOURCE_DIR}/../../../study/include")

# The test kernels' host code: what O2 compiles only outside GPU translation
# units (a FlatObject's cloning and relocation). An object library of its own,
# since the HIP libraries compile every source they hold as HIP.
add_library(gtrk_device_test_host OBJECT DeviceTestHost.cxx)
target_include_directories(gtrk_device_test_host PRIVATE ${fixtures})
target_link_libraries(gtrk_device_test_host PRIVATE O2::GTrkTracking)
set_target_properties(gtrk_device_test_host PROPERTIES POSITION_INDEPENDENT_CODE ON CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN YES)

function(gtrk_device_test name source backend)
  string(TOLOWER ${backend} lower)
  o2_add_test(${name}-${lower} SOURCES ${source} COMPONENT_NAME generic-tracking NAME generic-tracking-${name}-${lower}
              PUBLIC_LINK_LIBRARIES O2::GTrkTracking${backend} gtrk::device_test_${lower} O2::GTrkTracking O2::GPUTracking
              LABELS generic-tracking gpu ${ARGN} TARGETVARNAME target)
  if(target)
    target_include_directories(${target} PRIVATE ${fixtures})
    target_compile_definitions(${target} PRIVATE GTRK_TEST_${backend})
  endif()
endfunction()

foreach(backend IN ITEMS CUDA HIP)
  if(NOT ${backend}_ENABLED)
    continue()
  endif()
  string(TOLOWER ${backend} lower)

  gtrk_device_sources(sources ${backend} DeviceTestKernels.cu)
  add_library(gtrk_device_test_${lower} SHARED ${sources})
  add_library(gtrk::device_test_${lower} ALIAS gtrk_device_test_${lower})
  target_sources(gtrk_device_test_${lower} PRIVATE $<TARGET_OBJECTS:gtrk_device_test_host>)
  target_include_directories(gtrk_device_test_${lower} PRIVATE ${fixtures})
  target_link_libraries(gtrk_device_test_${lower} PRIVATE O2::GTrkTracking)
  gtrk_device_target(gtrk_device_test_${lower} ${backend} SHARED "_ZN2o24gtrk*testKernels${backend}E*")

  gtrk_device_test(device-runtime testDeviceRuntime.cxx ${backend})
  gtrk_device_test(gpu-stage testGpuStage.cxx ${backend})
  gtrk_device_test(flat-stage testFlatStage.cxx ${backend})
  gtrk_device_test(device-tracker testDeviceTracker.cxx ${backend})
  gtrk_device_test(slice-glue testSliceGlue.cxx ${backend})
  gtrk_device_test(device-storage testDeviceStorage.cxx ${backend})
  gtrk_device_test(road-glue testRoadGlue.cxx ${backend})
  gtrk_device_test(selection-glue testSelectionGlue.cxx ${backend})
  # O2's material LUT on the device, if GTRK_TEST_MATLUT names one.
  gtrk_device_test(device-material testDeviceMaterial.cxx ${backend})
endforeach()
