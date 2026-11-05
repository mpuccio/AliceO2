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
# Every unit that instantiates the header tracker -- the tracking device's, the
# study library's, the sink's and the benches' -- compiles at -O3 in
# RelWithDebInfo, where the tracker runs faster than at O2's -O2.
set(GTRK_TRACKER_OPTIMIZATION "$<$<AND:$<CONFIG:RelWithDebInfo>,$<COMPILE_LANGUAGE:CXX>>:-O3>")

add_subdirectory(base)
add_subdirectory(tracking)
add_subdirectory(workflow)
if(ENABLE_UPGRADES)
  add_subdirectory(tools)
  add_subdirectory(study)
endif()
