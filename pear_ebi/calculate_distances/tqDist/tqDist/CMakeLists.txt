CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

## make Visual Studio behave nicely
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/tqDist")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/tqDist")
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
    string( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
    set( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${PROJECT_SOURCE_DIR}/bin")
    set( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${PROJECT_SOURCE_DIR}/tqDist")
    set( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} "${PROJECT_SOURCE_DIR}/tqDist")
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

include_directories(${tqDist_SOURCE_DIR})

if(NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long -m64 -O3")
endif(NOT WIN32)

if(WIN32)
  set(NO_N4_128 1)
endif(WIN32)

if(WIN32)
  set(extension ".dll")
else(WIN32)
  set(extension ".so")
endif(WIN32)

if(NOT NO_N4_128)
  add_definitions(-DN4INT128)
endif(NOT NO_N4_128)

if(NOT NOEXTRACT)
  add_definitions(-DdoExtractAndContract)
endif(NOT NOEXTRACT)

if(CONTRACT_NUM)
  add_definitions("-DCONTRACT_MAX_EXTRA_SIZE=$(CONTRACT_NUM)")
endif(CONTRACT_NUM)

message(STATUS "CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

## LIBRARIES
add_library(tqDist_triplet SHARED
    HDT.cpp
    HDTCountingCCToC.cpp
    HDTCountingCTransform.cpp
    HDTCountingG.cpp
    HDTCountingIGToC.cpp
    HDTCountingLeaf.cpp
    HDTFactory.cpp
    HDTListUtils.cpp
    NewickParser.cpp
    RootedTree.cpp
    RootedTreeFactory.cpp
    AbstractDistanceCalculator.cpp
    TripletDistanceCalculator.cpp
    int_stuff.cpp)
set_target_properties(tqDist_triplet PROPERTIES PREFIX "lib" SUFFIX ${extension})

add_library(tqDist_quartet SHARED
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
set_target_properties(tqDist_quartet PROPERTIES COMPILE_FLAGS -DquartetsToo)
set_target_properties(tqDist_quartet PROPERTIES PREFIX "lib" SUFFIX ${extension})

## EXECUTABLES
add_executable(quartet_dist quartet_dist.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
set_target_properties(quartet_dist PROPERTIES COMPILE_FLAGS -DquartetsToo)

add_executable(all_pairs_quartet_dist all_pairs_quartet_distance.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
set_target_properties(all_pairs_quartet_dist PROPERTIES COMPILE_FLAGS -DquartetsToo)

add_executable(pairs_quartet_dist pairs_quartet_distance.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
set_target_properties(pairs_quartet_dist PROPERTIES COMPILE_FLAGS -DquartetsToo)

add_executable(triplet_dist triplet_dist.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  TripletDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)

add_executable(all_pairs_triplet_dist all_pairs_triplet_distance.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  TripletDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)

add_executable(pairs_triplet_dist pairs_triplet_distance.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  TripletDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)

## TESTS
add_executable(test_triplet test_triplet.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  TripletDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
add_test(
  NAME test_triplet 
  WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_triplet)

add_executable(test_quartet test_quartet.cpp
  HDT.cpp
  HDTCountingCCToC.cpp
  HDTCountingCTransform.cpp
  HDTCountingG.cpp
  HDTCountingIGToC.cpp
  HDTCountingLeaf.cpp
  HDTFactory.cpp
  HDTListUtils.cpp
  NewickParser.cpp
  RootedTree.cpp
  RootedTreeFactory.cpp
  AbstractDistanceCalculator.cpp
  QuartetDistanceCalculator.cpp
  int_stuff.cpp)
set_target_properties(test_quartet PROPERTIES COMPILE_FLAGS -DquartetsToo)
add_test(
  NAME test_quartet 
  WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_quartet)

## Installation
install(TARGETS triplet_dist RUNTIME DESTINATION bin)
install(TARGETS quartet_dist RUNTIME DESTINATION bin)
install(TARGETS pairs_triplet_dist RUNTIME DESTINATION bin)
install(TARGETS pairs_quartet_dist RUNTIME DESTINATION bin)
install(TARGETS all_pairs_triplet_dist RUNTIME DESTINATION bin)
install(TARGETS all_pairs_quartet_dist RUNTIME DESTINATION bin)
if(WIN32)
  install(FILES VCredist/msvcp120.dll VCredist/msvcr120.dll DESTINATION bin)
endif(WIN32)
