cmake_minimum_required(VERSION 3.0)

project (PROPANE CXX)
set (PROPANE_VERSION_MAJOR 2)
set (PROPANE_VERSION_MINOR 0)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

SET( CMAKE_CXX_FLAGS_FASTDEBUG "-Wall -march=native -g -O2" CACHE STRING 
    "Flags used by the C++ compiler during fast debugging builds." 
    FORCE )

SET( CMAKE_C_FLAGS_FASTDEBUG "-Wall -pedantic -g -O2" CACHE STRING
    "Flags used by the C compiler during fast debugging builds."
    FORCE )

SET( CMAKE_EXE_LINKER_FLAGS_FASTDEBUG
    "" CACHE STRING
    "Flags used for linking binaries during fast debugging builds."
    FORCE )

SET( CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG
    "" CACHE STRING
    "Flags used by the shared libraries linker during fast debugging builds."
    FORCE )
MARK_AS_ADVANCED (
    CMAKE_CXX_FLAGS_FASTDEBUG
    CMAKE_C_FLAGS_FASTDEBUG
    CMAKE_EXE_LINKER_FLAGS_FASTDEBUG
    CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE FastDebug
        CACHE STRING "Choose the type of build : None Debug Release RelWithDebInfo MinSizeRel FastDebug."
       FORCE)
endif()

message("* Current build type is : ${CMAKE_BUILD_TYPE}")

subdirs(src doc test)

message("* Build docs with 'make doc'")