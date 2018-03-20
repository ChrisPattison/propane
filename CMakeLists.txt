cmake_minimum_required(VERSION 3.7)

project (PROPANE CXX)
set (PROPANE_VERSION_MAJOR 2)
set (PROPANE_VERSION_MINOR 1)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Enable interprocedural optimization if supported
set(IPO_SUPPORTED NO)
if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.9")
    cmake_policy(VERSION 3.9)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT IPO_SUPPORTED)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION IPO_SUPPORTED)
endif()
message("* Using interprocedural optimizations : ${IPO_SUPPORTED}")


set(SHARED_FAST_FLAGS "-march=native")

# -------- Fast
if(NOT CMAKE_CXX_FLAGS_FAST)
    set(CMAKE_CXX_FLAGS_FAST "-Wall ${SHARED_FAST_FLAGS} -O3 -DNDEBUG" CACHE STRING
        "Flags used by the C++ compiler during fast release builds."
        FORCE)
endif()

if(NOT CMAKE_C_FLAGS_FAST)
    set( CMAKE_C_FLAGS_FAST "-Wall ${SHARED_FAST_FLAGS} -O3 -DNDEBUG" CACHE STRING
        "Flags used by the C compiler during fast release builds."
        FORCE )
endif()

if(NOT CMAKE_EXE_LINKER_FLAGS_FAST)
    set( CMAKE_EXE_LINKER_FLAGS_FAST
        "" CACHE STRING
        "Flags used for linking binaries during fast release builds."
        FORCE )
endif()

if(NOT CMAKE_SHARED_LINKER_FLAGS_FAST)
    set( CMAKE_SHARED_LINKER_FLAGS_FAST
        "" CACHE STRING
        "Flags used by the shared libraries linker during fast release builds."
        FORCE )
endif()
# -------- Fast Debug

if(NOT CMAKE_CXX_FLAGS_FASTDEBUG)
    set( CMAKE_CXX_FLAGS_FASTDEBUG "-Wall ${SHARED_FAST_FLAGS} -g -O2" CACHE STRING 
        "Flags used by the C++ compiler during fast debugging builds." 
        FORCE )
endif()

if(NOT CMAKE_C_FLAGS_FASTDEBUG)
    set( CMAKE_C_FLAGS_FASTDEBUG "-Wall ${SHARED_FAST_FLAGS} -g -O2" CACHE STRING
        "Flags used by the C compiler during fast debugging builds."
        FORCE )
endif()

if(NOT CMAKE_EXE_LINKER_FLAGS_FASTDEBUG)
    set( CMAKE_EXE_LINKER_FLAGS_FASTDEBUG
        "" CACHE STRING
        "Flags used for linking binaries during fast debugging builds."
        FORCE )
endif()

if(NOT CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG)
    set( CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG
        "" CACHE STRING
        "Flags used by the shared libraries linker during fast debugging builds."
        FORCE )
endif()

mark_as_advanced (
    CMAKE_CXX_FLAGS_FASTDEBUG
    CMAKE_C_FLAGS_FASTDEBUG
    CMAKE_EXE_LINKER_FLAGS_FASTDEBUG
    CMAKE_SHARED_LINKER_FLAGS_FASTDEBUG )

# -------- Build Type

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Fast
        CACHE STRING "Choose the type of build : None Debug Release RelWithDebInfo MinSizeRel FastDebug Fast."
       FORCE)
endif()

message("* Current build type is : ${CMAKE_BUILD_TYPE}")

subdirs(src doc test)

message("* Build docs with 'make doc'")