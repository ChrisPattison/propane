include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
string(TIMESTAMP BUILD_TIME UTC)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/version.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/version.hpp"
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(dSFMT)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/dSFMT)

add_executable(
    propane 
    main.cpp 
    graph.cpp 
    random_number_generator.cpp 
    parse.cpp 
    output.cpp
    parallel_population_annealing.cpp 
    parallel.cpp 
    population_annealing.cpp 
    string_util.cpp 
    parallel_types.cpp 
    compare.cpp
    log_lookup.cpp
)

install(TARGETS propane 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)

set_target_properties(propane PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
set(LIBS ${LIBS} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(propane PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(propane PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set(LIBS ${LIBS} dSFMT)

target_link_libraries(propane ${LIBS})
