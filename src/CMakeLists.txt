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
    psqa_solve 
    main.cpp 
    graph.cpp 
    random_number_generator.cpp 
    parse.cpp 
    output.cpp
    population_annealing.cpp 
    string_util.cpp 
    compare.cpp
    log_lookup.cpp
)

install(TARGETS psqa_solve 
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
)

set_target_properties(psqa_solve PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

set(Boost_USE_STATIC_LIBS ON)
find_package(Boost REQUIRED COMPONENTS program_options)
set(LIBS ${LIBS} ${Boost_LIBRARIES})
include_directories(${Boost_INCLUDE_DIRS})

set(LIBS ${LIBS} dSFMT)

target_link_libraries(psqa_solve ${LIBS})
