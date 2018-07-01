find_package(GTest)

if(GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})

    add_executable(
        tests
        population_annealing_test.cpp
        main.cpp
    )

    set_target_properties(tests PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
    )

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../src/dSFMT)

    set(LIBS psqa_lib ${GTEST_LIBRARIES})
    target_link_libraries(tests ${LIBS})
endif()
