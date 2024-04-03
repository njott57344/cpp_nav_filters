cmake_minimum_required(VERSION 3.16.3)

project(CppNavFiltUnitTests)
set(CMAKE_CXX_STANDARD 17)

if(BUILD_TESTING)

    add_executable(${PROJECT_NAME}
        cpp_nav_filt_unit_tests.cpp
    )

    target_link_libraries(${PROJECT_NAME} PUBLIC
        CppNavFiltSrc
    )

    target_link_libraries(${PROJECT_NAME} PRIVATE
        GTest::gtest_main
    )

    include(GoogleTest)
    gtest_discover_tests(${PROJECT_NAME})

endif()