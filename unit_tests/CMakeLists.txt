cmake_minimum_required(VERSION 3.16.3)

project(CppNavFiltUnitTests)
set(CMAKE_CXX_STANDARD 17)

set(Test1 UnitTestsGpsLeastSquares)
set(Test2 UnitTestsLooselyCoupledINS)

add_executable(Test1
    gps_least_squares_tests.cpp
)

target_link_libraries(Test1 PUBLIC
    CppNavFiltSrc
)