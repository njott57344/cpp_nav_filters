cmake_minimum_required(VERSION 3.16.3)

project(CppNavFiltDemo)
set(CMAKE_CXX_STANDARD 17)

set(DemoGpsLeastSquares GpsLeastSquares)

add_executable(DemoGpsLeastSquares
    gps_least_squares_demo.cpp
)

target_link_libraries(DemoGpsLeastSquares PUBLIC
    CppNavFiltSrc
)
