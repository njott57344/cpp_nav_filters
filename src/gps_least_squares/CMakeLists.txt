cmake_minimum_required(VERSION 3.16.3)

project(GpsLeastSquares)

add_library(GpsLeastSquaresLib src/gps_least_squares.cpp)

include_directories(
    include
)

add_library(${PROJECT_NAME}
    src/gps_least_squares.cpp
)

target_link_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/../build
)