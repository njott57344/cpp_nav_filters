cmake_minimum_required(VERSION 3.16.3)

project(CppNavFilt)

add_library(GpsLocalization
    src/common.cpp 
    src/gps_least_squares.cpp
)

include_directories(
    include
)

add_library(${PROJECT_NAME}
    src/common.cpp
    src/gps_least_squares.cpp
)

add_executable(SvPvtTest
    src/main.cpp
)

target_link_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/../build
)

target_link_libraries(SvPvtTest
    GpsLocalization
)