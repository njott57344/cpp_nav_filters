cmake_minimum_required(VERSION 3.16.3)

project(CppNavFilt)

include_directories(
    include
)

add_library(${PROJECT_NAME}
    src/common.cpp
    src/gps_least_squares.cpp
    src/frame_conversions.cpp
)

add_executable(SvPvtTest
    src/main.cpp
)

target_link_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/../build
)

target_link_libraries(SvPvtTest
    ${PROJECT_NAME}
)