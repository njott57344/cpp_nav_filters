cmake_minimum_required(VERSION 3.16.3)

project(Common)

add_library(CommonLib src/common.cpp)

include_directories(
    include
)

add_library(${PROJECT_NAME}
    src/common.cpp
)

target_link_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/../build
)

target_link_libraries(${PROJECT_NAME}
    CommonLib
)