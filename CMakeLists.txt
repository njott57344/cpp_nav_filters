cmake_minimum_required(VERSION 3.16.3)

project(CppNavFilt)
set(CMAKE_CXX_STANDARD 17)

# CMake Stuff for Sci Plot
include(FetchContent)

FetchContent_Declare(sciplot_content
    GIT_REPOSITORY https://github.com/sciplot/sciplot.git
    GIT_TAG master
)

FetchContent_GetProperties(sciplot_content)

if(NOT sciplot_content_POPULATED)
    FetchContent_Populate(sciplot_content)
endif()

include_directories(
    include # cpp nav filt headers
    ${sciplot_content_SOURCE_DIR}
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