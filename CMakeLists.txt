cmake_minimum_required(VERSION 3.16.3)

project(CppNavFilt)
set(CMAKE_CXX_STANDARD 17)

# CMake Stuff for Plotting
include(FetchContent)

FetchContent_Declare(matplotplusplus
    GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus
    GIT_TAG origin/master
)

FetchContent_GetProperties(matplotplusplus)

if(NOT matplotplusplus_POPULATED)
    FetchContent_Populate(matplotplusplus)
    add_subdirectory(${matplotplusplus_SOURCE_DIR} ${matplotplusplus_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Cmake Stuff for Matplot (c++ wrapper for pyton matplot)
find_package(PythonLibs REQUIRED)

include_directories(
    include # cpp nav filt headers
    ${PYTHON_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
    src/sv_manager.cpp
    src/gps_least_squares.cpp
    src/lc_ins.cpp
    src/cpp_nav_filt_lib.cpp
)

add_executable(SvPvtTest
    src/main.cpp
)

target_link_directories(${PROJECT_NAME}
    PRIVATE ${CMAKE_SOURCE_DIR}/../build
)

target_link_libraries(SvPvtTest PUBLIC matplot
    ${PROJECT_NAME}
    ${PYTHON_LIBRARIES}
)