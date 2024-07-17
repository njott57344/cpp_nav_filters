cmake_minimum_required(VERSION 3.16.3)

project(CppNavFiltSrc)
set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

set(Headers
    include/cpp_nav_filt_lib.h
    include/gps_least_squares.h
    include/sv_manager.h
    include/ins.h
)

set(Sources
    src/cpp_nav_filt_lib.cpp
    src/gps_least_squares.cpp
    src/sv_manager.cpp
    src/ins.cpp
)

add_library(${PROJECT_NAME} STATIC ${Sources} ${Headers})

target_include_directories(
    ${PROJECT_NAME} PUBLIC
    include # cpp nav filt headers
    ${PYTHON_INCLUDE_DIRS}
)