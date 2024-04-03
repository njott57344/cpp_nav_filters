cmake_minimum_required(VERSION 3.16.3)

project(CppNavFilt)
set(CMAKE_CXX_STANDARD 17)

# Update Git Submodules
include(cmake/UpdateSubmodules.cmake)

# Submodules
add_subdirectory(cpp_nav_filt_src)
add_subdirectory(demo)
add_subdirectory(external)
add_subdirectory(unit_tests)
