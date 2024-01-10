# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.16.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(sciplot_content-populate NONE)

include(ExternalProject)
ExternalProject_Add(sciplot_content-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/sciplot/sciplot.git" "GIT_TAG" "master"
                    SOURCE_DIR          "/home/njo0004/devel/cpp_nav_filters/build/_deps/sciplot_content-src"
                    BINARY_DIR          "/home/njo0004/devel/cpp_nav_filters/build/_deps/sciplot_content-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
)
