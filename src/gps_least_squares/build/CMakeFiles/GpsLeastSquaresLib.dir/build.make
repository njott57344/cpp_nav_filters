# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build

# Include any dependencies generated for this target.
include CMakeFiles/GpsLeastSquaresLib.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/GpsLeastSquaresLib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GpsLeastSquaresLib.dir/flags.make

CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o: CMakeFiles/GpsLeastSquaresLib.dir/flags.make
CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o: ../src/gps_least_squares.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o -c /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/src/gps_least_squares.cpp

CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/src/gps_least_squares.cpp > CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.i

CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/src/gps_least_squares.cpp -o CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.s

# Object files for target GpsLeastSquaresLib
GpsLeastSquaresLib_OBJECTS = \
"CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o"

# External object files for target GpsLeastSquaresLib
GpsLeastSquaresLib_EXTERNAL_OBJECTS =

libGpsLeastSquaresLib.a: CMakeFiles/GpsLeastSquaresLib.dir/src/gps_least_squares.cpp.o
libGpsLeastSquaresLib.a: CMakeFiles/GpsLeastSquaresLib.dir/build.make
libGpsLeastSquaresLib.a: CMakeFiles/GpsLeastSquaresLib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libGpsLeastSquaresLib.a"
	$(CMAKE_COMMAND) -P CMakeFiles/GpsLeastSquaresLib.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GpsLeastSquaresLib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GpsLeastSquaresLib.dir/build: libGpsLeastSquaresLib.a

.PHONY : CMakeFiles/GpsLeastSquaresLib.dir/build

CMakeFiles/GpsLeastSquaresLib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GpsLeastSquaresLib.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GpsLeastSquaresLib.dir/clean

CMakeFiles/GpsLeastSquaresLib.dir/depend:
	cd /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares/build/CMakeFiles/GpsLeastSquaresLib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GpsLeastSquaresLib.dir/depend

