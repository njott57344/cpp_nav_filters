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
CMAKE_SOURCE_DIR = /home/njo0004/devel/cpp_nav_filters

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/njo0004/devel/cpp_nav_filters/build

# Include any dependencies generated for this target.
include CMakeFiles/CppNavFilt.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CppNavFilt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CppNavFilt.dir/flags.make

CMakeFiles/CppNavFilt.dir/src/common.cpp.o: CMakeFiles/CppNavFilt.dir/flags.make
CMakeFiles/CppNavFilt.dir/src/common.cpp.o: ../src/common.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/njo0004/devel/cpp_nav_filters/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/CppNavFilt.dir/src/common.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CppNavFilt.dir/src/common.cpp.o -c /home/njo0004/devel/cpp_nav_filters/src/common.cpp

CMakeFiles/CppNavFilt.dir/src/common.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CppNavFilt.dir/src/common.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/njo0004/devel/cpp_nav_filters/src/common.cpp > CMakeFiles/CppNavFilt.dir/src/common.cpp.i

CMakeFiles/CppNavFilt.dir/src/common.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CppNavFilt.dir/src/common.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/njo0004/devel/cpp_nav_filters/src/common.cpp -o CMakeFiles/CppNavFilt.dir/src/common.cpp.s

CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o: CMakeFiles/CppNavFilt.dir/flags.make
CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o: ../src/gps_least_squares.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/njo0004/devel/cpp_nav_filters/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o -c /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares.cpp

CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares.cpp > CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.i

CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/njo0004/devel/cpp_nav_filters/src/gps_least_squares.cpp -o CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.s

CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o: CMakeFiles/CppNavFilt.dir/flags.make
CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o: ../src/frame_conversions.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/njo0004/devel/cpp_nav_filters/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o -c /home/njo0004/devel/cpp_nav_filters/src/frame_conversions.cpp

CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/njo0004/devel/cpp_nav_filters/src/frame_conversions.cpp > CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.i

CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/njo0004/devel/cpp_nav_filters/src/frame_conversions.cpp -o CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.s

# Object files for target CppNavFilt
CppNavFilt_OBJECTS = \
"CMakeFiles/CppNavFilt.dir/src/common.cpp.o" \
"CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o" \
"CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o"

# External object files for target CppNavFilt
CppNavFilt_EXTERNAL_OBJECTS =

libCppNavFilt.a: CMakeFiles/CppNavFilt.dir/src/common.cpp.o
libCppNavFilt.a: CMakeFiles/CppNavFilt.dir/src/gps_least_squares.cpp.o
libCppNavFilt.a: CMakeFiles/CppNavFilt.dir/src/frame_conversions.cpp.o
libCppNavFilt.a: CMakeFiles/CppNavFilt.dir/build.make
libCppNavFilt.a: CMakeFiles/CppNavFilt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/njo0004/devel/cpp_nav_filters/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX static library libCppNavFilt.a"
	$(CMAKE_COMMAND) -P CMakeFiles/CppNavFilt.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CppNavFilt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CppNavFilt.dir/build: libCppNavFilt.a

.PHONY : CMakeFiles/CppNavFilt.dir/build

CMakeFiles/CppNavFilt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CppNavFilt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CppNavFilt.dir/clean

CMakeFiles/CppNavFilt.dir/depend:
	cd /home/njo0004/devel/cpp_nav_filters/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/njo0004/devel/cpp_nav_filters /home/njo0004/devel/cpp_nav_filters /home/njo0004/devel/cpp_nav_filters/build /home/njo0004/devel/cpp_nav_filters/build /home/njo0004/devel/cpp_nav_filters/build/CMakeFiles/CppNavFilt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CppNavFilt.dir/depend

