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
CMAKE_SOURCE_DIR = /home/njo0004/devel/cpp_nav_filters/src/common

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/njo0004/devel/cpp_nav_filters/src/common/build

# Include any dependencies generated for this target.
include CMakeFiles/Common.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Common.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Common.dir/flags.make

CMakeFiles/Common.dir/src/common.cpp.o: CMakeFiles/Common.dir/flags.make
CMakeFiles/Common.dir/src/common.cpp.o: ../src/common.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/njo0004/devel/cpp_nav_filters/src/common/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Common.dir/src/common.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Common.dir/src/common.cpp.o -c /home/njo0004/devel/cpp_nav_filters/src/common/src/common.cpp

CMakeFiles/Common.dir/src/common.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Common.dir/src/common.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/njo0004/devel/cpp_nav_filters/src/common/src/common.cpp > CMakeFiles/Common.dir/src/common.cpp.i

CMakeFiles/Common.dir/src/common.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Common.dir/src/common.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/njo0004/devel/cpp_nav_filters/src/common/src/common.cpp -o CMakeFiles/Common.dir/src/common.cpp.s

# Object files for target Common
Common_OBJECTS = \
"CMakeFiles/Common.dir/src/common.cpp.o"

# External object files for target Common
Common_EXTERNAL_OBJECTS =

libCommon.a: CMakeFiles/Common.dir/src/common.cpp.o
libCommon.a: CMakeFiles/Common.dir/build.make
libCommon.a: CMakeFiles/Common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/njo0004/devel/cpp_nav_filters/src/common/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libCommon.a"
	$(CMAKE_COMMAND) -P CMakeFiles/Common.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Common.dir/build: libCommon.a

.PHONY : CMakeFiles/Common.dir/build

CMakeFiles/Common.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Common.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Common.dir/clean

CMakeFiles/Common.dir/depend:
	cd /home/njo0004/devel/cpp_nav_filters/src/common/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/njo0004/devel/cpp_nav_filters/src/common /home/njo0004/devel/cpp_nav_filters/src/common /home/njo0004/devel/cpp_nav_filters/src/common/build /home/njo0004/devel/cpp_nav_filters/src/common/build /home/njo0004/devel/cpp_nav_filters/src/common/build/CMakeFiles/Common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Common.dir/depend
