# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build

# Include any dependencies generated for this target.
include unsupported/test/CMakeFiles/matrix_square_root_5.dir/depend.make

# Include the progress variables for this target.
include unsupported/test/CMakeFiles/matrix_square_root_5.dir/progress.make

# Include the compile flags for this target's objects.
include unsupported/test/CMakeFiles/matrix_square_root_5.dir/flags.make

unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o: unsupported/test/CMakeFiles/matrix_square_root_5.dir/flags.make
unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o: ../unsupported/test/matrix_square_root.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o -c /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/test/matrix_square_root.cpp

unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.i"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/test/matrix_square_root.cpp > CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.i

unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.s"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/test/matrix_square_root.cpp -o CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.s

# Object files for target matrix_square_root_5
matrix_square_root_5_OBJECTS = \
"CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o"

# External object files for target matrix_square_root_5
matrix_square_root_5_EXTERNAL_OBJECTS =

unsupported/test/matrix_square_root_5: unsupported/test/CMakeFiles/matrix_square_root_5.dir/matrix_square_root.cpp.o
unsupported/test/matrix_square_root_5: unsupported/test/CMakeFiles/matrix_square_root_5.dir/build.make
unsupported/test/matrix_square_root_5: unsupported/test/CMakeFiles/matrix_square_root_5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable matrix_square_root_5"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrix_square_root_5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
unsupported/test/CMakeFiles/matrix_square_root_5.dir/build: unsupported/test/matrix_square_root_5

.PHONY : unsupported/test/CMakeFiles/matrix_square_root_5.dir/build

unsupported/test/CMakeFiles/matrix_square_root_5.dir/clean:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test && $(CMAKE_COMMAND) -P CMakeFiles/matrix_square_root_5.dir/cmake_clean.cmake
.PHONY : unsupported/test/CMakeFiles/matrix_square_root_5.dir/clean

unsupported/test/CMakeFiles/matrix_square_root_5.dir/depend:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731 /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/unsupported/test /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/unsupported/test/CMakeFiles/matrix_square_root_5.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unsupported/test/CMakeFiles/matrix_square_root_5.dir/depend

