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
include test/CMakeFiles/array_replicate_4.dir/depend.make

# Include the progress variables for this target.
include test/CMakeFiles/array_replicate_4.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/array_replicate_4.dir/flags.make

test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o: test/CMakeFiles/array_replicate_4.dir/flags.make
test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o: ../test/array_replicate.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o -c /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/test/array_replicate.cpp

test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/array_replicate_4.dir/array_replicate.cpp.i"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/test/array_replicate.cpp > CMakeFiles/array_replicate_4.dir/array_replicate.cpp.i

test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/array_replicate_4.dir/array_replicate.cpp.s"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/test/array_replicate.cpp -o CMakeFiles/array_replicate_4.dir/array_replicate.cpp.s

# Object files for target array_replicate_4
array_replicate_4_OBJECTS = \
"CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o"

# External object files for target array_replicate_4
array_replicate_4_EXTERNAL_OBJECTS =

test/array_replicate_4: test/CMakeFiles/array_replicate_4.dir/array_replicate.cpp.o
test/array_replicate_4: test/CMakeFiles/array_replicate_4.dir/build.make
test/array_replicate_4: test/CMakeFiles/array_replicate_4.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable array_replicate_4"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/array_replicate_4.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/array_replicate_4.dir/build: test/array_replicate_4

.PHONY : test/CMakeFiles/array_replicate_4.dir/build

test/CMakeFiles/array_replicate_4.dir/clean:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test && $(CMAKE_COMMAND) -P CMakeFiles/array_replicate_4.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/array_replicate_4.dir/clean

test/CMakeFiles/array_replicate_4.dir/depend:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731 /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/test /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/test/CMakeFiles/array_replicate_4.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/array_replicate_4.dir/depend

