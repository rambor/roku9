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
include doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/depend.make

# Include the progress variables for this target.
include doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/progress.make

# Include the compile flags for this target's objects.
include doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/flags.make

doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o: doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/flags.make
doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o: ../doc/examples/Tutorial_BlockOperations_block_assignment.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o -c /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/doc/examples/Tutorial_BlockOperations_block_assignment.cpp

doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.i"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/doc/examples/Tutorial_BlockOperations_block_assignment.cpp > CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.i

doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.s"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/doc/examples/Tutorial_BlockOperations_block_assignment.cpp -o CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.s

# Object files for target Tutorial_BlockOperations_block_assignment
Tutorial_BlockOperations_block_assignment_OBJECTS = \
"CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o"

# External object files for target Tutorial_BlockOperations_block_assignment
Tutorial_BlockOperations_block_assignment_EXTERNAL_OBJECTS =

doc/examples/Tutorial_BlockOperations_block_assignment: doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/Tutorial_BlockOperations_block_assignment.cpp.o
doc/examples/Tutorial_BlockOperations_block_assignment: doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/build.make
doc/examples/Tutorial_BlockOperations_block_assignment: doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable Tutorial_BlockOperations_block_assignment"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/link.txt --verbose=$(VERBOSE)
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && ./Tutorial_BlockOperations_block_assignment >/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples/Tutorial_BlockOperations_block_assignment.out

# Rule to build all files generated by this target.
doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/build: doc/examples/Tutorial_BlockOperations_block_assignment

.PHONY : doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/build

doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/clean:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples && $(CMAKE_COMMAND) -P CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/cmake_clean.cmake
.PHONY : doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/clean

doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/depend:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731 /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/doc/examples /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/examples/CMakeFiles/Tutorial_BlockOperations_block_assignment.dir/depend

