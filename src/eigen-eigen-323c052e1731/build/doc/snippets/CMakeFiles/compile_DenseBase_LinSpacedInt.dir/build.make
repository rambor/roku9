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
include doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/depend.make

# Include the progress variables for this target.
include doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/progress.make

# Include the compile flags for this target's objects.
include doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/flags.make

doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o: doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/flags.make
doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o: doc/snippets/compile_DenseBase_LinSpacedInt.cpp
doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o: ../doc/snippets/DenseBase_LinSpacedInt.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o -c /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets/compile_DenseBase_LinSpacedInt.cpp

doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.i"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets/compile_DenseBase_LinSpacedInt.cpp > CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.i

doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.s"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && /usr/local/Cellar/gcc/8.3.0_2/bin/g++-8 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets/compile_DenseBase_LinSpacedInt.cpp -o CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.s

# Object files for target compile_DenseBase_LinSpacedInt
compile_DenseBase_LinSpacedInt_OBJECTS = \
"CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o"

# External object files for target compile_DenseBase_LinSpacedInt
compile_DenseBase_LinSpacedInt_EXTERNAL_OBJECTS =

doc/snippets/compile_DenseBase_LinSpacedInt: doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/compile_DenseBase_LinSpacedInt.cpp.o
doc/snippets/compile_DenseBase_LinSpacedInt: doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/build.make
doc/snippets/compile_DenseBase_LinSpacedInt: doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable compile_DenseBase_LinSpacedInt"
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/compile_DenseBase_LinSpacedInt.dir/link.txt --verbose=$(VERBOSE)
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && ./compile_DenseBase_LinSpacedInt >/Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets/DenseBase_LinSpacedInt.out

# Rule to build all files generated by this target.
doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/build: doc/snippets/compile_DenseBase_LinSpacedInt

.PHONY : doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/build

doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/clean:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets && $(CMAKE_COMMAND) -P CMakeFiles/compile_DenseBase_LinSpacedInt.dir/cmake_clean.cmake
.PHONY : doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/clean

doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/depend:
	cd /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731 /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/doc/snippets /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets /Users/xos81802/CLionProjects/aligner/src/eigen-eigen-323c052e1731/build/doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : doc/snippets/CMakeFiles/compile_DenseBase_LinSpacedInt.dir/depend

