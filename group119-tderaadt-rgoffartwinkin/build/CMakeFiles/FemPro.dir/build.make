# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build

# Include any dependencies generated for this target.
include CMakeFiles/FemPro.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/FemPro.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/FemPro.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/FemPro.dir/flags.make

CMakeFiles/FemPro.dir/src/fem.c.o: CMakeFiles/FemPro.dir/flags.make
CMakeFiles/FemPro.dir/src/fem.c.o: ../src/fem.c
CMakeFiles/FemPro.dir/src/fem.c.o: CMakeFiles/FemPro.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/FemPro.dir/src/fem.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/FemPro.dir/src/fem.c.o -MF CMakeFiles/FemPro.dir/src/fem.c.o.d -o CMakeFiles/FemPro.dir/src/fem.c.o -c /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/fem.c

CMakeFiles/FemPro.dir/src/fem.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/FemPro.dir/src/fem.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/fem.c > CMakeFiles/FemPro.dir/src/fem.c.i

CMakeFiles/FemPro.dir/src/fem.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/FemPro.dir/src/fem.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/fem.c -o CMakeFiles/FemPro.dir/src/fem.c.s

CMakeFiles/FemPro.dir/src/main.c.o: CMakeFiles/FemPro.dir/flags.make
CMakeFiles/FemPro.dir/src/main.c.o: ../src/main.c
CMakeFiles/FemPro.dir/src/main.c.o: CMakeFiles/FemPro.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object CMakeFiles/FemPro.dir/src/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT CMakeFiles/FemPro.dir/src/main.c.o -MF CMakeFiles/FemPro.dir/src/main.c.o.d -o CMakeFiles/FemPro.dir/src/main.c.o -c /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/main.c

CMakeFiles/FemPro.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/FemPro.dir/src/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/main.c > CMakeFiles/FemPro.dir/src/main.c.i

CMakeFiles/FemPro.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/FemPro.dir/src/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/src/main.c -o CMakeFiles/FemPro.dir/src/main.c.s

# Object files for target FemPro
FemPro_OBJECTS = \
"CMakeFiles/FemPro.dir/src/fem.c.o" \
"CMakeFiles/FemPro.dir/src/main.c.o"

# External object files for target FemPro
FemPro_EXTERNAL_OBJECTS =

FemPro: CMakeFiles/FemPro.dir/src/fem.c.o
FemPro: CMakeFiles/FemPro.dir/src/main.c.o
FemPro: CMakeFiles/FemPro.dir/build.make
FemPro: CMakeFiles/FemPro.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking C executable FemPro"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/FemPro.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/FemPro.dir/build: FemPro
.PHONY : CMakeFiles/FemPro.dir/build

CMakeFiles/FemPro.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/FemPro.dir/cmake_clean.cmake
.PHONY : CMakeFiles/FemPro.dir/clean

CMakeFiles/FemPro.dir/depend:
	cd /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build /home/pepsiduck/Bureau/FEM_Project/group119-tderaadt-rgoffartwinkin/build/CMakeFiles/FemPro.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/FemPro.dir/depend

