# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

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
CMAKE_SOURCE_DIR = /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build

# Include any dependencies generated for this target.
include CMakeFiles/exampleB3a.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/exampleB3a.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/exampleB3a.dir/flags.make

CMakeFiles/exampleB3a.dir/exampleB3a.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/exampleB3a.cc.o: ../exampleB3a.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/exampleB3a.dir/exampleB3a.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/exampleB3a.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/exampleB3a.cc

CMakeFiles/exampleB3a.dir/exampleB3a.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/exampleB3a.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/exampleB3a.cc > CMakeFiles/exampleB3a.dir/exampleB3a.cc.i

CMakeFiles/exampleB3a.dir/exampleB3a.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/exampleB3a.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/exampleB3a.cc -o CMakeFiles/exampleB3a.dir/exampleB3a.cc.s

CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.requires

CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.provides: CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.provides

CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.provides.build: CMakeFiles/exampleB3a.dir/exampleB3a.cc.o


CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o: ../src/B3DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3DetectorConstruction.cc

CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3DetectorConstruction.cc > CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.i

CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3DetectorConstruction.cc -o CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.s

CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o


CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o: ../src/B3Hit.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3Hit.cc

CMakeFiles/exampleB3a.dir/src/B3Hit.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3Hit.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3Hit.cc > CMakeFiles/exampleB3a.dir/src/B3Hit.cc.i

CMakeFiles/exampleB3a.dir/src/B3Hit.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3Hit.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3Hit.cc -o CMakeFiles/exampleB3a.dir/src/B3Hit.cc.s

CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o


CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o: ../src/B3PhysicsList.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsList.cc

CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsList.cc > CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.i

CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsList.cc -o CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.s

CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o


CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o: ../src/B3PhysicsListMessenger.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsListMessenger.cc

CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsListMessenger.cc > CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.i

CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PhysicsListMessenger.cc -o CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.s

CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o


CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o: ../src/B3PrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PrimaryGeneratorAction.cc

CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PrimaryGeneratorAction.cc > CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.i

CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3PrimaryGeneratorAction.cc -o CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.s

CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o


CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o: ../src/B3StackingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3StackingAction.cc

CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3StackingAction.cc > CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.i

CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3StackingAction.cc -o CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.s

CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o


CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o: ../src/B3SteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3SteppingAction.cc

CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3SteppingAction.cc > CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.i

CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3SteppingAction.cc -o CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.s

CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o


CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o: ../src/B3aActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aActionInitialization.cc

CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aActionInitialization.cc > CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.i

CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aActionInitialization.cc -o CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.s

CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o


CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o: ../src/B3aEventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aEventAction.cc

CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aEventAction.cc > CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.i

CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aEventAction.cc -o CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.s

CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o


CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o: ../src/B3aHistoManager.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aHistoManager.cc

CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aHistoManager.cc > CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.i

CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aHistoManager.cc -o CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.s

CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o


CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o: ../src/B3aRun.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRun.cc

CMakeFiles/exampleB3a.dir/src/B3aRun.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3aRun.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRun.cc > CMakeFiles/exampleB3a.dir/src/B3aRun.cc.i

CMakeFiles/exampleB3a.dir/src/B3aRun.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3aRun.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRun.cc -o CMakeFiles/exampleB3a.dir/src/B3aRun.cc.s

CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o


CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o: CMakeFiles/exampleB3a.dir/flags.make
CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o: ../src/B3aRunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o"
	/usr/bin/g++-7  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o -c /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRunAction.cc

CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.i"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRunAction.cc > CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.i

CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.s"
	/usr/bin/g++-7 $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/src/B3aRunAction.cc -o CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.s

CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.requires:

.PHONY : CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.requires

CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.provides: CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/exampleB3a.dir/build.make CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.provides.build
.PHONY : CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.provides

CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.provides.build: CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o


# Object files for target exampleB3a
exampleB3a_OBJECTS = \
"CMakeFiles/exampleB3a.dir/exampleB3a.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o" \
"CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o"

# External object files for target exampleB3a
exampleB3a_EXTERNAL_OBJECTS =

exampleB3a: CMakeFiles/exampleB3a.dir/exampleB3a.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o
exampleB3a: CMakeFiles/exampleB3a.dir/build.make
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4Tree.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4GMocren.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4visHepRep.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4RayTracer.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4VRML.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4OpenGL.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4gl2ps.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4interfaces.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4persistency.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4error_propagation.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4readout.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4physicslists.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4parmodels.so
exampleB3a: /usr/local/lib/libxerces-c.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4FR.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4vis_management.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4modeling.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libXmu.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libXext.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libXt.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libSM.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libICE.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libX11.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libGLU.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libGL.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.9.5
exampleB3a: /usr/lib/x86_64-linux-gnu/libQt5PrintSupport.so.5.9.5
exampleB3a: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.9.5
exampleB3a: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.9.5
exampleB3a: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.9.5
exampleB3a: /usr/local/lib/libxerces-c.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4run.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4event.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4tracking.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4processes.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4analysis.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4zlib.so
exampleB3a: /usr/lib/x86_64-linux-gnu/libexpat.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4digits_hits.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4track.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4particles.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4geometry.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4materials.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4graphics_reps.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4intercoms.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4global.so
exampleB3a: /home/mitt-unix/Bureau/Desktop/Geant4/G4-install/lib/libG4clhep.so
exampleB3a: CMakeFiles/exampleB3a.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX executable exampleB3a"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/exampleB3a.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/exampleB3a.dir/build: exampleB3a

.PHONY : CMakeFiles/exampleB3a.dir/build

CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/exampleB3a.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3DetectorConstruction.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3Hit.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3PhysicsList.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3PhysicsListMessenger.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3PrimaryGeneratorAction.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3StackingAction.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3SteppingAction.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3aActionInitialization.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3aEventAction.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3aHistoManager.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3aRun.cc.o.requires
CMakeFiles/exampleB3a.dir/requires: CMakeFiles/exampleB3a.dir/src/B3aRunAction.cc.o.requires

.PHONY : CMakeFiles/exampleB3a.dir/requires

CMakeFiles/exampleB3a.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/exampleB3a.dir/cmake_clean.cmake
.PHONY : CMakeFiles/exampleB3a.dir/clean

CMakeFiles/exampleB3a.dir/depend:
	cd /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build /home/mitt-unix/Bureau/Desktop/Geant4/PlasticPET/src/B3a/build/CMakeFiles/exampleB3a.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/exampleB3a.dir/depend

