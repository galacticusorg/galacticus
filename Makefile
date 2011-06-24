# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Preprocessor:
PREPROCESSOR = cpp

# Fortran compiler:
F03COMPILER = gfortran

# C compiler:
CCOMPILER = gcc

# C++ compiler:
CPPCOMPILER = g++

# Linker for Condor standard universe executables. Uncomment the second line to link for submission to the Condor standard universe.
CONDORLINKER = 
#CONDORLINKER = condor_compile

# Module type (used for checking if module interfaces have changed):
MODULETYPE = GCC-f95-on-LINUX

# Fortran compiler flags:
F03FLAGS = -ffree-line-length-none -frecursive -J./work/build/ -I./work/build/ ${GALACTICUS_FLAGS} -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules
# Error checking flags
F03FLAGS += -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow
# Add bounds checking.
#F03FLAGS += -fbounds-check
# Add profiling.
F03FLAGS += -g
# A copy of the flags prior to any optimizations.
F03FLAGS_NOOPT := $(F03FLAGS)
# Optimization flags.
#F03FLAGS += -O3 -ffinite-math-only -fno-math-errno -march=native
# For OpenMP compilation.
F03FLAGS += -fopenmp

# C compiler flags:
CFLAGS = -I./source/ -I./work/build/ ${GALACTICUS_FLAGS}
CFLAGS += -g

# C++ compiler flags:
CPPFLAGS = -I./source/ -I./work/build/ ${GALACTICUS_FLAGS}
CPPFLAGS += -g

# Libraries:
LIBS = -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys -lfgsl_gfortran -lgsl -lgslcblas -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lm -lz -lstdc++

# List of additional Makefiles which contain dependency information
MAKE_DEPS = ./work/build/Makefile_Module_Deps ./work/build/Makefile_Use_Deps ./work/build/Makefile_Include_Deps

# General suffix rules: i.e. rules for making a file of one suffix from files of another suffix.

# Object (*.o) files are built by compiling Fortran 90 (*.F90) source files. Ensure that any modules they make are "touch"ed so
# that their modification time is after that of the corresponding object file (to avoid circular remake problems).
vpath %.F90 source
./work/build/%.o : %.F90 ./work/build/%.m ./work/build/%.d Makefile
	@mlist=`cat ./work/build/$*.m` ; \
	for mod in $$mlist ; \
	do \
	 if [ -f $$mod ] ; then mv $$mod $$mod~; fi \
	done
	$(F03COMPILER) -c $< -o ./work/build/$*.o $(F03FLAGS)
	@mlist=`cat ./work/build/$*.m` ; \
	for mod in $$mlist ; \
	do \
	 if [ -f $$mod~ ] ; then \
	  if perl -w ./scripts/build/Compare_Module_Files.pl -compiler $(MODULETYPE) $$mod $$mod~ ; then \
	   mv $$mod~ $$mod ; \
	  else \
	   rm $$mod~ ; \
	  fi \
	 fi \
	done

# Object (*.o) can also be built from C source files.
vpath %.c source
./work/build/%.o : %.c ./work/build/%.d Makefile
	$(CCOMPILER) -c $< -o ./work/build/$*.o $(CFLAGS)

# Object (*.o) can also be built from C++ source files.
vpath %.cpp source
./work/build/%.o : %.cpp ./work/build/%.d Makefile
	$(CPPCOMPILER) -c $< -o ./work/build/$*.o $(CPPFLAGS)

# Special rules required for building some sources (unfortunate, but necessary....)
# bivar.F90 doesn't like to be compiled with any optimization:
./work/build/Bivar/bivar.o : ./source/Bivar/bivar.F90 Makefile
	$(F03COMPILER) -c $< -o ./work/build/Bivar/bivar.o $(F03FLAGS_NOOPT)

# Rule for running *.Inc files through the preprocessor.
./work/build/%.inc : ./work/build/%.Inc Makefile
	$(PREPROCESSOR) -C $< -o ./work/build/$*.tmp
	mv -f ./work/build/$*.tmp ./work/build/$*.inc

# Dependency files (*.d) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.d : ./source/%.F90
	@echo ./work/build/$*.o > ./work/build/$*.d
./work/build/%.d : ./source/%.c
	@echo ./work/build/$*.o > ./work/build/$*.d
./work/build/%.d : ./source/%.cpp
	@echo ./work/build/$*.o > ./work/build/$*.d

# Module list files are created empty by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.m : ./source/%.F90
	@touch ./work/build/$*.m

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the associated dependency (*.d)
# file.
%.exe : ./work/build/%.o ./work/build/%.d `cat ./work/build/$*.d` $(MAKE_DEPS)
	 $(CONDORLINKER) $(F03COMPILER) `cat $*.d` -o $*.exe $(F03FLAGS) $(LIBS)
	 ./scripts/build/Find_Executable_Size.pl $*.exe $*.size

# Ensure that we don't delete object files which make considers to be intermediate
.PRECIOUS: %.o %.d %.dd %.m %.make

# Include depenencies on "include" files.
-include ./work/build/Makefile_Include_Deps 

# Include module dependencies.
-include ./work/build/Makefile_Module_Deps

# Include module use dependencies.
-include ./work/build/Makefile_Use_Deps

# Include rules to build include files generated from directives.
-include ./work/build/Makefile_Directives

# Rules for memory management routines.
./work/build/Allocatable_Arrays.xml: ./scripts/build/Find_Allocatable_Arrays.pl source/*.[fF]90 $(wildcard source/*.Inc)
	./scripts/build/Find_Allocatable_Arrays.pl `pwd`

./work/build/utility.memory_management.precontain.inc: ./scripts/build/Make_Memory_Usage_Routines.pl ./work/build/Allocatable_Arrays.xml
	./scripts/build/Make_Memory_Usage_Routines.pl

./work/build/utility.memory_management.postcontain.inc:
	@touch ./work/build/utility.memory_management.postcontain.inc

# Rules for version routines.
./work/build/galacticus.output.version.revision.inc: $(wildcard .bzr/branch/*)
	@if [ -d .bzr ] ; then awk '{print "integer, parameter :: bazaarRevision="$$1}' .bzr/branch/last-revision > ./work/build/galacticus.output.version.revision.inc; else echo "integer, parameter :: bazaarRevision=-1" > ./work/build/galacticus.output.version.revision.inc; fi

# Rules for cleaning up.
clean: tidy
	rm -f *.exe

tidy:
	rm -rf work/build/*

# Rule for making all executables.
-include ./work/build/Makefile_All_Execs
all: deps $(all_exes)

# Rules for building dependency Makefiles.
./work/build/Makefile_Module_Deps: ./scripts/build/Find_Module_Dependencies.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp $(wildcard source/*.Inc)
	./scripts/build/Find_Module_Dependencies.pl `pwd`

./work/build/Makefile_Use_Deps: ./scripts/build/Find_Use_Dependencies.pl ./work/build/Makefile_Include_Deps source/*.[fF]90 source/*.h source/*.c source/*.cpp $(wildcard source/*.Inc)
	./scripts/build/Find_Use_Dependencies.pl `pwd` $(MAKE)

./work/build/Makefile_Directives: ./scripts/build/Code_Directive_Parser.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	./scripts/build/Code_Directive_Parser.pl `pwd`

./work/build/Makefile_Include_Deps: ./scripts/build/Find_Include_Dependencies.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	./scripts/build/Find_Include_Dependencies.pl `pwd`

./work/build/Makefile_All_Execs: ./scripts/build/Find_Programs.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	./scripts/build/Find_Programs.pl `pwd`

deps: $(MAKE_DEPS) ./work/build/Makefile_All_Execs
