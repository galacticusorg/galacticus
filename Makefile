# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Preprocessor:
PREPROCESSOR ?= cpp

# Fortran compiler:
FCCOMPILER ?= mpif90

# C compiler:
CCOMPILER ?= mpicc

# C++ compiler:
CPPCOMPILER ?= g++

# Linker for Condor standard universe executables. Uncomment the second line to link for submission to the Condor standard universe.
CONDORLINKER = 
#CONDORLINKER = condor_compile

# Module type (used for checking if module interfaces have changed):
MODULETYPE ?= GCC-f95-on-LINUX

# Fortran compiler flags:
FCFLAGS += -ffree-line-length-none -frecursive -J./work/build/ -I./work/build/ ${GALACTICUS_FCFLAGS} -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules -fintrinsic-modules-path /usr/lib64/openmpi/lib
# Fortran77 compiler flags:
F77FLAGS = -g
# Error checking flags
FCFLAGS += -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core
# Add bounds checking.
#FCFLAGS += -fbounds-check
# Add profiling.
FCFLAGS += -g
# A copy of the flags prior to any optimizations.
FCFLAGS_NOOPT := $(FCFLAGS)
# Optimization flags.
FCFLAGS += -O3 -ffinite-math-only -fno-math-errno -march=native
# For OpenMP compilation.
FCFLAGS += -fopenmp

# C compiler flags:
CFLAGS = -I./source/ -I./work/build/ -I/opt/gsl-trunk/include ${GALACTICUS_CFLAGS}
CFLAGS += -g

# C++ compiler flags:
CPPFLAGS += -I./source/ -I./work/build/ ${GALACTICUS_CPPFLAGS}
CPPFLAGS += -g

# List of additional Makefiles which contain dependency information
MAKE_DEPS = ./work/build/Makefile_Module_Deps ./work/build/Makefile_Use_Deps ./work/build/Makefile_Include_Deps

# Get versions of build tools.
FCCOMPILER_VERSION = `$(FCCOMPILER) -v 2>&1`
CCOMPILER_VERSION = `$(CCOMPILER) -v 2>&1`
CPPCOMPILER_VERSION = `$(CPPCOMPILER) -v 2>&1`

# General suffix rules: i.e. rules for making a file of one suffix from files of another suffix.

# Object (*.o) files are built by compiling Fortran 90 (*.F90) source files. Ensure that any modules they make are "touch"ed so
# that their modification time is after that of the corresponding object file (to avoid circular remake problems).
vpath %.F90 source
./work/build/%.o : %.F90 ./work/build/%.m ./work/build/%.d ./work/build/%.fl Makefile
	@mlist=`cat ./work/build/$*.m` ; \
	for mod in $$mlist ; \
	do \
	 if [ -f $$mod ] ; then mv $$mod $$mod~; fi \
	done
	$(FCCOMPILER) -c $< -o ./work/build/$*.o $(FCFLAGS)
	@mlist=`cat ./work/build/$*.m` ; \
	for mod in $$mlist ; \
	do \
	 if [ -f $$mod~ ] ; then \
	  file $$mod | grep -q ASCII ; \
	  if [ $$? -eq 0 ]; then \
	   if perl -w ./scripts/build/Compare_Module_Files.pl -compiler $(MODULETYPE) $$mod $$mod~ ; then \
	    mv $$mod~ $$mod ; \
	   else \
	    rm $$mod~ ; \
	   fi \
	 else \
	   gunzip -c $$mod > $$mod.gu ; \
	   gunzip -c $$mod~ > $$mod~.gu ; \
	   if perl -w ./scripts/build/Compare_Module_Files.pl -compiler $(MODULETYPE) $$mod.gu $$mod~.gu ; then \
	    mv $$mod~ $$mod ; \
	   else \
	    rm $$mod~ ; \
	   fi ; \
	   rm $$mod.gu $$mod~.gu ; \
          fi \
	 fi \
	done

# Object (*.o) files are built by compiling C (*.c) source files.
vpath %.c source
./work/build/%.o : %.c ./work/build/%.d ./work/build/%.fl Makefile
	$(CCOMPILER) -c $< -o ./work/build/$*.o $(CFLAGS)

# Object (*.o) can also be built from C++ source files.
vpath %.cpp source
./work/build/%.o : %.cpp ./work/build/%.d ./work/build/%.fl Makefile
	$(CPPCOMPILER) -c $< -o ./work/build/$*.o $(CPPFLAGS)

# Object (*.o) files are built by compiling Fortran (*.f) source files.
%.o : %.f %.d Makefile
	$(FCCOMPILER) -c $< -o $*.o $(F77FLAGS)

# Special rules required for building some sources (unfortunate, but necessary....)
# bivar.F90 doesn't like to be compiled with any optimization:
./work/build/Bivar/bivar.o : ./source/Bivar/bivar.F90 Makefile
	$(FCCOMPILER) -c $< -o ./work/build/Bivar/bivar.o $(FCFLAGS_NOOPT)

# Rule for running *.Inc files through the preprocessor.
./work/build/%.Inc : ./source/%.Inc
	cp -f ./source/$*.Inc ./work/build/$*.Inc
./work/build/%.inc : ./work/build/%.Inc Makefile
	$(PREPROCESSOR) -nostdinc -C $< -o ./work/build/$*.tmp
	mv -f ./work/build/$*.tmp ./work/build/$*.inc

# Dependency files (*.d) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.d : ./source/%.F90
	@echo ./work/build/$*.o > ./work/build/$*.d
./work/build/%.d : ./source/%.c
	@echo ./work/build/$*.o > ./work/build/$*.d
./work/build/%.d : ./source/%.cpp
	@echo ./work/build/$*.o > ./work/build/$*.d
%.d : %.f
	@echo $*.o > $*.d
%.d :
	@mkdir -p `dirname $*.d`
	@touch $*.d


# Library files (*.fl) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps file, but this acts as a fallback rule.
./work/build/%.fl : ./source/%.F90
	@touch ./work/build/$*.fl
./work/build/%.fl : ./source/%.c
	@touch ./work/build/$*.fl
./work/build/%.fl : ./source/%.cpp
	@touch ./work/build/$*.fl

# GraphViz files (*.gv) are created with just a node entry by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.F90.gv : ./source/%.F90
	@echo \"$*.F90\" > ./work/build/$*.F90.gv
./work/build/%.c.gv : ./source/%.c
	@echo \"$*.c\" > ./work/build/$*.c.gv
./work/build/%.cpp.gv : ./source/%.cpp
	@echo \"$*.cpp\" > ./work/build/$*.cpp.gv

# Create a PostScript files showing tree diagrams of source file dependencies.
%.tps : ./work/build/%.gv
	@echo digraph Tree \{ > $*.dot
	@cat ./work/build/$*.gv >> $*.dot
	@echo \} >> $*.dot
	dot -Tps $*.dot -o $*.tps
	@rm $*.dot

# Module list files are created empty by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.m : ./source/%.F90
	@touch ./work/build/$*.m

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the associated dependency (*.d)
# file.
%.exe : ./work/build/%.o ./work/build/%.d `cat ./work/build/$*.d` $(MAKE_DEPS)
	 $(CONDORLINKER) $(FCCOMPILER) `cat $*.d` -o $*.exe $(FCFLAGS) `scripts/build/Library_Dependencies.pl $*.exe $(FCFLAGS)`
	 ./scripts/build/Find_Executable_Size.pl $*.exe $*.size
	 ./scripts/build/Find_Parameter_Dependencies.pl `pwd` $*.exe

# Ensure that we don't delete object files which make considers to be intermediate
.PRECIOUS: %.o %.d %.dd %.m %.make %.Inc

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
./work/build/galacticus.output.version.revision.inc: $(wildcard .hg/branch)
	@if [ -f .hg/branch ] ; then hg tip | awk 'BEGIN {FS=":";r=-1;h=""} {if ((NR == 1 && NF == 3 ) || $$1 == "changeset") {r=$$2;h=$$3}} END {print "integer, parameter :: hgRevision="r"\ncharacter(len=12), parameter :: hgHash=\""h"\""}' > ./work/build/galacticus.output.version.revision.inc; else printf 'integer, parameter :: hgRevision=-1\ncharacter(len=12), parameter :: hgHash=""\n' > ./work/build/galacticus.output.version.revision.inc; fi

# Rules for build information routines.
./work/build/galacticus.output.build.environment.inc:
	@echo PREPROCESSOR=\"$(PREPROCESSOR)\" > ./work/build/galacticus.output.build.environment.inc
	@echo FCCOMPILER=\"$(FCCOMPILER)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CCOMPILER=\"$(CCOMPILER)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CPPCOMPILER=\"$(CPPCOMPILER)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo MODULETYPE=\"$(MODULETYPE)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo FCFLAGS=\"$(FCFLAGS)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo FCFLAGS_NOOPT=\"$(FCFLAGS_NOOPT)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CFLAGS=\"$(CFLAGS)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CPPFLAGS=\"$(CPPFLAGS)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo LIBS=\"$(LIBS)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo FCCOMPILER_VERSION=\"$(FCCOMPILER_VERSION)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CCOMPILER_VERSION=\"$(CCOMPILER_VERSION)\" >> ./work/build/galacticus.output.build.environment.inc
	@echo CPPCOMPILER_VERSION=\"$(CPPCOMPILER_VERSION)\" >> ./work/build/galacticus.output.build.environment.inc

# Rules for unique label function creation.
dfiles := $(patsubst source/%.F90,work/build/%.d,$(wildcard source/*.F90))
mfiles := $(patsubst source/%.F90,work/build/%.m,$(wildcard source/*.F90))
./work/build/utility.input_parameters.unique_labels.inc:
	@touch ./work/build/utility.input_parameters.unique_labels.inc
./work/build/utility.input_parameters.unique_labels.visibilities.inc: $(dfiles) $(mfiles)
	scripts/build/Make_Unique_Label_Functions.pl `pwd`

# Rules for changeset creation.
Galacticus.exe: ./work/build/galacticus.hg.patch ./work/build/galacticus.hg.bundle
./work/build/galacticus.hg.patch:
	hg diff > ./work/build/galacticus.hg.patch | true
./work/build/galacticus.hg.bundle:
	hg bundle -t none ./work/build/galacticus.hg.bundle https://abensonca@bitbucket.org/abensonca/galacticus | true

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
	@mkdir -p work/build
	./scripts/build/Find_Module_Dependencies.pl `pwd`

./work/build/Makefile_Use_Deps: ./scripts/build/Find_Use_Dependencies.pl ./work/build/Makefile_Include_Deps source/*.[fF]90 source/*.h source/*.c source/*.cpp $(wildcard source/*.Inc)
	@mkdir -p work/build
	./scripts/build/Find_Use_Dependencies.pl `pwd` $(MAKE)

./work/build/Makefile_Directives: ./scripts/build/Code_Directive_Parser.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	@mkdir -p work/build
	./scripts/build/Code_Directive_Parser.pl `pwd`

./work/build/Makefile_Include_Deps: ./scripts/build/Find_Include_Dependencies.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	@mkdir -p work/build
	./scripts/build/Find_Include_Dependencies.pl `pwd`

./work/build/Makefile_All_Execs: ./scripts/build/Find_Programs.pl source/*.[fF]90 source/*.h source/*.c source/*.cpp
	@mkdir -p work/build
	./scripts/build/Find_Programs.pl `pwd`

deps: $(MAKE_DEPS) ./work/build/Makefile_All_Execs
