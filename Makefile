# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Preprocessor:
PREPROCESSOR = cpp

# Fortran compiler:
F03COMPILER = gfortran

# Module type (used for checking if module interfaces have changed):
MODULETYPE = GCC-f95-on-LINUX

# Fortran compiler flags:
F03FLAGS = -ffree-line-length-none -frecursive -J./work/build/ -I./work/build/ -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -Wall -fbounds-check -g -fbacktrace -ffpe-trap=invalid
## Currently OpenMP seems to cause problems with procedure pointers in gfortran.
##F03FLAGS += -fopenmp # For OpenMP compilation.

# Libraries:
LIBS = -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys -lfgsl_gfortran -lgsl -lgslcblas -lm -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran

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

# Rule for running *.Inc files through the preprocessor.
./work/build/%.inc : ./work/build/%.Inc Makefile
	$(PREPROCESSOR) $< ./work/build/$*.inc

# Dependency files (*.d) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.d : ./source/%.F90
	@echo ./work/build/$*.o > ./work/build/$*.d

# Module list files are created empty by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Module_Deps files, but this acts as a fallback rule.
./work/build/%.m : ./source/%.F90
	@touch ./work/build/$*.m

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the associated dependency (*.d)
# file.
%.exe : ./work/build/%.o ./work/build/%.d `cat ./work/build/$*.d` $(MAKE_DEPS)
	 $(F03COMPILER) `cat $*.d` -o $*.exe $(F03FLAGS) $(LIBS)
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
./work/build/Allocatable_Arrays.xml: ./scripts/build/Find_Allocatable_Arrays.pl source/*.[fF]90 $(wildcard source/*.inc90)
	./scripts/build/Find_Allocatable_Arrays.pl `pwd`

./work/build/Type_Definitions.xml: ./scripts/build/Find_Type_Dependencies.pl source/*.[fF]90 $(wildcard source/*.inc90)
	./scripts/build/Find_Type_Dependencies.pl `pwd`

./work/build/utility.memory_management.use.inc: ./scripts/build/Make_Memory_Usage_Routines.pl ./work/build/Allocatable_Arrays.xml ./work/build/Type_Definitions.xml
	./scripts/build/Make_Memory_Usage_Routines.pl
	./scripts/build/Find_Use_Dependencies.pl .

./work/build/utility.memory_management.precontain.inc:
	@touch ./work/build/utility.memory_usage.precontain.inc

./work/build/utility.memory_management.postcontain.inc:
	@touch ./work/build/utility.memory_usage.postcontain.inc

# Rules for cleaning up.
clean: tidy
	rm -f *.exe

tidy:
	rm -rf work/build/*

# Rule for making all executables.
-include ./work/build/Makefile_All_Execs
all: deps $(all_exes)

# Rules for building dependency Makefiles.
./work/build/Makefile_Module_Deps: ./scripts/build/Find_Module_Dependencies.pl source/*.[fF]90 $(wildcard source/*.inc90)
	./scripts/build/Find_Module_Dependencies.pl `pwd`

./work/build/Makefile_Use_Deps: ./scripts/build/Find_Use_Dependencies.pl source/*.[fF]90 $(wildcard source/*.inc90)
	./scripts/build/Find_Use_Dependencies.pl `pwd` $(MAKE)

./work/build/Makefile_Directives: ./scripts/build/Code_Directive_Parser.pl source/*.[fF]90 $(wildcard source/*.inc90)
	./scripts/build/Code_Directive_Parser.pl `pwd`

./work/build/Makefile_Include_Deps: ./scripts/build/Find_Include_Dependencies.pl source/*.[fF]90
	./scripts/build/Find_Include_Dependencies.pl `pwd`

./work/build/Makefile_All_Execs: ./scripts/build/Find_Programs.pl source/*.[fF]90
	./scripts/build/Find_Programs.pl `pwd`

deps: $(MAKE_DEPS) ./work/build/Makefile_All_Execs
