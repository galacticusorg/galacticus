# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Build option.
GALACTICUS_BUILD_OPTION ?= default
ifeq '$(GALACTICUS_BUILD_OPTION)' 'default'
export BUILDPATH = ./work/build
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
export BUILDPATH = ./work/buildMPI
endif

# Preprocessor:
PREPROCESSOR ?= cpp

# Fortran compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
FCCOMPILER ?= mpif90
else
FCCOMPILER ?= gfortran
endif

# C compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
CCOMPILER ?= mpicc
else
CCOMPILER ?= gcc
endif

# C++ compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
CPPCOMPILER ?= mpic++
else
CPPCOMPILER ?= g++
endif

# Linker for Condor standard universe executables. Uncomment the second line to link for submission to the Condor standard universe.
CONDORLINKER = 
#CONDORLINKER = condor_compile

# Module type (used for checking if module interfaces have changed):
MODULETYPE ?= GCC-f95-on-LINUX

# Fortran compiler flags:
FCFLAGS += -ffree-line-length-none -frecursive -DBUILDPATH=\'$(BUILDPATH)\' -J$(BUILDPATH)/ -I$(BUILDPATH)/ ${GALACTICUS_FCFLAGS} -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules -fintrinsic-modules-path /usr/lib64/openmpi/lib -pthread
# Fortran77 compiler flags:
F77FLAGS = -g -DBUILDPATH=\'$(BUILDPATH)\'
# Error checking flags
FCFLAGS += -Wall -g -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core
# Add bounds checking.
#FCFLAGS += -fbounds-check
# Add profiling.
FCFLAGS += -g
# A copy of the flags prior to any optimizations.
FCFLAGS_NOOPT := $(FCFLAGS)
# Optimization flags.
FCFLAGS += -O3 -ffinite-math-only -fno-math-errno
# For OpenMP compilation.
FCFLAGS += -fopenmp

# C compiler flags:
CFLAGS = -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ -I/opt/gsl-trunk/include ${GALACTICUS_CFLAGS}
CFLAGS += -g

# C++ compiler flags:
CPPFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ ${GALACTICUS_CPPFLAGS}
CPPFLAGS += -g

# Detect MPI compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
FCFLAGS += -DUSEMPI
CFLAGS += -DUSEMPI
CPPFLAGS += -DUSEMPI 
endif

# Detect YEPPP libraries.
# ifdef YEPROOT
# FCFLAGS += -I$(YEPROOT)/bindings/fortran/modules/$(YEPPLATFORM)-gfortran/ -L$(YEPBINARIES) -DYEPPP
# endif

# List of additional Makefiles which contain dependency information
MAKE_DEPS = $(BUILDPATH)/Makefile_Module_Deps $(BUILDPATH)/Makefile_Use_Deps $(BUILDPATH)/Makefile_Include_Deps

# Get versions of build tools.
FCCOMPILER_VERSION = `$(FCCOMPILER) -v 2>&1`
CCOMPILER_VERSION = `$(CCOMPILER) -v 2>&1`
CPPCOMPILER_VERSION = `$(CPPCOMPILER) -v 2>&1`

# General suffix rules: i.e. rules for making a file of one suffix from files of another suffix.

# Object (*.o) files are built by preprocessing and then compiling Fortran 90 (*.F90) source
# files. Note that .F90 source files should not have names which coincide with the name of a
# module - this will lead to circular dependency problems as Make becomes confused about how to
# build the module file.
vpath %.F90 source
$(BUILDPATH)/%.p.F90 : source/%.F90
	scripts/build/preprocess.pl source/$*.F90 $(BUILDPATH)/$*.p.F90
$(BUILDPATH)/%.o : $(BUILDPATH)/%.p.F90 $(BUILDPATH)/%.m $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mlist=`cat $(BUILDPATH)/$*.m` ; \
	for mod in $$mlist ; \
	do \
	 if [ -f $$mod ] ; then mv $$mod $$mod~; fi \
	done
	$(FCCOMPILER) -c $(BUILDPATH)/$*.p.F90 -o $(BUILDPATH)/$*.o $(FCFLAGS) 2>&1 | ./scripts/build/postprocess.pl $(BUILDPATH)/$*.p.F90
	@mlist=`cat $(BUILDPATH)/$*.m` ; \
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
$(BUILDPATH)/%.o : %.c $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CFLAGS)

# Object (*.o) can also be built from C++ source files.
vpath %.cpp source
$(BUILDPATH)/%.o : %.cpp $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CPPCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CPPFLAGS)

# Object (*.o) files are built by compiling Fortran (*.f) source files.
vpath %.f source
$(BUILDPATH)/%.o : %.f $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(F77FLAGS)

# Special rules required for building some sources (unfortunate, but necessary....)
# bivar.F90 doesn't like to be compiled with any optimization:
$(BUILDPATH)/Bivar/bivar.o : ./source/Bivar/bivar.F90 Makefile
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/Bivar/bivar.o $(FCFLAGS_NOOPT)
# pfq.new.f
$(BUILDPATH)/pFq/pfq.new.o : ./source/pFq/pfq.new.f Makefile
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/pFq/pfq.new.o $(FCFLAGS)

# Rule for running *.Inc files through the preprocessor. We strip out single quote characters in comment lines to avoid spurious
# complaints from the preprocessor.
$(BUILDPATH)/%.Inc : ./source/%.Inc
	scripts/build/preprocess.pl ./source/$*.Inc $(BUILDPATH)/$*.Inc
$(BUILDPATH)/%.inc : $(BUILDPATH)/%.Inc Makefile
	perl -MRegexp::Common -ne '$$l=$$_;if ( $$l =~ m/($$RE{comment}{Fortran})/ ) {($$m = $$1) =~ s/(?<!\\)'\''//g; $$l =~ s/$$RE{comment}{Fortran}/$$m/}; print $$l' $< | $(PREPROCESSOR) -nostdinc -C -o $(BUILDPATH)/$*.tmp
	mv -f $(BUILDPATH)/$*.tmp $(BUILDPATH)/$*.inc

# Dependency files (*.d) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
$(BUILDPATH)/%.d : ./source/%.F90
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d
$(BUILDPATH)/%.d : ./source/%.f
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d
$(BUILDPATH)/%.d : ./source/%.c
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d
$(BUILDPATH)/%.d : ./source/%.cpp
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d
%.d : %.f
	@echo $*.o > $*.d
%.d :
	@mkdir -p `dirname $*.d`
	@touch $*.d


# Library files (*.fl) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps file, but this acts as a fallback rule.
$(BUILDPATH)/%.fl : ./source/%.F90
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.f
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.c
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.cpp
	@touch $(BUILDPATH)/$*.fl

# GraphViz files (*.gv) are created with just a node entry by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Deps Makefile_Module_Deps files, but this acts as a fallback rule.
$(BUILDPATH)/%.F90.gv : ./source/%.F90
	@echo \"$*.F90\" > $(BUILDPATH)/$*.F90.gv
$(BUILDPATH)/%.c.gv : ./source/%.c
	@echo \"$*.c\" > $(BUILDPATH)/$*.c.gv
$(BUILDPATH)/%.cpp.gv : ./source/%.cpp
	@echo \"$*.cpp\" > $(BUILDPATH)/$*.cpp.gv

# Create a PostScript files showing tree diagrams of source file dependencies.
%.tps : $(BUILDPATH)/%.gv
	@echo digraph Tree \{ > $*.dot
	@cat $(BUILDPATH)/$*.gv >> $*.dot
	@echo \} >> $*.dot
	dot -Tps $*.dot -o $*.tps
	@rm $*.dot

# Module list files are created empty by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Module_Deps files, but this acts as a fallback rule.
$(BUILDPATH)/%.m : ./source/%.F90
	@touch $(BUILDPATH)/$*.m

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the associated dependency (*.d)
# file.
%.exe : $(BUILDPATH)/%.o $(BUILDPATH)/%.d `cat $(BUILDPATH)/$*.d` $(MAKE_DEPS)
	 $(CONDORLINKER) $(FCCOMPILER) `cat $*.d` -o $*.exe $(FCFLAGS) `scripts/build/Library_Dependencies.pl $*.exe $(FCFLAGS)`
	 ./scripts/build/Find_Executable_Size.pl $*.exe $*.size
	 ./scripts/build/Find_Parameter_Dependencies.pl `pwd` $*.exe

# Ensure that we don't delete object files which make considers to be intermediate
.PRECIOUS: %.o %.d %.dd %.m %.make %.Inc $(BUILDPATH)/%.p.F90

# Cancel all builtin rules.
.SUFFIXES:

# Include depenencies on "include" files.
-include $(BUILDPATH)/Makefile_Include_Deps 

# Include module dependencies.
-include $(BUILDPATH)/Makefile_Module_Deps

# Include rules to build include files generated from directives.
-include $(BUILDPATH)/Makefile_Directives

# Include module use dependencies. Include this after Makefile_Directives, as Makefile_Directives will
# specify dependencies for Makefile_Use_Deps
-include $(BUILDPATH)/Makefile_Use_Deps

# Rules for memory management routines.
$(BUILDPATH)/Allocatable_Arrays.xml: ./scripts/build/Find_Allocatable_Arrays.pl source/*.[fF]90 $(wildcard source/*.Inc)
	./scripts/build/Find_Allocatable_Arrays.pl `pwd`

$(BUILDPATH)/utility.memory_management.precontain.inc: ./scripts/build/Make_Memory_Usage_Routines.pl $(BUILDPATH)/Allocatable_Arrays.xml
	./scripts/build/Make_Memory_Usage_Routines.pl

$(BUILDPATH)/utility.memory_management.postcontain.inc:
	@touch $(BUILDPATH)/utility.memory_management.postcontain.inc

# Rules for version routines.
$(BUILDPATH)/galacticus.output.version.revision.inc: $(wildcard .hg/branch)
	@if [ -f .hg/branch ] ; then hg tip | awk 'BEGIN {FS=":";r=-1;h=""} {if ((NR == 1 && NF == 3 ) || $$1 == "parent") {r=$$2;h=$$3}} END {print "integer, parameter :: hgRevision="r"\ncharacter(len=12), parameter :: hgHash=\""h"\""}' > $(BUILDPATH)/galacticus.output.version.revision.inc; else printf 'integer, parameter :: hgRevision=-1\ncharacter(len=12), parameter :: hgHash=""\n' > $(BUILDPATH)/galacticus.output.version.revision.inc; fi

# Rules for build information routines.
$(BUILDPATH)/galacticus.output.build.environment.inc:
	@echo PREPROCESSOR=\"$(PREPROCESSOR)\" > $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCCOMPILER=\"$(FCCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CCOMPILER=\"$(CCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPCOMPILER=\"$(CPPCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo MODULETYPE=\"$(MODULETYPE)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCFLAGS=\"$(FCFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCFLAGS_NOOPT=\"$(FCFLAGS_NOOPT)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CFLAGS=\"$(CFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPFLAGS=\"$(CPPFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo LIBS=\"$(LIBS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCCOMPILER_VERSION=\"$(FCCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CCOMPILER_VERSION=\"$(CCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPCOMPILER_VERSION=\"$(CPPCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc

# Rules for unique label function creation.
dfiles := $(patsubst source/%.F90,$(BUILDPATH)/%.d,$(wildcard source/*.F90))
mfiles := $(patsubst source/%.F90,$(BUILDPATH)/%.m,$(wildcard source/*.F90))
$(BUILDPATH)/utility.input_parameters.unique_labels.inc:
	@touch $(BUILDPATH)/utility.input_parameters.unique_labels.inc
$(BUILDPATH)/utility.input_parameters.unique_labels.visibilities.inc: $(dfiles) $(mfiles)
	./scripts/build/Make_Unique_Label_Functions.pl `pwd`

# Rules for changeset creation.
Galacticus.exe: $(BUILDPATH)/galacticus.hg.patch $(BUILDPATH)/galacticus.hg.bundle
$(BUILDPATH)/galacticus.hg.patch:
	hg diff > $(BUILDPATH)/galacticus.hg.patch | true
$(BUILDPATH)/galacticus.hg.bundle:
	hg bundle -t none $(BUILDPATH)/galacticus.hg.bundle https://abensonca@bitbucket.org/abensonca/galacticus | true

# Rules for cleaning up.
clean: tidy
	rm -f *.exe

tidy:
	rm -rf $(BUILDPATH)/*

# Rule for making all executables.
-include $(BUILDPATH)/Makefile_All_Execs
all: deps $(all_exes)

# Rules for building dependency Makefiles.
$(BUILDPATH)/Makefile_Module_Deps: ./scripts/build/Find_Module_Dependencies.pl source/*.[fF]90 source/*.h source/*.c $(wildcard source/*.cpp) $(wildcard source/*.Inc)
	@mkdir -p $(BUILDPATH)
	./scripts/build/Find_Module_Dependencies.pl `pwd`

$(BUILDPATH)/Makefile_Use_Deps: ./scripts/build/Find_Use_Dependencies.pl $(BUILDPATH)/Code_Directive_Locations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Deps source/*.[fF]90 source/*.h source/*.c $(wildcard source/*.cpp) $(wildcard source/*.Inc)
	@mkdir -p $(BUILDPATH)
	./scripts/build/Find_Use_Dependencies.pl `pwd` $(MAKE)

$(BUILDPATH)/Makefile_Directives: ./scripts/build/Code_Directive_Parser.pl source/*.[fF]90 source/*.h source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/Code_Directive_Parser.pl `pwd`

$(BUILDPATH)/Makefile_Include_Deps: ./scripts/build/Find_Include_Dependencies.pl source/*.[fF]90 source/*.h source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/Find_Include_Dependencies.pl `pwd`

$(BUILDPATH)/Makefile_All_Execs: ./scripts/build/Find_Programs.pl source/*.[fF]90 source/*.h source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/Find_Programs.pl `pwd`

deps: $(MAKE_DEPS) $(BUILDPATH)/Makefile_All_Execs

# Rules for FFTLog library.
source/FFTlog/fftlog.f:
	mkdir -p source/FFTlog
	wget http://casa.colorado.edu/~ajsh/FFTLog/fftlog.tar.gz -O - | tar xvz -C source/FFTlog -f -
	if [ ! -e source/FFTlog/fftlog.f ]; then \
	 echo "subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)" > source/FFTlog/fftlog.f; \
	 echo " stop 'FFTlog was not downloaded'" > source/FFTlog/fftlog.f; \
	 echo "end subroutine fhti" >> source/FFTlog/fftlog.f; \
	 touch source/FFTlog/cdgamma.f; \
	 touch source/FFTlog/drfftb.f; \
	 touch source/FFTlog/drfftf.f; \
	 touch source/FFTlog/drffti.f; \
	fi

