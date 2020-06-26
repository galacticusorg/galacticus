# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Build option.
GALACTICUS_BUILD_OPTION ?= default
ifeq '$(GALACTICUS_BUILD_OPTION)' 'default'
export BUILDPATH = ./work/build
export SUFFIX =
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
export BUILDPATH = ./work/buildMPI
export SUFFIX =
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'gprof'
export BUILDPATH = ./work/buildGProf
export SUFFIX = _gprof
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'odeprof'
export BUILDPATH = ./work/buildODEProf
export SUFFIX = _odeProf
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
export BUILDPATH = ./work/buildCompileProf
export SUFFIX =
endif

# Preprocessor:
PREPROCESSOR ?= cpp

# Fortran compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
FCCOMPILER ?= mpif90
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
FCCOMPILER ?= /usr/bin/time -f 'Compile statistics (file time memory): $* %e %M' gfortran
else
FCCOMPILER ?= gfortran
endif

# C compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
CCOMPILER ?= mpicc
else
CCOMPILER ?= gcc
endif
export CCOMPILER

# C++ compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
CPPCOMPILER ?= mpic++
else
CPPCOMPILER ?= g++
endif

# Linker for Condor standard universe executables. Uncomment the second line to link for submission to the Condor standard universe.
CONDORLINKER = 
#CONDORLINKER = condor_compile

# Fortran compiler flags:
FCFLAGS += -ffree-line-length-none -frecursive -DBUILDPATH=\'$(BUILDPATH)\' -J$(BUILDPATH)/moduleBuild/ -I$(BUILDPATH)/ ${GALACTICUS_FCFLAGS} -fintrinsic-modules-path /usr/local/finclude -fintrinsic-modules-path /usr/local/include/gfortran -fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/lib/gfortran/modules -fintrinsic-modules-path /usr/include/gfortran -fintrinsic-modules-path /usr/include -fintrinsic-modules-path /usr/finclude -fintrinsic-modules-path /usr/lib64/gfortran/modules -fintrinsic-modules-path /usr/lib64/openmpi/lib -pthread
# Fortran77 compiler flags:
F77FLAGS = ${GALACTICUS_F77FLAGS} -DBUILDPATH=\'$(BUILDPATH)\'
# Error checking flags
FCFLAGS += -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core
# Add bounds checking.
#FCFLAGS += -fbounds-check
# A copy of the flags prior to any optimizations.
FCFLAGS_NOOPT := $(FCFLAGS)
# Optimization flags.
FCFLAGS += -O3 -ffinite-math-only -fno-math-errno
# For OpenMP compilation.
FCFLAGS += -fopenmp

# C compiler flags:
CFLAGS = -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ -fopenmp ${GALACTICUS_CFLAGS}
export CFLAGS

# C++ compiler flags:
CPPFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ -fopenmp ${GALACTICUS_CPPFLAGS}

# Detect GProf compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'gprof'
FCFLAGS       += -pg
FCFLAGS_NOOPT += -pg
CFLAGS        += -pg
CPPFLAGS      += -pg
else
FCFLAGS       += -g
FCFLAGS_NOOPT += -g
CFLAGS        += -g
CPPFLAGS      += -g
endif

# Detect ODE profiling compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'odeprof'
FCFLAGS       += -DPROFILE
FCFLAGS_NOOPT += -DPROFILE
CFLAGS        += -DPROFILE
CPPFLAGS      += -DPROFILE
endif

# Detect MPI compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
FCFLAGS       += -DUSEMPI
FCFLAGS_NOOPT += -DUSEMPI
CFLAGS        += -DUSEMPI
CPPFLAGS      += -DUSEMPI 
endif

# Detect YEPPP libraries.
# ifdef YEPROOT
# FCFLAGS += -I$(YEPROOT)/bindings/fortran/modules/$(YEPPLATFORM)-gfortran/ -L$(YEPBINARIES) -DYEPPP
# endif

# OFD locks option.
OFDLOCKS ?= enabled
ifeq '$(OFDLOCKS)' 'enabled'
CFLAGS   += -DOFDLOCKS
CPPFLAGS += -DOFDLOCKS 
else
CFLAGS   += -DNOOFDLOCKS
CPPFLAGS += -DNOOFDLOCKS 
endif

# Object debugging.
GALACTICUS_OBJECTS_DEBUG ?= no
ifeq '$(GALACTICUS_OBJECTS_DEBUG)' 'yes'
FCFLAGS += -DOBJECTDEBUG
endif

# List of additional Makefiles which contain dependency information
MAKE_DEPS = $(BUILDPATH)/Makefile_Module_Dependencies $(BUILDPATH)/Makefile_Use_Dependencies $(BUILDPATH)/Makefile_Include_Dependencies

# List of files which must always be checked for update.
UPDATE_DEPS = $(BUILDPATH)/allocatableArrays.xml.up

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
$(BUILDPATH)/%.p.F90.up : source/%.F90 $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml
	./scripts/build/preprocess.pl source/$*.F90 $(BUILDPATH)/$*.p.F90
$(BUILDPATH)/%.p.F90 : | $(BUILDPATH)/%.p.F90.up
	@true
$(BUILDPATH)/%.o : $(BUILDPATH)/%.p.F90 $(BUILDPATH)/%.m $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $(BUILDPATH)/$*.p.F90 -o $(BUILDPATH)/$*.o $(FCFLAGS) 2>&1 | ./scripts/build/postprocess.pl $(BUILDPATH)/$*.p.F90
	@mlist=`cat $(BUILDPATH)/$*.m` ; \
	for mod in $$mlist ; \
	do \
         if [ -f $$mod ] ; then \
	  if cmp -s $$mod $(BUILDPATH)/moduleBuild/`basename $$mod`; then \
	   rm $(BUILDPATH)/moduleBuild/`basename $$mod`; \
	  else \
	   mv $(BUILDPATH)/moduleBuild/`basename $$mod` $(BUILDPATH); \
	  fi \
         else \
	  mv $(BUILDPATH)/moduleBuild/`basename $$mod` $(BUILDPATH); \
         fi \
	done

# Rules for building HDF5 C interoperability types data file.
$(BUILDPATH)/hdf5FCInterop.dat  : $(BUILDPATH)/hdf5FCInterop.exe $(BUILDPATH)/hdf5FCInteropC.exe
	$(BUILDPATH)/hdf5FCInterop.exe  >  $(BUILDPATH)/hdf5FCInterop.dat
	$(BUILDPATH)/hdf5FCInteropC.exe >> $(BUILDPATH)/hdf5FCInterop.dat
$(BUILDPATH)/hdf5FCInterop.exe  : source/hdf5FCInterop.F90
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) source/hdf5FCInterop.F90 -o $(BUILDPATH)/hdf5FCInterop.exe $(FCFLAGS)
$(BUILDPATH)/hdf5FCInteropC.exe : source/hdf5FCInteropC.c
	$(CCOMPILER) source/hdf5FCInteropC.c -o $(BUILDPATH)/hdf5FCInteropC.exe $(CFLAGS)

# Configuration of file locking implementation.
-include $(BUILDPATH)/Makefile_Config_OFD
$(BUILDPATH)/Makefile_Config_OFD: source/flock_config.c
	@mkdir -p $(BUILDPATH)
	$(CCOMPILER) -c source/flock_config.c -o $(BUILDPATH)/flock_config.o $(CFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DOFDAVAIL"   >  $(BUILDPATH)/Makefile_Config_OFD ; \
	 echo "CFLAGS   += -DOFDAVAIL"   >> $(BUILDPATH)/Makefile_Config_OFD ; \
	 echo "CPPFLAGS += -DOFDAVAIL"   >> $(BUILDPATH)/Makefile_Config_OFD ; \
	else \
	 echo "FCFLAGS  += -DOFDUNAVAIL" >  $(BUILDPATH)/Makefile_Config_OFD ; \
	 echo "CFLAGS   += -DOFDUNAVAIL" >> $(BUILDPATH)/Makefile_Config_OFD ; \
	 echo "CPPFLAGS += -DOFDUNAVAIL" >> $(BUILDPATH)/Makefile_Config_OFD ; \
	fi

# Configuration for availability of FFTW3.
-include $(BUILDPATH)/Makefile_Config_FFTW3
$(BUILDPATH)/Makefile_Config_FFTW3: source/fftw3_config.F90
	@mkdir -p $(BUILDPATH)
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c source/fftw3_config.F90 -o $(BUILDPATH)/fftw3_config.o $(FCFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS += -DFFTW3AVAIL"   > $(BUILDPATH)/Makefile_Config_FFTW3 ; \
	else \
	 echo "FCFLAGS += -DFFTW3UNAVAIL" > $(BUILDPATH)/Makefile_Config_FFTW3 ; \
	fi

# Configuration for availability of ANN.
-include $(BUILDPATH)/Makefile_Config_ANN
$(BUILDPATH)/Makefile_Config_ANN: source/ann_config.cpp
	@mkdir -p $(BUILDPATH)
	$(CPPCOMPILER) -c source/ann_config.cpp -o $(BUILDPATH)/ann_config.o $(CPPFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DANNAVAIL"   >  $(BUILDPATH)/Makefile_Config_ANN ; \
	 echo "CPPFLAGS += -DANNAVAIL"   >> $(BUILDPATH)/Makefile_Config_ANN ; \
	else \
	 echo "FCFLAGS  += -DANNUNAVAIL" >  $(BUILDPATH)/Makefile_Config_ANN ; \
	 echo "CPPFLAGS += -DANNUNAVAIL" >> $(BUILDPATH)/Makefile_Config_ANN ; \
	fi

# Configuration for availability of libmatheval.
-include $(BUILDPATH)/Makefile_Config_MathEval
$(BUILDPATH)/Makefile_Config_MathEval: source/libmatheval_config.cpp
	@mkdir -p $(BUILDPATH)
	$(CPPCOMPILER) -c source/libmatheval_config.cpp -o $(BUILDPATH)/libmatheval_config.o $(CPPFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DMATHEVALAVAIL"   >  $(BUILDPATH)/Makefile_Config_MathEval ; \
	 echo "CPPFLAGS += -DMATHEVALAVAIL"   >> $(BUILDPATH)/Makefile_Config_MathEval ; \
	else \
	 echo "FCFLAGS  += -DMATHEVALUNAVAIL" >  $(BUILDPATH)/Makefile_Config_MathEval ; \
	 echo "CPPFLAGS += -DMATHEVALUNAVAIL" >> $(BUILDPATH)/Makefile_Config_MathEval ; \
	fi

# Object (*.o) files are built by compiling C (*.c) source files.
vpath %.c source
$(BUILDPATH)/%.o : %.c $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CFLAGS)

# Object (*.o) can also be built from C++ source files.
vpath %.cpp source
$(BUILDPATH)/%.o : %.cpp $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CPPCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CPPFLAGS)

# Rules for FFTLog library.
source/FFTlog/cdgamma.f source/FFTlog/drfftb.f source/FFTlog/drffti.f source/FFTlog/drfftf.f: source/FFTlog/fftlog.f
source/FFTlog/fftlog.f:
	mkdir -p source/FFTlog
	mkdir -p $(BUILDPATH)/FFTlog
	wget http://jila.colorado.edu/~ajsh/FFTLog/fftlog.tgz -O - | tar xvz -C source/FFTlog -f -
	if [ ! -e source/FFTlog/fftlog.f ]; then \
	 echo "      subroutine fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)" >  source/FFTlog/fftlog.f; \
	 echo "      stop 'FFTlog was not downloaded - to try again" >> source/FFTlog/fftlog.f; \
	 echo "     & remove the source/FFTlog directory'"           >> source/FFTlog/fftlog.f; \
	 echo "      end subroutine fhti"                            >> source/FFTlog/fftlog.f; \
	 touch source/FFTlog/cdgamma.f; \
	 touch source/FFTlog/drfftb.f; \
	 touch source/FFTlog/drfftf.f; \
	 touch source/FFTlog/drffti.f; \
	else \
	 cd source/FFTlog; \
	 patch < ../drfftb.f.patch; \
	 patch < ../drfftf.f.patch; \
	 patch < ../drffti.f.patch; \
	 cd -; \
	 ./scripts/build/useDependencies.pl `pwd`; \
	fi
	echo $(BUILDPATH)/FFTlog/cdgamma.o > $(BUILDPATH)/FFTlog/cdgamma.d
	echo $(BUILDPATH)/FFTlog/drfftb.o  > $(BUILDPATH)/FFTlog/drfftb.d
	echo $(BUILDPATH)/FFTlog/drfftf.o  > $(BUILDPATH)/FFTlog/drfftf.d
	echo $(BUILDPATH)/FFTlog/drffti.o  > $(BUILDPATH)/FFTlog/drffti.d
	echo $(BUILDPATH)/FFTlog/fftlog.o  > $(BUILDPATH)/FFTlog/fftlog.d

$(BUILDPATH)/FFTlog/%.o: ./source/FFTlog/%.f Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/FFTlog/$*.o $(F77FLAGS) -Wno-argument-mismatch -std=legacy

# Object (*.o) files are built by compiling Fortran (*.f) source files.
vpath %.f source
$(BUILDPATH)/%.o : %.f $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(F77FLAGS)

# Special rules required for building some sources (unfortunate, but necessary....)
# bivar.F90 doesn't like to be compiled with any optimization:
$(BUILDPATH)/Bivar/bivar.o : ./source/Bivar/bivar.F90 Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/Bivar/bivar.o $(FCFLAGS_NOOPT)
	@if [ -f $(BUILDPATH)/bivar.mod ] ; then \
	 if cmp -s $(BUILDPATH)/bivar.mod $(BUILDPATH)/moduleBuild/bivar.mod; then \
	  rm $(BUILDPATH)/moduleBuild/bivar.mod; \
	 else \
	  mv $(BUILDPATH)/moduleBuild/bivar.mod $(BUILDPATH); \
	 fi \
        else \
	 mv $(BUILDPATH)/moduleBuild/bivar.mod $(BUILDPATH); \
        fi \
# pfq.new.f
$(BUILDPATH)/pFq/pfq.new.o : ./source/pFq/pfq.new.f Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/pFq/pfq.new.o $(FCFLAGS)

# Rule for running *.Inc files through the preprocessor. We strip out single quote characters in comment lines to avoid spurious
# complaints from the preprocessor.
$(BUILDPATH)/%.Inc.up : ./source/%.Inc $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml
	./scripts/build/preprocess.pl ./source/$*.Inc $(BUILDPATH)/$*.Inc
$(BUILDPATH)/%.Inc : | $(BUILDPATH)/%.Inc.up
	@true
$(BUILDPATH)/%.inc : $(BUILDPATH)/%.Inc Makefile
	perl -MRegexp::Common -ne '$$l=$$_;$$l =~ s/($$RE{comment}{Fortran}{-keep})/\/\*$$4\*\/$$5/; print $$l' $< | cpp -nostdinc -C | perl -MRegexp::Common -ne '$$l=$$_;$$l =~ s/($$RE{comment}{C}{-keep})/!$$4/; print $$l' > $(BUILDPATH)/$*.tmp
	mv -f $(BUILDPATH)/$*.tmp $(BUILDPATH)/$*.inc

# Dependency files (*.d) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Dependencies Makefile_Module_Dependencies files, but this acts as a fallback rule.
$(BUILDPATH)/%.d : ./source/%.F90
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d~
	@if cmp -s $(BUILDPATH)/$*.d $(BUILDPATH)/$*.d~ ; then \
	 rm $(BUILDPATH)/$*.d~ ; \
	else \
	 mv $(BUILDPATH)/$*.d~ $(BUILDPATH)/$*.d ; \
	fi
$(BUILDPATH)/%.d : ./source/%.f
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d~
	@if cmp -s $(BUILDPATH)/$*.d $(BUILDPATH)/$*.d~ ; then \
	 rm $(BUILDPATH)/$*.d~ ; \
	else \
	 mv $(BUILDPATH)/$*.d~ $(BUILDPATH)/$*.d ; \
	fi
$(BUILDPATH)/%.d : ./source/%.c
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d~
	@if cmp -s $(BUILDPATH)/$*.d $(BUILDPATH)/$*.d~ ; then \
	 rm $(BUILDPATH)/$*.d~ ; \
	else \
	 mv $(BUILDPATH)/$*.d~ $(BUILDPATH)/$*.d ; \
	fi
$(BUILDPATH)/%.d : ./source/%.cpp
	@echo $(BUILDPATH)/$*.o > $(BUILDPATH)/$*.d~
	@if cmp -s $(BUILDPATH)/$*.d $(BUILDPATH)/$*.d~ ; then \
	 rm $(BUILDPATH)/$*.d~ ; \
	else \
	 mv $(BUILDPATH)/$*.d~ $(BUILDPATH)/$*.d ; \
	fi
%.d : %.f
	@echo $*.o > $*.d~
	@if cmp -s $*.d $*.d~ ; then \
	 rm $*.d~ ; \
	else \
	 mv $*.d~ $*.d ; \
	fi
%.d :
	@if [ ! -f $*.d ]; then \
	 mkdir -p `dirname $*.d`; \
	 touch $*.d; \
	fi

# In some instances a C header file may not be required due to being preprocessed out, but a dependency for it still appears in
# the Makefiles. In these instances we just create an empty file, since it won't be used anyway.
%.h :
	@if [ ! -f $*.h ]; then \
	 mkdir -p `dirname $*.h`; \
	 touch $*.h; \
	fi

# Rules for making update (".up") files if no explicit rule is given.
%.up :
	touch $*.up

# Library files (*.fl) are created as empty files by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Dependencies file, but this acts as a fallback rule.
$(BUILDPATH)/%.fl : ./source/%.F90
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.f
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.c
	@touch $(BUILDPATH)/$*.fl
$(BUILDPATH)/%.fl : ./source/%.cpp
	@touch $(BUILDPATH)/$*.fl

# GraphViz files (*.gv) are created with just a node entry by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Use_Dependencies Makefile_Module_Dependencies files, but this acts as a fallback rule.
$(BUILDPATH)/%.F90.gv : ./source/%.F90
	@echo \"$*.F90\" > $(BUILDPATH)/$*.F90.gv
$(BUILDPATH)/%.c.gv : ./source/%.c
	@echo \"$*.c\" > $(BUILDPATH)/$*.c.gv
$(BUILDPATH)/%.cpp.gv : ./source/%.cpp
	@echo \"$*.cpp\" > $(BUILDPATH)/$*.cpp.gv

# Create a PostScript files showing tree diagrams of source file dependencies.
%.tpdf : $(BUILDPATH)/%.gv
	@echo digraph Tree \{ > $*.dot
	@cat $(BUILDPATH)/$*.gv >> $*.dot
	@echo \} >> $*.dot
	dot -Tpdf $*.dot -o $*.pdf
	@rm $*.dot

# Module list files are created empty by default. Normally this rule is overruled by a specific set of rules in the
# Makefile_Module_Dependencies files, but this acts as a fallback rule.
$(BUILDPATH)/%.m : ./source/%.F90
	@touch $(BUILDPATH)/$*.m

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the associated dependency (*.d)
# file.
%.exe: $(BUILDPATH)/%.o $(BUILDPATH)/%.d `cat $(BUILDPATH)/$*.d` $(MAKE_DEPS) $(UPDATE_DEPS)
	./scripts/build/parameterDependencies.pl `pwd` $*.exe
	$(FCCOMPILER) -c $(BUILDPATH)/$*.parameters.F90 -o $(BUILDPATH)/$*.parameters.o $(FCFLAGS)
	./scripts/build/sourceDigests.pl `pwd` $*.exe
	$(CCOMPILER) -c $(BUILDPATH)/$*.md5s.c -o $(BUILDPATH)/$*.md5s.o $(CFLAGS)
	$(CONDORLINKER) $(FCCOMPILER) `cat $*.d` $(BUILDPATH)/$*.parameters.o $(BUILDPATH)/$*.md5s.o -o $*.exe$(SUFFIX) $(FCFLAGS) `scripts/build/libraryDependencies.pl $*.exe $(FCFLAGS)`

# Ensure that we don't delete object files which make considers to be intermediate
.PRECIOUS: $(BUILDPATH)/%.p.F90 $(BUILDPATH)/%.p.F90.up $(BUILDPATH)/%.Inc $(BUILDPATH)/%.Inc.up $(BUILDPATH)/%.d

# Cancel all builtin rules.
.SUFFIXES:

# Include depenencies on "include" files.
-include $(BUILDPATH)/Makefile_Include_Dependencies 

# Include module dependencies.
-include $(BUILDPATH)/Makefile_Module_Dependencies

# Include rules to build include files generated from directives.
-include $(BUILDPATH)/Makefile_Directives

# Include module use dependencies. Include this after Makefile_Directives, as Makefile_Directives will
# specify dependencies for Makefile_Use_Dependencies
-include $(BUILDPATH)/Makefile_Use_Dependencies

# Rules for memory management routines.
$(BUILDPATH)/allocatableArrays.xml.up: ./scripts/build/allocatableArrays.pl source/*.[fF]90 $(wildcard source/*.Inc)
	./scripts/build/allocatableArrays.pl `pwd`
$(BUILDPATH)/allocatableArrays.xml: | $(BUILDPATH)/allocatableArrays.xml.up

$(BUILDPATH)/utility.memory_management.preContain.inc: ./scripts/build/memoryManagementFunctions.pl $(BUILDPATH)/allocatableArrays.xml
	./scripts/build/memoryManagementFunctions.pl

$(BUILDPATH)/utility.memory_management.postContain.inc:
	@touch $(BUILDPATH)/utility.memory_management.postContain.inc

$(BUILDPATH)/openMPCriticalSections.count.inc $(BUILDPATH)/openMPCriticalSections.enumerate.inc: $(BUILDPATH)/openMPCriticalSections.xml
	@touch $(BUILDPATH)/openMPCriticalSections.count.inc $(BUILDPATH)/openMPCriticalSections.enumerate.inc
$(BUILDPATH)/openMPCriticalSections.xml: ./scripts/build/enumerateOpenMPCriticalSections.pl
	./scripts/build/enumerateOpenMPCriticalSections.pl `pwd`

# Rules for version routines.
$(BUILDPATH)/galacticus.output.version.revision.inc: $(wildcard .git/refs/heads/master)
	@if [ -f .git/refs/heads/master ] ; then git rev-parse HEAD | awk '{print "character(len=42), parameter :: gitHash=\""$$1"\""}' > $(BUILDPATH)/galacticus.output.version.revision.inc; else printf 'character(len=42), parameter :: gitHash="(unknown)"\n' > $(BUILDPATH)/galacticus.output.version.revision.inc; fi
	@if [ -f .git/refs/heads/master ] ; then git branch | awk '{if ($$1 == "*") print "character(len=128), parameter :: gitBranch=\""$$2"\""}' >> $(BUILDPATH)/galacticus.output.version.revision.inc; else printf 'character(len=128), parameter :: gitBranch="(unknown)"\n' >> $(BUILDPATH)/galacticus.output.version.revision.inc; fi
	@date -u '+%a %b %d %k:%M:%S UTC %Y' | awk '{print "character(len=32), parameter :: buildTime=\""$$0"\""}' >> $(BUILDPATH)/galacticus.output.version.revision.inc

# Rules for build information routines.
$(BUILDPATH)/galacticus.output.build.environment.inc:
	@echo PREPROCESSOR=\"$(PREPROCESSOR)\" > $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCCOMPILER=\"$(FCCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CCOMPILER=\"$(CCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPCOMPILER=\"$(CPPCOMPILER)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCFLAGS=\"$(FCFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCFLAGS_NOOPT=\"$(FCFLAGS_NOOPT)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CFLAGS=\"$(CFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPFLAGS=\"$(CPPFLAGS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo LIBS=\"$(LIBS)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo FCCOMPILER_VERSION=\"$(FCCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CCOMPILER_VERSION=\"$(CCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc
	@echo CPPCOMPILER_VERSION=\"$(CPPCOMPILER_VERSION)\" >> $(BUILDPATH)/galacticus.output.build.environment.inc

# Rules for changeset creation.
Galacticus.exe: $(BUILDPATH)/galacticus.git.patch $(BUILDPATH)/galacticus.git.bundle
$(BUILDPATH)/galacticus.git.patch:
	git diff > $(BUILDPATH)/galacticus.git.patch || echo unknown > $(BUILDPATH)/galacticus.git.patch
$(BUILDPATH)/galacticus.git.bundle:
	git bundle create $(BUILDPATH)/galacticus.git.bundle HEAD ^origin > /dev/null 2>&1 || echo unknown > $(BUILDPATH)/galacticus.git.bundle

# Rules for cleaning up.
clean: tidy
	rm -f *.exe

tidy:
	rm -rf $(BUILDPATH)/*

# Rule for making all executables.
-include $(BUILDPATH)/Makefile_All_Execs
all: deps $(all_exes)

# Rules for building dependency Makefiles.
$(BUILDPATH)/Makefile_Module_Dependencies: ./scripts/build/moduleDependencies.pl $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies source/*.[fF]90 $(wildcard source/*.h) source/*.c $(wildcard source/*.cpp) $(wildcard source/*.Inc)
	@mkdir -p $(BUILDPATH)
	./scripts/build/moduleDependencies.pl `pwd`

$(BUILDPATH)/Makefile_Use_Dependencies: ./scripts/build/useDependencies.pl $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies source/*.[fF]90 $(wildcard source/*.h) source/*.c $(wildcard source/*.cpp) $(wildcard source/*.Inc)
	@mkdir -p $(BUILDPATH)
	./scripts/build/useDependencies.pl `pwd`

$(BUILDPATH)/Makefile_Directives: ./scripts/build/codeDirectivesParse.pl source/*.[fF]90 $(wildcard source/*.h) source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/codeDirectivesParse.pl `pwd`
	./scripts/build/stateStorables.pl `pwd`
	./scripts/build/deepCopyActions.pl `pwd`

$(BUILDPATH)/Makefile_Include_Dependencies: ./scripts/build/includeDependencies.pl source/*.[fF]90 $(wildcard source/*.h) source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/includeDependencies.pl `pwd`

$(BUILDPATH)/Makefile_All_Execs: ./scripts/build/findExecutables.pl source/*.[fF]90 $(wildcard source/*.h) source/*.c $(wildcard source/*.cpp)
	@mkdir -p $(BUILDPATH)
	./scripts/build/findExecutables.pl `pwd`

deps: $(MAKE_DEPS) $(BUILDPATH)/Makefile_All_Execs

# Rules for XSpec code.
aux/XSpec/%.o: ./aux/XSpec/%.f Makefile
	$(FCCOMPILER) -c $< -o aux/XSpec/$*.o $(FCFLAGS)
