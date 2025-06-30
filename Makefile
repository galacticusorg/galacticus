# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Build option.
GALACTICUS_BUILD_OPTION ?= default
ifeq '$(GALACTICUS_BUILD_OPTION)' 'default'
export BUILDPATH ?= ./work/build
export SUFFIX ?=
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
export BUILDPATH ?= ./work/buildMPI
export SUFFIX ?=
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'lib'
export BUILDPATH ?= ./work/buildLib
export SUFFIX ?=_lib
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'gprof'
export BUILDPATH ?= ./work/buildGProf
export SUFFIX ?= _gprof
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'odeprof'
export BUILDPATH ?= ./work/buildODEProf
export SUFFIX ?= _odeProf
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
export BUILDPATH ?= ./work/build
export SUFFIX ?=
endif

# Preprocessor:
PREPROCESSOR ?= cpp

# Profiling options.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
SHELL = ./scripts/build/profiler.sh
endif

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
FCFLAGS += -ffree-line-length-none -frecursive -DBUILDPATH=\'$(BUILDPATH)\' -J$(BUILDPATH)/moduleBuild/ -I$(BUILDPATH)/ ${GALACTICUS_FCFLAGS} -pthread
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
FCFLAGS  += -fopenmp
CFLAGS   += -fopenmp
CPPFLAGS += -fopenmp
# Detect static compilation
STATIC=$(findstring -static,${FCFLAGS})
ifeq '${STATIC}' '-static'
FCFLAGS  += -DSTATIC
CFLAGS   += -DSTATIC
CPPFLAGS += -DSTATIC
endif

# C compiler flags:
CFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ ${GALACTICUS_CFLAGS}
export CFLAGS

# C++ compiler flags:
CPPFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' -I./source/ -I$(BUILDPATH)/ ${GALACTICUS_CPPFLAGS}

# Detect library compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'lib'
FCFLAGS       += -fPIC
FCFLAGS_NOOPT += -fPIC
F77FLAGS      += -fPIC
CFLAGS        += -fPIC
CPPFLAGS      += -fPIC
endif

# Detect GProf compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'gprof'
FCFLAGS       += -pg
FCFLAGS_NOOPT += -pg
F77FLAGS      += -pg
CFLAGS        += -pg
CPPFLAGS      += -pg
else
FCFLAGS       += -g
FCFLAGS_NOOPT += -g
F77FLAGS      += -g
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
MAKE_DEPS = $(BUILDPATH)/Makefile_Module_Dependencies $(BUILDPATH)/Makefile_Use_Dependencies $(BUILDPATH)/Makefile_Include_Dependencies $(BUILDPATH)/Makefile_Library_Dependencies

# Get versions of build tools.
FCCOMPILER_VERSION = `$(FCCOMPILER) -v 2>&1`
CCOMPILER_VERSION = `$(CCOMPILER) -v 2>&1`
CPPCOMPILER_VERSION = `$(CPPCOMPILER) -v 2>&1`

# Find all source files.
ALLSOURCES    = $(wildcard source/*.[fF]90              source/*.h source/*.c source/*.cpp)
ALLSOURCESINC = $(wildcard source/*.[fF]90 source/*.Inc source/*.h source/*.c source/*.cpp)

# General suffix rules: i.e. rules for making a file of one suffix from files of another suffix.

# Object (*.o) files are built by preprocessing and then compiling Fortran 90 (*.F90) source
# files. Note that .F90 source files should not have names which coincide with the name of a
# module - this will lead to circular dependency problems as Make becomes confused about how to
# build the module file.
vpath %.F90 source
$(BUILDPATH)/%.p.F90.up : source/%.F90 $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml
	./scripts/build/preprocess.pl source/$*.F90 $(BUILDPATH)/$*.p.F90
$(BUILDPATH)/%.p.F90 : $(BUILDPATH)/%.p.F90.up
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

# Rule for building include file with preprocessor directives to detect OS. For some reason gfortran does not define these automatically.
$(BUILDPATH)/os.inc:
	$(CCOMPILER) -dM -E - < /dev/null | grep -e __APPLE__ -e __linux__ -e __aarch64__ > $(BUILDPATH)/os.inc

# Rules for building HDF5 C interoperability types data file.
$(BUILDPATH)/hdf5FCInterop.dat  : $(BUILDPATH)/hdf5FCInterop.exe $(BUILDPATH)/hdf5FCInteropC.exe
	$(BUILDPATH)/hdf5FCInterop.exe  >  $(BUILDPATH)/hdf5FCInterop.dat
	$(BUILDPATH)/hdf5FCInteropC.exe >> $(BUILDPATH)/hdf5FCInterop.dat
$(BUILDPATH)/hdf5FCInterop.exe  : source/hdf5FCInterop.F90
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) source/hdf5FCInterop.F90 -o $(BUILDPATH)/hdf5FCInterop.exe $(FCFLAGS)
$(BUILDPATH)/hdf5FCInteropC.exe : source/hdf5FCInteropC.c
	$(CCOMPILER) source/hdf5FCInteropC.c -o $(BUILDPATH)/hdf5FCInteropC.exe $(CFLAGS)

# Configuration of proc filesystem.
-include $(BUILDPATH)/Makefile_Config_Proc
$(BUILDPATH)/Makefile_Config_Proc: source/proc_config.c
	@mkdir -p $(BUILDPATH)
	@touch $(BUILDPATH)/Makefile_Config_Proc
	$(CCOMPILER) source/proc_config.c -o $(BUILDPATH)/proc_config $(CFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 $(BUILDPATH)/proc_config > /dev/null 2>&1 ; \
	 if [ $$? -eq 0 ] ; then \
	  echo "FCFLAGS  += -DPROCPS"   >  $(BUILDPATH)/Makefile_Config_Proc ; \
	  echo "CFLAGS   += -DPROCPS"   >> $(BUILDPATH)/Makefile_Config_Proc ; \
	  echo "CPPFLAGS += -DPROCPS"   >> $(BUILDPATH)/Makefile_Config_Proc ; \
	 fi \
	fi

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

# Configuration for availability of qhull.
-include $(BUILDPATH)/Makefile_Config_QHull
$(BUILDPATH)/Makefile_Config_QHull: source/qhull_config.cpp
	@mkdir -p $(BUILDPATH)
	$(CPPCOMPILER) -c source/qhull_config.cpp -o $(BUILDPATH)/qhull_config.o $(CPPFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DQHULLAVAIL"   >  $(BUILDPATH)/Makefile_Config_QHull ; \
	 echo "CPPFLAGS += -DQHULLAVAIL"   >> $(BUILDPATH)/Makefile_Config_QHull ; \
	else \
	 echo "FCFLAGS  += -DQHULLUNAVAIL" >  $(BUILDPATH)/Makefile_Config_QHull ; \
	 echo "CPPFLAGS += -DQHULLUNAVAIL" >> $(BUILDPATH)/Makefile_Config_QHull ; \
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

# Configuration for availability of libgit2. For static builds, we
# make libgit2 unavailable, as typically the `gssapi_krb5` library
# (linked by libgit2) does not have a static version available. In
# this case we will fall back to working through the `git` command
# line.
-include $(BUILDPATH)/Makefile_Config_Git2
ifeq '${STATIC}' '-static'
$(BUILDPATH)/Makefile_Config_Git2:
	@mkdir -p $(BUILDPATH)
	echo "FCFLAGS  += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2
	echo "CFLAGS   += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2
	echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
else ifeq '${USEGIT2}' 'no'
$(BUILDPATH)/Makefile_Config_Git2:
	@mkdir -p $(BUILDPATH)
	echo "FCFLAGS  += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2
	echo "CFLAGS   += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2
	echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
else
$(BUILDPATH)/Makefile_Config_Git2: source/libgit2_config.c
	@mkdir -p $(BUILDPATH)
	$(CCOMPILER) -c source/libgit2_config.c -o $(BUILDPATH)/libgit2_config.o $(CFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DGIT2AVAIL"   >  $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CFLAGS   += -DGIT2AVAIL"   >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CPPFLAGS += -DGIT2AVAIL"   >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	else \
	 echo "FCFLAGS  += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CFLAGS   += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	fi
endif

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
	if command -v wget &> /dev/null; then \
	 wget --no-check-certificate https://github.com/emsig/fftlog/archive/refs/heads/main.zip -O  source/FFTlog/main.zip; \
	else \
	 curl --insecure -L https://github.com/emsig/fftlog/archive/refs/heads/main.zip --output source/FFTlog/main.zip;\
	fi
	cd source/FFTlog; \
	unzip main.zip; \
	mv fftlog-main/src/*.f .; \
	cd -
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

# Object (*.o) files are built by compiling Fortran (*.[fF]) source files.
vpath %.f source
$(BUILDPATH)/%.o : %.f $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(F77FLAGS)
vpath %.F source
$(BUILDPATH)/%.o : %.F $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
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
$(BUILDPATH)/%.Inc : $(BUILDPATH)/%.Inc.up
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
$(BUILDPATH)/%.d : ./source/%.F
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
%.d : %.F
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
$(BUILDPATH)/%.fl : ./source/%.F
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
%.exe: $(BUILDPATH)/%.o $(BUILDPATH)/%.d `cat $(BUILDPATH)/$*.d` $(MAKE_DEPS)
	./scripts/build/parameterDependencies.pl `pwd` $*.exe
	$(FCCOMPILER) -c $(BUILDPATH)/$*.parameters.F90 -o $(BUILDPATH)/$*.parameters.o $(FCFLAGS)
	./scripts/build/sourceDigests.pl `pwd` $*.exe
	$(CCOMPILER) -c $(BUILDPATH)/$*.md5s.c -o $(BUILDPATH)/$*.md5s.o $(CFLAGS)
	$(CONDORLINKER) $(FCCOMPILER) `cat $*.d` $(BUILDPATH)/$*.parameters.o $(BUILDPATH)/$*.md5s.o -o $*.exe$(SUFFIX) $(FCFLAGS) `scripts/build/libraryDependencies.pl $*.exe $(FCFLAGS)` 2>&1 | ./scripts/build/postprocessLinker.pl

# Library.
-include $(BUILDPATH)/Makefile_Library_Dependencies 
$(BUILDPATH)/Makefile_Library_Dependencies:
	./scripts/build/libraryInterfacesDependencies.pl
$(BUILDPATH)/libgalacticus.Inc: $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/stateStorables.xml
	./scripts/build/libraryInterfaces.pl
$(BUILDPATH)/libgalacticus.p.Inc.up : $(BUILDPATH)/libgalacticus.Inc $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml
	./scripts/build/preprocess.pl $(BUILDPATH)/libgalacticus.Inc $(BUILDPATH)/libgalacticus.p.Inc
$(BUILDPATH)/libgalacticus.p.Inc : $(BUILDPATH)/libgalacticus.p.Inc.up
	@true
$(BUILDPATH)/libgalacticus.inc : $(BUILDPATH)/libgalacticus.p.Inc Makefile
	perl -MRegexp::Common -ne '$$l=$$_;$$l =~ s/($$RE{comment}{Fortran}{-keep})/\/\*$$4\*\/$$5/; print $$l' $(BUILDPATH)/libgalacticus.p.Inc | cpp -nostdinc -C | perl -MRegexp::Common -ne '$$l=$$_;$$l =~ s/($$RE{comment}{C}{-keep})/!$$4/; print $$l' > $(BUILDPATH)/libgalacticus.tmp
	mv -f $(BUILDPATH)/libgalacticus.tmp $(BUILDPATH)/libgalacticus.inc
libgalacticus.so: $(BUILDPATH)/libgalacticus.o $(BUILDPATH)/libgalacticus_classes.d
	./scripts/build/parameterDependencies.pl `pwd` libgalacticus.o
	$(FCCOMPILER) -c $(BUILDPATH)/libgalacticus.parameters.F90 -o $(BUILDPATH)/libgalacticus.parameters.o $(FCFLAGS)
	./scripts/build/sourceDigests.pl `pwd` libgalacticus.o
	$(CCOMPILER) -c $(BUILDPATH)/libgalacticus.md5s.c -o $(BUILDPATH)/libgalacticus.md5s.o $(CFLAGS)
	$(FCCOMPILER) -shared `sort -u $(BUILDPATH)/libgalacticus.d $(BUILDPATH)/libgalacticus_classes.d` $(BUILDPATH)/libgalacticus.parameters.o $(BUILDPATH)/libgalacticus.md5s.o -o libgalacticus.so $(FCFLAGS) `scripts/build/libraryDependencies.pl libgalacticus.o $(FCFLAGS)`

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

$(BUILDPATH)/openMPCriticalSections.count.inc $(BUILDPATH)/openMPCriticalSections.enumerate.inc: $(BUILDPATH)/openMPCriticalSections.xml
	@touch $(BUILDPATH)/openMPCriticalSections.count.inc $(BUILDPATH)/openMPCriticalSections.enumerate.inc
$(BUILDPATH)/openMPCriticalSections.xml: ./scripts/build/enumerateOpenMPCriticalSections.pl
	./scripts/build/enumerateOpenMPCriticalSections.pl `pwd`

# Dependency on dependencies.
$(BUILDPATH)/utility.dependencies.p.F90.up : aux/dependencies.yml

# Rules for version routines.
$(BUILDPATH)/output.version.revision.inc: $(wildcard .git/refs/heads/master)
	@if [ -f .git/refs/heads/master ] ; then git rev-parse HEAD | awk '{print "character(len=40), parameter :: gitHash=\""$$1"\""}' > $(BUILDPATH)/output.version.revision.inc; else printf 'character(len=40), parameter :: gitHash="unknown"\n' > $(BUILDPATH)/output.version.revision.inc; fi
	@if [ -f .git/refs/heads/master ] ; then git branch | awk '{if ($$1 == "*") print "character(len=128), parameter :: gitBranch=\""$$2"\""}' >> $(BUILDPATH)/output.version.revision.inc; else printf 'character(len=128), parameter :: gitBranch="(unknown)"\n' >> $(BUILDPATH)/output.version.revision.inc; fi
	@date -u '+%a %b %d %k:%M:%S UTC %Y' | awk '{print "character(len=32), parameter :: buildTime=\""$$0"\""}' >> $(BUILDPATH)/output.version.revision.inc

# Rules for build information routines.
$(BUILDPATH)/output.build.environment.inc:
	@echo PREPROCESSOR=\"$(PREPROCESSOR)\" > $(BUILDPATH)/output.build.environment.inc
	@echo FCCOMPILER=\"$(FCCOMPILER)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CCOMPILER=\"$(CCOMPILER)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CPPCOMPILER=\"$(CPPCOMPILER)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo FCFLAGS=\"$(FCFLAGS)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo FCFLAGS_NOOPT=\"$(FCFLAGS_NOOPT)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CFLAGS=\"$(CFLAGS)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CPPFLAGS=\"$(CPPFLAGS)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo LIBS=\"$(LIBS)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo FCCOMPILER_VERSION=\"$(FCCOMPILER_VERSION)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CCOMPILER_VERSION=\"$(CCOMPILER_VERSION)\" >> $(BUILDPATH)/output.build.environment.inc
	@echo CPPCOMPILER_VERSION=\"$(CPPCOMPILER_VERSION)\" >> $(BUILDPATH)/output.build.environment.inc

# Rules for changeset creation.
Galacticus.exe: $(BUILDPATH)/galacticus.git.patch $(BUILDPATH)/galacticus.git.bundle
$(BUILDPATH)/galacticus.git.patch:
	git diff > $(BUILDPATH)/galacticus.git.patch 2>&1 || echo unknown > $(BUILDPATH)/galacticus.git.patch
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
$(BUILDPATH)/Makefile_Module_Dependencies: ./scripts/build/moduleDependencies.pl $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies $(ALLSOURCESINC)
	@mkdir -p $(BUILDPATH)
	./scripts/build/moduleDependencies.pl `pwd`

$(BUILDPATH)/Makefile_Use_Dependencies: ./scripts/build/useDependencies.pl $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies $(BUILDPATH)/Makefile_Library_Dependencies $(ALLSOURCESINC)
	@mkdir -p $(BUILDPATH)
	./scripts/build/useDependencies.pl `pwd`

$(BUILDPATH)/Makefile_Directives: ./scripts/build/codeDirectivesParse.pl $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/codeDirectivesParse.pl `pwd`
	./scripts/build/stateStorables.pl `pwd`
	./scripts/build/deepCopyActions.pl `pwd`

$(BUILDPATH)/Makefile_Include_Dependencies: ./scripts/build/includeDependencies.pl $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/includeDependencies.pl `pwd`

$(BUILDPATH)/Makefile_All_Execs: ./scripts/build/findExecutables.pl $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/findExecutables.pl `pwd`

deps: $(MAKE_DEPS) $(BUILDPATH)/Makefile_All_Execs

# Rules for XSpec code.
aux/XSpec/%.o: ./aux/XSpec/%.f Makefile
	$(FCCOMPILER) -c $< -o aux/XSpec/$*.o $(FCFLAGS)
