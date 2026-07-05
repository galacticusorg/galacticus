# Main Makefile for building Galacticus.
#
# Andrew Benson (06-Feb-2010)

# Detect Operating System
UNAME_S := $(shell uname -s)

# Conditional assignment for macOS
ifeq ($(UNAME_S),Darwin)
    # For MacOS we must set LC_ALL=C to avoid problems with non-UTF8 characters and sed.
    export LC_ALL=C
endif

# Build option.
GALACTICUS_BUILD_OPTION ?= default
ifdef  BUILDPATH
 override BUILDPATH := $(patsubst %/,%,$(BUILDPATH))
endif
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
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'perf'
export BUILDPATH ?= ./work/buildPerf
export SUFFIX ?= _perf
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'odeprof'
export BUILDPATH ?= ./work/buildODEProf
export SUFFIX ?= _odeProf
else ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
export BUILDPATH ?= ./work/build
export SUFFIX ?=
endif

# Convenience flag: non-empty when this is a shared-library build. Use via `ifneq ($(IS_LIB_BUILD),)` to gate
# library-only logic (PIC flags, library interface generation, dependency prereqs).
IS_LIB_BUILD := $(filter lib,$(GALACTICUS_BUILD_OPTION))

# Preprocessor:
PREPROCESSOR ?= cpp

# Make the Galacticus python tree importable for any python3 process the
# Makefile launches.  This replaces the per-module
# `sys.path.insert(0, $GALACTICUS_EXEC_PATH/python)` shims that the
# scripts under scripts/build/ used to carry, and complements (does not
# require) `pip install -e .` -- either is sufficient on its own.
export PYTHONPATH := $(CURDIR)/python$(if $(PYTHONPATH),:$(PYTHONPATH))

# Profiling options.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'compileprof'
SHELL = ./scripts/build/profiler.sh
endif

# Fortran compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
ifdef MPIFCCOMPILER
FCCOMPILER = $(MPIFCCOMPILER)
else
FCCOMPILER ?= mpif90
endif
else
FCCOMPILER ?= gfortran
endif

# C compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
ifdef MPICCOMPILER
CCOMPILER = $(MPICCOMPILER)
else
CCOMPILER ?= mpicc
endif
else
CCOMPILER ?= gcc
endif
export CCOMPILER

# C++ compiler:
ifeq '$(GALACTICUS_BUILD_OPTION)' 'MPI'
ifdef MPICPPCOMPILER
CPPCOMPILER = $(MPICPPCOMPILER)
else
CPPCOMPILER ?= mpic++
endif
else
CPPCOMPILER ?= g++
endif

# Fortran compiler flags:
FCFLAGS += -ffree-line-length-none -frecursive -DBUILDPATH=\'$(BUILDPATH)\' -J$(BUILDPATH)/moduleBuild/ -I$(BUILDPATH)/ ${GALACTICUS_FCFLAGS} -pthread
# Fortran77 compiler flags:
F77FLAGS = ${GALACTICUS_F77FLAGS} -DBUILDPATH=\'$(BUILDPATH)\'
# Error checking flags
FCFLAGS += -Wall -fbacktrace -ffpe-trap=invalid,zero,overflow -fdump-core
# Add bounds checking.
#FCFLAGS += -fbounds-check
# A copy of the flags prior to any optimizations. No longer used to compile anything: it is
# recorded in output.build.environment.inc and reported as build provenance metadata by
# source/output/build.F90.
FCFLAGS_NOOPT := $(FCFLAGS)
# Optimization flags.
FCFLAGS += -O3 -ffinite-math-only -fno-math-errno
# For OpenMP compilation.
FCFLAGS  += -fopenmp
CFLAGS   += -fopenmp
CPPFLAGS += -fopenmp
# Link-time optimization option. Enabled by default; set `LTO=disabled` to turn off. This is needed, for
# example, on Apple Silicon, where the DWARF that LTO emits into the object files is large and only loosely
# matches the final, link-time recompiled code, which makes `dsymutil` balloon in memory and get OOM-killed
# (even `dsymutil --num-threads 1` is insufficient).
LTO ?= enabled
ifeq '$(LTO)' 'enabled'
# Use `-flto=jobserver` (not `-flto=auto`) so the parallelism of the link-time recompilation (the LTRANS
# phase) is governed by GNU make's job server rather than by the number of CPUs the compiler detects. This
# fixes two problems on shared/HPC nodes:
#   1. `-flto=auto` detects *every* CPU on the node (e.g. 640) regardless of how many were requested via
#      `make -jN` (e.g. 64), badly oversubscribing a node we only partially own.
#   2. With `make -jN all`, many executables link simultaneously; with `-flto=auto` each link independently
#      spawns N LTRANS jobs, so the total process count explodes (N links x N jobs).
# The job server is a single token pool shared across the whole build, so total LTRANS parallelism is capped
# at the `-jN` requested, no matter how many executables link at once. This requires a `+` prefix on the link
# recipe (see the `%.exe` rule) so make exposes the job server to the compiler. With no job server present
# (e.g. a plain `make Galacticus.exe` with no `-j`), LTO simply runs serially.
FCFLAGS  += -flto=jobserver
CFLAGS   += -flto=jobserver
CPPFLAGS += -flto=jobserver
endif
# Detect static compilation
STATIC=$(findstring -static,${FCFLAGS})
ifeq '${STATIC}' '-static'
FCFLAGS  += -DSTATIC
CFLAGS   += -DSTATIC
CPPFLAGS += -DSTATIC
endif

# The -O3 + LTO middle-end emits false-positive -Wstringop-overread warnings on in-place Fortran
# character substring assignments (e.g. `s = s(2:len_trim(s))`). This flag is rejected by the Fortran
# front-end (f951), so it must be applied only at the LTO link step via FCFLAGS_LINK (see the %.exe rule).
FCFLAGS_LINK  += -Wno-stringop-overread

# C compiler flags. The source tree is hierarchical, so add an include path for every source
# subdirectory (header files such as md5.h or gsl_odeiv2.h live in subdirectories).
CFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' $(addprefix -I,$(SOURCEDIRS)) -I$(BUILDPATH)/ ${GALACTICUS_CFLAGS}
export CFLAGS

# C++ compiler flags:
CPPFLAGS += -DBUILDPATH=\'$(BUILDPATH)\' $(addprefix -I,$(SOURCEDIRS)) -I$(BUILDPATH)/ ${GALACTICUS_CPPFLAGS}

# Detect library compile.
ifneq ($(IS_LIB_BUILD),)
FCFLAGS       += -fPIC
FCFLAGS_NOOPT += -fPIC
F77FLAGS      += -fPIC
CFLAGS        += -fPIC
CPPFLAGS      += -fPIC
# Place GNU Fortran trampolines (used when internal procedures are passed as actual arguments) on the heap
# rather than the stack. Stack-based trampolines require an executable stack, which causes the resulting
# shared library to be marked as such; newer kernels/loaders then refuse to `dlopen()` it, breaking the
# Python interface with "cannot enable executable stack as shared object requires: Invalid argument". Heap
# trampolines avoid this without an executable stack. This flag requires GCC 14+ (Galacticus requires GCC
# 16+) and is applied only to library builds.
FCFLAGS       += -ftrampoline-impl=heap
FCFLAGS_NOOPT += -ftrampoline-impl=heap
# The `-z noexecstack` linker option (used below when linking libgalacticus.so) is specific to GNU ld.
# Apple's linker does not understand it, and Mach-O has no executable-stack marking to begin with, so set
# it only on non-Darwin (Linux) systems.
ifneq ($(UNAME_S),Darwin)
LINKNOEXECSTACK := -Wl,-z,noexecstack
endif
endif

# Add debugging symbols.
FCFLAGS       += -g
FCFLAGS_NOOPT += -g
F77FLAGS      += -g
CFLAGS        += -g
CPPFLAGS      += -g

# Detect GProf compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'gprof'
FCFLAGS       += -pg
FCFLAGS_NOOPT += -pg
F77FLAGS      += -pg
CFLAGS        += -pg
CPPFLAGS      += -pg
endif

# Detect perf compile.
ifeq '$(GALACTICUS_BUILD_OPTION)' 'perf'
FCFLAGS       += -fno-omit-frame-pointer
FCFLAGS_NOOPT += -fno-omit-frame-pointer
F77FLAGS      += -fno-omit-frame-pointer
CFLAGS        += -fno-omit-frame-pointer
CPPFLAGS      += -fno-omit-frame-pointer
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

# List of additional Makefiles which contain dependency information. The library interface dependencies are only needed for
# library builds, so are added below conditionally - generating them is slow and unnecessary for regular (or MPI) builds.
MAKE_DEPS = $(BUILDPATH)/Makefile_Module_Dependencies $(BUILDPATH)/Makefile_Use_Dependencies $(BUILDPATH)/Makefile_Include_Dependencies
ifneq ($(IS_LIB_BUILD),)
MAKE_DEPS += $(BUILDPATH)/Makefile_Library_Dependencies
endif

# Get versions of build tools.
FCCOMPILER_VERSION = `$(FCCOMPILER) -v 2>&1`
CCOMPILER_VERSION = `$(CCOMPILER) -v 2>&1`
CPPCOMPILER_VERSION = `$(CPPCOMPILER) -v 2>&1`

# Determine if MD5 locking is needed.
COUNT_TARGETS := $(words $(MAKECMDGOALS))
LOCKMD5 ?= yes
ifeq ($(COUNT_TARGETS),1)
  LOCKMD5=no
endif

# Recursively expand a wildcard using only Make built-in functions. This deliberately avoids
# $(shell ...): the 'compileprof' build option sets SHELL to a profiler wrapper that prints timing
# text (containing colons) to stdout, which would corrupt any $(shell ...)-derived file lists and
# break Makefile parsing.
#   $(call rwildcard,<dir>,<pattern>)  e.g. $(call rwildcard,source,*.F90)
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2)$(filter $(subst *,%,$2),$d))

# Find all source files, recursing through the full source directory hierarchy.
ALLSOURCES    = $(call rwildcard,source,*.F90) $(call rwildcard,source,*.f90) $(call rwildcard,source,*.h) $(call rwildcard,source,*.c) $(call rwildcard,source,*.cpp)
ALLSOURCESINC = $(ALLSOURCES) $(call rwildcard,source,*.Inc)

# All source directories (the source tree is hierarchical), for use in vpath and include search
# paths. Enumerate *every* directory under source/ -- not just those that directly contain source
# files -- because some sources include a vendored header by a parent-relative path (e.g.
# `#include <gslODEInitVal2/gsl_odeiv2.h>`), which needs -I pointing at an intermediate directory
# (here source/external/) that contains no source files of its own. `$(wildcard $d/.)` is non-empty
# only for directories, so this recurses through the tree using only Make built-ins (no $(shell),
# which would be corrupted by the 'compileprof' profiler SHELL).
rsubdirs = $(foreach d,$(wildcard $1/*),$(if $(wildcard $d/.),$d $(call rsubdirs,$d)))
SOURCEDIRS := source $(call rsubdirs,source)

# General suffix rules: i.e. rules for making a file of one suffix from files of another suffix.

# Object (*.o) files are built by preprocessing and then compiling Fortran 90 (*.F90) source
# files. Note that .F90 source files should not have names which coincide with the name of a
# module - this will lead to circular dependency problems as Make becomes confused about how to
# build the module file.
vpath %.F90 $(SOURCEDIRS)
# Preprocessing depends on stateStorables.xml and deepCopyActions.xml because the preprocessor's
# process hooks (StateStorable, DeepCopyActions, FunctionClass) read them at run time — generated
# state-store and deep-copy code embeds knowledge of every cataloged type. Both catalogs are
# written only-if-changed, so this triggers a re-preprocess sweep only when catalog content
# actually changes (e.g. a stateStorable type is added or removed), and the `.up` sentinel then
# limits recompilation to files whose preprocessed output really differs.
$(BUILDPATH)/%.p.F90.up : source/%.F90 $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml $(BUILDPATH)/stateStorables.xml $(BUILDPATH)/deepCopyActions.xml
	./scripts/build/preprocess.py source/$*.F90 $(BUILDPATH)/$*.p.F90
$(BUILDPATH)/%.p.F90 : $(BUILDPATH)/%.p.F90.up
	@true
# Determine whether we are compiling for Apple Silicon (macOS on AArch64). We query the C compiler's
# predefined macros (as is done for the os.inc rule below) rather than the host's uname, so that this
# reflects the actual compilation target; note that gfortran does not define these macros itself,
# which is why the C compiler is used.
APPLE_SILICON := $(shell defs=`$(CCOMPILER) -dM -E - < /dev/null 2>/dev/null`; echo "$$defs" | grep -q __APPLE__ && echo "$$defs" | grep -q __aarch64__ && echo yes)
# Work around a gfortran (GCC) internal compiler error in the AArch64 (Apple Silicon) back-end. With
# gfortran 16 the constructors of some output analysis classes trigger:
#   internal compiler error: in aarch64_function_arg_alignment, at config/aarch64/aarch64.cc
# during the RTL "expand" pass, while laying out procedure arguments. This is a compiler bug (not an
# error in our code) - see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=124146. The ICE is present at
# -O1 and above, so we compile the affected object files at -O0 to avoid it (this flag is appended
# after the global -O3, and so takes precedence). This is gated on Apple Silicon so that other
# architectures retain full optimization. Add further object files to this list if they are found to
# trigger the same ICE; these overrides (and the gate) can be removed once the upstream bug is fixed.
ifeq ($(APPLE_SILICON),yes)
FCFLAGS_AARCH64_ICE_OBJECTS = \
	$(BUILDPATH)/output.analyses.volume_function_1d.o
$(FCFLAGS_AARCH64_ICE_OBJECTS): FCFLAGS += -O0
endif

# The compiler writes diagnostics to a file which is then fed to postprocess.py, rather than piping
# compiler output directly into postprocess.py. In a pipeline the compiler's exit status would be
# discarded (the recipe's status would be postprocess.py's, and /bin/sh has no pipefail), so a hard
# compiler failure that emits no "Error:" line — an internal compiler error, a segfault, an OOM
# kill — would be reported as success, leaving a missing or stale object file that only surfaces
# later as a confusing link error. Capturing and testing both statuses makes such failures fail the
# recipe immediately. (libraryInterfacesDependencies.py emits a mirror of this recipe for the
# library wrapper units — keep the two in sync.)
$(BUILDPATH)/%.o : $(BUILDPATH)/%.p.F90 $(BUILDPATH)/%.m $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $(BUILDPATH)/$*.p.F90 -o $(BUILDPATH)/$*.o $(FCFLAGS) > $(BUILDPATH)/$*.o.diag 2>&1; \
	compileStatus=$$?; \
	./scripts/build/postprocess.py $(BUILDPATH)/$*.p.F90 < $(BUILDPATH)/$*.o.diag; \
	postprocessStatus=$$?; \
	rm -f $(BUILDPATH)/$*.o.diag; \
	[ $$compileStatus -eq 0 ] && [ $$postprocessStatus -eq 0 ]
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
$(BUILDPATH)/hdf5FCInterop.exe  : source/system/hdf5FCInterop.F90
	@mkdir -p $(BUILDPATH)/moduleBuild
	+$(FCCOMPILER) source/system/hdf5FCInterop.F90 -o $(BUILDPATH)/hdf5FCInterop.exe $(FCFLAGS)
$(BUILDPATH)/hdf5FCInteropC.exe : source/system/hdf5FCInteropC.c
	+$(CCOMPILER) source/system/hdf5FCInteropC.c -o $(BUILDPATH)/hdf5FCInteropC.exe $(CFLAGS)

# Canned recipe for feature-availability probes:
#   $(call CONFIG_PROBE,<name>,<compiler-var>,<source>,<flags-var>,<flag-vars>,<symbol>)
# Attempts to compile the probe source $(3) with the compiler named by $(2) (using the flags
# variable named by $(4)), then writes $(BUILDPATH)/Makefile_Config_$(1) appending -D$(6)AVAIL (on
# success) or -D$(6)UNAVAIL (on failure) to each flags variable listed in $(5). Compiler and flags
# are passed by NAME (dereferenced via the computed variable reference `$($(2))`) because flag
# values may contain commas, which would break $(call) argument splitting. To add a new probe, add
# a `-include` line, a rule with the probe source as prerequisite, and one $(call) line — see the
# probes below. Probes with a different shape (Proc runs the compiled binary; Git2 has static/opt-
# out variants and a restricted include path) remain hand-written.
define CONFIG_PROBE
@mkdir -p $(BUILDPATH) $(BUILDPATH)/moduleBuild
$($(2)) -c $(3) -o $(BUILDPATH)/$(1)_config.o $($(4)) > /dev/null 2>&1 ; \
if [ $$? -eq 0 ] ; then availability=AVAIL ; else availability=UNAVAIL ; fi ; \
rm -f $(BUILDPATH)/Makefile_Config_$(1) ; \
for flagVariable in $(5) ; do \
 echo "$$flagVariable += -D$(6)$$availability" >> $(BUILDPATH)/Makefile_Config_$(1) ; \
done
endef

# Configuration of proc filesystem.
-include $(BUILDPATH)/Makefile_Config_Proc
$(BUILDPATH)/Makefile_Config_Proc: source/system/proc_config.c
	@mkdir -p $(BUILDPATH)
	@touch $(BUILDPATH)/Makefile_Config_Proc
	$(CCOMPILER) source/system/proc_config.c -o $(BUILDPATH)/proc_config $(CFLAGS) > /dev/null 2>&1 ; \
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
$(BUILDPATH)/Makefile_Config_OFD: source/system/flock_config.c
	$(call CONFIG_PROBE,OFD,CCOMPILER,source/system/flock_config.c,CFLAGS,FCFLAGS CFLAGS CPPFLAGS,OFD)

# Configuration for availability of FFTW3.
-include $(BUILDPATH)/Makefile_Config_FFTW3
$(BUILDPATH)/Makefile_Config_FFTW3: source/external/FFTW/fftw3_config.F90
	$(call CONFIG_PROBE,FFTW3,FCCOMPILER,source/external/FFTW/fftw3_config.F90,FCFLAGS,FCFLAGS,FFTW3)

# The ANN, qhull, and libmatheval probes compile against only the system/user headers
# (GALACTICUS_CPPFLAGS), *not* the full CPPFLAGS. The latter adds `-I$(BUILDPATH)/`, which can
# contain zero-byte header stubs left by an earlier UNAVAIL build (created by the generic `%.h`
# rule for headers that were preprocessed out — e.g. `Qhull.h`). Such a stub satisfies the probe's
# `#include`, making it spuriously report AVAIL even though the real headers are absent, and the
# subsequent real compile then fails (or the availability flip-flops between builds). This mirrors
# the libgit2 probe below, which hit the identical trap first. None of these probe sources need
# vendored or generated headers, so the restricted include path loses nothing.

# Configuration for availability of ANN.
-include $(BUILDPATH)/Makefile_Config_ANN
$(BUILDPATH)/Makefile_Config_ANN: source/external/ANN/ann_config.cpp
	$(call CONFIG_PROBE,ANN,CPPCOMPILER,source/external/ANN/ann_config.cpp,GALACTICUS_CPPFLAGS,FCFLAGS CPPFLAGS,ANN)

# Configuration for availability of qhull.
-include $(BUILDPATH)/Makefile_Config_QHull
$(BUILDPATH)/Makefile_Config_QHull: source/external/Qhull/qhull_config.cpp
	$(call CONFIG_PROBE,QHull,CPPCOMPILER,source/external/Qhull/qhull_config.cpp,GALACTICUS_CPPFLAGS,FCFLAGS CPPFLAGS,QHULL)

# Configuration for availability of libmatheval.
-include $(BUILDPATH)/Makefile_Config_MathEval
$(BUILDPATH)/Makefile_Config_MathEval: source/system/libmatheval_config.cpp
	$(call CONFIG_PROBE,MathEval,CPPCOMPILER,source/system/libmatheval_config.cpp,GALACTICUS_CPPFLAGS,FCFLAGS CPPFLAGS,MATHEVAL)

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
	echo "CFLAGS   += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
	echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
else ifeq '${USEGIT2}' 'no'
$(BUILDPATH)/Makefile_Config_Git2:
	@mkdir -p $(BUILDPATH)
	echo "FCFLAGS  += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2
	echo "CFLAGS   += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
	echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2
else
$(BUILDPATH)/Makefile_Config_Git2: source/system/libgit2_config.c
	@mkdir -p $(BUILDPATH)
# Probe only against the system/user libgit2 headers (via GALACTICUS_CFLAGS),
# *not* the full CFLAGS. The latter adds `-I$(BUILDPATH)/`, which can contain a
# zero-byte `git2.h` stub left by an earlier GIT2UNAVAIL build (created by the
# generic `%.h` rule for the preprocessed-out include in git2.c). That stub
# would shadow the real header and make the probe spuriously fail, trapping the
# build in GIT2UNAVAIL even when a working libgit2 is installed.
	$(CCOMPILER) -c source/system/libgit2_config.c -o $(BUILDPATH)/libgit2_config.o $(GALACTICUS_CFLAGS) > /dev/null 2>&1 ; \
	if [ $$? -eq 0 ] ; then \
	 echo "FCFLAGS  += -DGIT2AVAIL"   >  $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CFLAGS   += -DGIT2AVAIL"   >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CPPFLAGS += -DGIT2AVAIL"   >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	else \
	 echo "FCFLAGS  += -DGIT2UNAVAIL" >  $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CFLAGS   += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	 echo "CPPFLAGS += -DGIT2UNAVAIL" >> $(BUILDPATH)/Makefile_Config_Git2 ; \
	fi
endif

# Object (*.o) files are built by compiling C (*.c) source files.
vpath %.c $(SOURCEDIRS)
$(BUILDPATH)/%.o : %.c $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CFLAGS)

# Object (*.o) can also be built from C++ source files.
vpath %.cpp $(SOURCEDIRS)
$(BUILDPATH)/%.o : %.cpp $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	$(CPPCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(CPPFLAGS)

# Rules for the QHull library. Use the C++17 standard for these files since they are not compatible with later C++ standards
# (triggering 'template-id not allowed for constructor' errors). `Makefile` is a prerequisite, matching the generic `%.o`
# rules: qhull.cpp preprocesses to an empty translation unit under -DQHULLUNAVAIL, so a stale object compiled under a
# different availability result must be rebuilt when the flags change — otherwise the link fails with an undefined
# reference to `convexHullVolumeC` (or silently links a stub).
$(BUILDPATH)/external/Qhull/qhull.o : source/external/Qhull/qhull.cpp Makefile
	@mkdir -p $(BUILDPATH)/external/Qhull
	$(CPPCOMPILER) -c source/external/Qhull/qhull.cpp -o $(BUILDPATH)/external/Qhull/qhull.o $(CPPFLAGS) -std=gnu++17

# Rules for FFTLog library. The FFTLog source is vendored directly in
# source/external/FFTlog (it was previously downloaded from GitHub at build
# time, which made builds dependent on GitHub being reachable). See
# source/external/FFTlog/README.md for provenance, license (CC0), and a
# description of the patches that have been applied to the vendored source.
# Dependency (*.d) files for these object files are generated by the rule below,
# which mirrors the generic "$(BUILDPATH)/%.d : ./source/%.f" rule but also
# ensures the build sub-directory exists.
$(BUILDPATH)/external/FFTlog/%.d : ./source/external/FFTlog/%.f
	@mkdir -p $(BUILDPATH)/external/FFTlog
	@echo $(BUILDPATH)/external/FFTlog/$*.o > $(BUILDPATH)/external/FFTlog/$*.d~
	@if cmp -s $(BUILDPATH)/external/FFTlog/$*.d $(BUILDPATH)/external/FFTlog/$*.d~ ; then \
	 rm $(BUILDPATH)/external/FFTlog/$*.d~ ; \
	else \
	 mv $(BUILDPATH)/external/FFTlog/$*.d~ $(BUILDPATH)/external/FFTlog/$*.d ; \
	fi

$(BUILDPATH)/external/FFTlog/%.o: ./source/external/FFTlog/%.f Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	@mkdir -p $(BUILDPATH)/external/FFTlog
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/external/FFTlog/$*.o $(F77FLAGS) -Wno-argument-mismatch -std=legacy

# Object (*.o) files are built by compiling Fortran (*.[fF]) source files.
vpath %.f $(SOURCEDIRS)
$(BUILDPATH)/%.o : %.f $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(F77FLAGS)
vpath %.F $(SOURCEDIRS)
$(BUILDPATH)/%.o : %.F $(BUILDPATH)/%.d $(BUILDPATH)/%.fl Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/$*.o $(F77FLAGS)

# Special rules required for building some sources (unfortunate, but necessary....)
# pfq.new.f
$(BUILDPATH)/external/pFq/pfq.new.o : ./source/external/pFq/pfq.new.f Makefile
	@mkdir -p $(BUILDPATH)/moduleBuild
	$(FCCOMPILER) -c $< -o $(BUILDPATH)/external/pFq/pfq.new.o $(FCFLAGS)

# Rule for running *.Inc files through the preprocessor. We strip out single quote characters in comment lines to avoid spurious
# complaints from the preprocessor.
#
# The preprocessed intermediate under $(BUILDPATH) is named `%.p.Inc` (matching the libgalacticus
# chain), NOT `%.Inc`: an intermediate named `%.Inc` differs from the final processed `%.inc` only
# by case, so on case-insensitive filesystems (macOS APFS) the `mv` in the `%.inc` recipe below
# overwrote the intermediate in place, and every incremental build re-preprocessed its own cpp
# output under a perpetually-bumped timestamp.
$(BUILDPATH)/%.p.Inc.up : ./source/%.Inc $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml $(BUILDPATH)/stateStorables.xml $(BUILDPATH)/deepCopyActions.xml
	./scripts/build/preprocess.py ./source/$*.Inc $(BUILDPATH)/$*.p.Inc
$(BUILDPATH)/%.p.Inc : $(BUILDPATH)/%.p.Inc.up
	@true
$(BUILDPATH)/%.inc : $(BUILDPATH)/%.p.Inc Makefile
# Guard against a missing input: `sed` is not the last command in the pipeline, so its failure
# would otherwise be masked and an empty husk of cpp boilerplate written in place of the real
# processed include. A phantom `%.p.Inc` can arise here because the `%.p.Inc <- %.p.Inc.up`
# sentinel above "makes" its target without creating a file when the `.up` came from a fallback
# rule.
	@test -f $< || { echo "Makefile: cannot preprocess $@: input $< does not exist" >&2 ; exit 1 ; }
	sed -E s/'^([[:space:]]*)!(.*)'/'\1\/\*\2\*\/'/ $< | cpp -nostdinc -C | sed -E s/'^([[:space:]]*)\/\*(.*)\*\/'/'\1!\2'/ > $(BUILDPATH)/$*.tmp
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

# Executables (*.exe) are built by linking together all of the object files (*.o) specified in the
# associated dependency (*.d) file. There is no generic `%.exe` pattern rule: findExecutables.py
# discovers every program in the source tree and writes an explicit rule for each into
# $(BUILDPATH)/Makefile_All_Execs (included below). Those generated rules invoke the canned recipes
# defined here, so the link procedure has a single owner: change it here, not in the generator.

# Canned recipe: prepare link-time metadata objects for target $(1) (the name passed to the
# metadata scripts, e.g. `tests.foo.exe` or `libgalacticus.o`) with build-artifact stem $(2) (the
# path under $(BUILDPATH), e.g. `tests/foo` or `libgalacticus`). Enumerates every input parameter
# the target can accept (compiled into $(2).parameters.o) and computes per-type source digests
# (compiled into $(2).md5s.o). MD5 sidecar locking is enabled only under parallel make (see
# LOCKMD5); the `if` chain and the sourceDigests invocation must remain one logical shell line so
# the `useLocks` shell variable survives.
define LINK_METADATA
./scripts/build/parameterDependencies.py `pwd` $(1)
$(FCCOMPILER) -c $(BUILDPATH)/$(2).parameters.F90 -o $(BUILDPATH)/$(2).parameters.o $(FCFLAGS)
@if echo "$(MAKEFLAGS)" | grep -q -E -- ' -j1( |$$)'; then \
 useLocks=no; \
elif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j( |$$)'; then \
 useLocks=$(LOCKMD5); \
elif echo "$(MAKEFLAGS)" | grep -q -E -- ' -j[0-9]+( |$$)'; then \
 useLocks=$(LOCKMD5); \
else \
 useLocks=no; \
fi; \
./scripts/build/sourceDigests.py `pwd` $(1) $$useLocks
$(CCOMPILER) -c $(BUILDPATH)/$(2).md5s.c -o $(BUILDPATH)/$(2).md5s.o $(CFLAGS)
endef

# Canned recipe: link executable $(1).exe$(SUFFIX) from build-artifact stem $(2). The link step
# writes diagnostics to a file rather than piping them directly into postprocessLinker.py: in a
# pipeline the linker's exit status would be discarded, so a failed link (`undefined reference`,
# `ld returned 1 exit status`) would be recorded by make as a successful build. Capturing and
# testing both statuses makes link failures fail the recipe immediately. The `+` prefix exposes
# make's jobserver to the compiler so `-flto=jobserver` can parallelize the LTO link.
define LINK_EXECUTABLE
$(call LINK_METADATA,$(2).exe,$(2))
+$(FCCOMPILER) `cat $(BUILDPATH)/$(2).d` $(BUILDPATH)/$(2).parameters.o $(BUILDPATH)/$(2).md5s.o -o $(1).exe$(SUFFIX) $(FCFLAGS) $(FCFLAGS_LINK) `./scripts/build/libraryDependencies.py $(2).exe $(FCFLAGS)` > $(BUILDPATH)/$(2).link.diag 2>&1; \
linkStatus=$$?; \
./scripts/build/postprocessLinker.py < $(BUILDPATH)/$(2).link.diag; \
postprocessStatus=$$?; \
rm -f $(BUILDPATH)/$(2).link.diag; \
[ $$linkStatus -eq 0 ] && [ $$postprocessStatus -eq 0 ]
endef

# Library. These rules generate Fortran interface wrappers and their dependencies for the shared library build; the generator
# scripts (libraryInterfaces.py, libraryInterfacesDependencies.py) are slow, so we only activate them when actually performing
# a library build.
ifneq ($(IS_LIB_BUILD),)
-include $(BUILDPATH)/Makefile_Library_Dependencies
$(BUILDPATH)/Makefile_Library_Dependencies: $(BUILDPATH)/libgalacticus.Inc ./scripts/build/libraryInterfacesDependencies.py
	./scripts/build/libraryInterfacesDependencies.py
$(BUILDPATH)/libgalacticus.Inc: $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/stateStorables.xml ./source/libraryClasses.xml ./scripts/build/libraryInterfaces.py ./python/LibraryInterfaces/Pipeline.py ./python/LibraryInterfaces/Emitters.py ./python/LibraryInterfaces/ArgSpec.py ./python/LibraryInterfaces/Hierarchy.py
	./scripts/build/libraryInterfaces.py
$(BUILDPATH)/libgalacticus.p.Inc.up : $(BUILDPATH)/libgalacticus.Inc $(BUILDPATH)/hdf5FCInterop.dat $(BUILDPATH)/openMPCriticalSections.xml $(BUILDPATH)/stateStorables.xml $(BUILDPATH)/deepCopyActions.xml
	./scripts/build/preprocess.py $(BUILDPATH)/libgalacticus.Inc $(BUILDPATH)/libgalacticus.p.Inc
$(BUILDPATH)/libgalacticus.p.Inc : $(BUILDPATH)/libgalacticus.p.Inc.up
	@true
# The processed include is named `libgalacticus.preprocessed.inc`, NOT `libgalacticus.inc`: the
# latter differs from the generated `libgalacticus.Inc` only by case, and on case-insensitive
# filesystems (macOS APFS) the two are the same file. That collision let a stale or placeholder
# processed include masquerade as an up-to-date `libgalacticus.Inc` (silently skipping
# libraryInterfaces.py in mixed executable-then-library build trees), and made the `mv` below
# overwrite `libgalacticus.Inc` in place on macOS, so every incremental library rebuild
# reprocessed its own cpp output under a perpetually-bumped timestamp.
$(BUILDPATH)/libgalacticus.preprocessed.inc : $(BUILDPATH)/libgalacticus.p.Inc Makefile
	sed -E s/'^([[:space:]]*)!(.*)'/'\1\/\*\2\*\/'/ $(BUILDPATH)/libgalacticus.p.Inc | cpp -nostdinc -C | sed -E s/'^([[:space:]]*)\/\*(.*)\*\/'/'\1!\2'/ > $(BUILDPATH)/libgalacticus.tmp
	mv -f $(BUILDPATH)/libgalacticus.tmp $(BUILDPATH)/libgalacticus.preprocessed.inc
libgalacticus.so: $(BUILDPATH)/libgalacticus.o $(BUILDPATH)/libgalacticus_classes.d
	$(call LINK_METADATA,libgalacticus.o,libgalacticus)
# Link with a non-executable stack (`-z noexecstack`). Without this the shared library can be marked as
# requiring an executable stack (e.g. because an input object lacks a `.note.GNU-stack` section, or because
# GNU Fortran emits stack-based trampolines). Newer kernels/loaders refuse to `dlopen()` such a library,
# causing the Python interface to fail with "cannot enable executable stack as shared object requires:
# Invalid argument". This flag is applied only to the library link, leaving executable builds unchanged.
	+$(FCCOMPILER) -shared $(LINKNOEXECSTACK) `sort -u $(BUILDPATH)/libgalacticus.d $(BUILDPATH)/libgalacticus_classes.d` $(BUILDPATH)/libgalacticus.parameters.o $(BUILDPATH)/libgalacticus.md5s.o -o libgalacticus.so $(FCFLAGS) $(FCFLAGS_LINK) `scripts/build/libraryDependencies.py libgalacticus.o $(FCFLAGS)`
else
# Non-library builds never compile source/libgalacticus.F90, but includeDependencies.py still scans
# it and so orders Makefile_Use_Dependencies after its processed include. The interface generators
# (libraryInterfaces.py etc.) are slow and deliberately skipped here, so satisfy that ordering
# prerequisite with an explicit empty placeholder. Previously no rule existed in this branch and
# make fell back to the generic `%.inc <- %.Inc <- %.Inc.up` pattern chain, which manufactured a
# husk by running `sed` on a nonexistent file with its failure masked mid-pipeline -- the same
# silent-failure shape this build system's status-propagation work exists to eliminate.
$(BUILDPATH)/libgalacticus.preprocessed.inc :
	@touch $(BUILDPATH)/libgalacticus.preprocessed.inc
endif

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
# The enumeration must be regenerated whenever a named `!$omp critical` section is added to,
# removed from, or renamed in any Fortran source — so the rule depends on all Fortran sources. The
# script touches the `.up` sentinel on every run but rewrites the .xml only when its content
# changed, so a regeneration that finds no new critical sections does not cascade into
# re-preprocessing every source file (the `%.p.F90.up` rules depend on the .xml).
$(BUILDPATH)/openMPCriticalSections.xml.up: ./scripts/build/enumerateOpenMPCriticalSections.py $(call rwildcard,source,*.F90) $(call rwildcard,source,*.f90) $(call rwildcard,source,*.F) $(call rwildcard,source,*.f)
	@mkdir -p $(BUILDPATH)
	./scripts/build/enumerateOpenMPCriticalSections.py `pwd`
$(BUILDPATH)/openMPCriticalSections.xml: $(BUILDPATH)/openMPCriticalSections.xml.up
	@true

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
$(BUILDPATH)/Makefile_Module_Dependencies: ./scripts/build/moduleDependencies.py $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies $(ALLSOURCESINC)
	@mkdir -p $(BUILDPATH)
	./scripts/build/moduleDependencies.py `pwd`

# For library builds, useDependencies.py must scan the generated library wrapper sources under $(BUILDPATH)/libgalacticus/, so
# we make it depend on the library include generation. For non-library builds we skip this, since the wrapper sources are not
# needed and the generators (libraryInterfaces.py, libraryInterfacesDependencies.py) are slow.
ifneq ($(IS_LIB_BUILD),)
USE_DEPS_LIBRARY_PREREQS = $(BUILDPATH)/Makefile_Library_Dependencies $(BUILDPATH)/libgalacticus.Inc
endif
$(BUILDPATH)/Makefile_Use_Dependencies: ./scripts/build/useDependencies.py $(BUILDPATH)/directiveLocations.xml $(BUILDPATH)/Makefile_Directives $(BUILDPATH)/Makefile_Include_Dependencies $(USE_DEPS_LIBRARY_PREREQS) $(ALLSOURCESINC)
	@mkdir -p $(BUILDPATH)
	./scripts/build/useDependencies.py `pwd`

# The three directive-catalog scripts run in a single recipe, in this order, because
# stateStorables.py and deepCopyActions.py both read the directiveLocations.xml written by
# codeDirectivesParse.py. All three are prerequisites (previously only codeDirectivesParse.py was,
# so edits to the other two did not regenerate their catalogs). The scripts' other outputs
# (directiveLocations.xml, stateStorables.xml, deepCopyActions.xml, the per-directive XML files)
# are byproducts of this recipe; they are not given rules of their own.
#
# Makefile_Directives MUST be a direct target with this recipe, and codeDirectivesParse.py MUST
# write it unconditionally (not only-if-changed). It is `-include`d, and make only RE-EXECUTES
# itself — the restart that re-reads a regenerated included makefile — when that makefile is remade
# BY ITS OWN RULE running and updating it. On a clean build the re-read Makefile_Directives is what
# tells make how to build the generated `include`-directive files (e.g. the node-component include)
# and adds the ordering line making Makefile_Use_Dependencies depend on them. An earlier
# formulation put the real work on a `directiveCatalogs.stamp` prerequisite with an `@true` recipe
# here, and wrote Makefile_Directives only-if-changed: make then never saw Makefile_Directives as
# "remade by its rule", never restarted, and so ran useDependencies.py once — before those includes
# existed — missing every `use` they contribute and producing link-time "undefined reference"
# failures for node-component data modules, etc.
$(BUILDPATH)/Makefile_Directives: ./scripts/build/codeDirectivesParse.py ./scripts/build/stateStorables.py ./scripts/build/deepCopyActions.py $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/codeDirectivesParse.py `pwd`
	./scripts/build/stateStorables.py `pwd`
	./scripts/build/deepCopyActions.py `pwd`

$(BUILDPATH)/Makefile_Include_Dependencies: ./scripts/build/includeDependencies.py $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/includeDependencies.py `pwd`

$(BUILDPATH)/Makefile_All_Execs: ./scripts/build/findExecutables.py $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/findExecutables.py `pwd`

deps: $(MAKE_DEPS) $(BUILDPATH)/Makefile_All_Execs

# Typed parameter catalog: a machine-readable description of every functionClass
# implementation's input parameters (name, inferred type, default, nesting) and
# the objects it builds.  Regenerated whenever the source changes; consumed by
# parameter-file validation tooling.  Not on the compile path, so build it
# explicitly with `make parameters-catalog`.
$(BUILDPATH)/parameters.catalog.json: ./scripts/build/parameterCatalog.py ./python/Galacticus/Parameters/__init__.py ./python/Galacticus/Parameters/catalog.py ./python/Galacticus/Parameters/inference.py $(ALLSOURCES)
	@mkdir -p $(BUILDPATH)
	./scripts/build/parameterCatalog.py `pwd` $(BUILDPATH)/parameters.catalog.json

parameters-catalog: $(BUILDPATH)/parameters.catalog.json

# Editor-assistance schema for parameter files, generated from the catalog and
# committed (so editors can reference it). Regenerate with `make parameters-schema`
# after changing functionClass implementations or enumerations; CI checks it is
# up to date.
parameters-schema:
	./scripts/build/parameterSchema.py `pwd` schema/parameters.xsd
