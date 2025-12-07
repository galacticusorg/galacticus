# Galacticus Docker image
# Uses Docker multi-stage build to build Galacticus.

ARG TAG=latest
FROM ghcr.io/galacticusorg/buildenv:${TAG} AS build
ARG REPO=galacticusorg/galacticus
ARG BRANCH=master

# Set build options.
## * The flags are also set in galacticus/buildenv:latest so we don't really need to reset them here.
## * We force use of the BFD linker here. The GCC in galacticus/buildenv:latest uses the gold linker by default. But, the gold
##   linker seems to not correctly allow us to get values of some GSL constants (e.g. gsl_root_fsolver_brent) in Fortran.
ENV GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH/finclude -fintrinsic-modules-path $INSTALL_PATH/include -fintrinsic-modules-path $INSTALL_PATH/include/gfortran -fintrinsic-modules-path $INSTALL_PATH/lib/gfortran/modules -L$INSTALL_PATH/lib -L$INSTALL_PATH/lib64 -fuse-ld=bfd"
ENV GALACTICUS_CFLAGS="-fuse-ld=bfd"
ENV GALACTICUS_CPPFLAGS="-fuse-ld=bfd"
ENV GALACTICUS_EXEC_PATH=/opt/galacticus
ENV GALACTICUS_DATA_PATH=/opt/datasets

RUN     pwd && ls

# Clone datasets.
RUN     cd /opt &&\
	git clone --depth 1 -b ${BRANCH} https://github.com/${REPO}.git galacticus &&\
	git clone --depth 1 https://github.com/galacticusorg/datasets.git datasets

# Build Galacticus.
RUN     cd /opt/galacticus &&\
	make -j4 Galacticus.exe &&\
	rm -rf work/build

# Build external tools.
RUN     cd /opt/galacticus &&\
        classVersion=`awk '{if ($1 == "class:") print $2}' ${GALACTICUS_EXEC_PATH}/aux/dependencies.yml` &&\
        cambVersion=`awk '{if ($1 == "camb:") print $2}' ${GALACTICUS_EXEC_PATH}/aux/dependencies.yml` &&\
        forutilsVersion=`awk '{if ($1 == "forutils:") print $2}' ${GALACTICUS_EXEC_PATH}/aux/dependencies.yml` &&\
        fspsVersion=`awk '{if ($1 == "fsps:") print $2}' ${GALACTICUS_EXEC_PATH}/aux/dependencies.yml` &&\
        cloudyVersion=`awk '{if ($1 == "cloudy:") print $2}' ${GALACTICUS_EXEC_PATH}/aux/dependencies.yml` &&\
	./Galacticus.exe parameters/buildTools.xml &&\
	rm /opt/datasets/dynamic/c${cloudyVersion}.tar.gz /opt/datasets/dynamic/CAMB_${cambVersion}.tar.gz /opt/datasets/dynamic/class_public-${classVersion}.tar.gz /opt/datasets/dynamic/FSPS_${fspsVersion}.tar.gz
	
