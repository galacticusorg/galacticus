# Galacticus Singularity image
# Uses SingularityHub to build Galacticus.
# Version: 2020-08-06

Bootstrap:docker
From:galacticusorg/buildenv:latest

%environment
	export INSTALL_PATH=/usr/local
	export PATH=$INSTALL_PATH/gcc-10/bin:$INSTALL_PATH/bin:$PATH
	export LD_LIBRARY_PATH=$INSTALL_PATH/lib64:$INSTALL_PATH/lib:$INSTALL_PATH/gcc-10/lib64:$INSTALL_PATH/gcc-10/lib:/usr/lib/x86_64-linux-gnu
	export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
	export GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH/finclude -fintrinsic-modules-path $INSTALL_PATH/include -fintrinsic-modules-path $INSTALL_PATH/include/gfortran -fintrinsic-modules-path $INSTALL_PATH/lib/gfortran/modules -L$INSTALL_PATH/lib -L$INSTALL_PATH/lib64 -fuse-ld=bfd"
	export GALACTICUS_CFLAGS="-fuse-ld=bfd"
	export GALACTICUS_CPPFLAGS="-fuse-ld=bfd"
	export GALACTICUS_EXEC_PATH=/opt/galacticus
	export GALACTICUS_DATA_PATH=/opt/datasets

%runscript
	echo "Building Galacticus container..."

%post
	export INSTALL_PATH=/usr/local
	export PATH=$INSTALL_PATH/gcc-10/bin:$INSTALL_PATH/bin:$PATH
	export LD_LIBRARY_PATH=$INSTALL_PATH/lib64:$INSTALL_PATH/lib:$INSTALL_PATH/gcc-10/lib64:$INSTALL_PATH/gcc-10/lib:/usr/lib/x86_64-linux-gnu
	export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
	export GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH/finclude -fintrinsic-modules-path $INSTALL_PATH/include -fintrinsic-modules-path $INSTALL_PATH/include/gfortran -fintrinsic-modules-path $INSTALL_PATH/lib/gfortran/modules -L$INSTALL_PATH/lib -L$INSTALL_PATH/lib64 -fuse-ld=bfd"
	export GALACTICUS_CFLAGS="-fuse-ld=bfd"
	export GALACTICUS_CPPFLAGS="-fuse-ld=bfd"
	export GALACTICUS_EXEC_PATH=/opt/galacticus
	export GALACTICUS_DATA_PATH=/opt/datasets
	echo ENVIRONMENT
	env
	echo $GALACTICUS_EXEC_PATH
	cd /opt
	git clone --depth 1 https://github.com/galacticusorg/galacticus.git galacticus
	git clone --depth 1 https://github.com/galacticusorg/datasets.git datasets
	cd /opt/galacticus
	make -j2 Galacticus.exe
	./Galacticus.exe parameters/buildTools.xml
	rm /opt/datasets/dynamic/c17.02.tar.gz /opt/datasets/dynamic/CAMB.tar.gz
