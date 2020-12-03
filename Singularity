# Galacticus Singularity image
# Uses SingularityHub to build Galacticus.
# Version: 2020-12-03

Bootstrap:docker
From:galacticusorg/galacticus:latest

%environment
	export INSTALL_PATH=/usr/local
	export PATH=$INSTALL_PATH/gcc-11/bin:$INSTALL_PATH/bin:$PATH
	export LD_LIBRARY_PATH=$INSTALL_PATH/lib64:$INSTALL_PATH/lib:$INSTALL_PATH/gcc-11/lib64:$INSTALL_PATH/gcc-11/lib:/usr/lib/x86_64-linux-gnu
	export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
	export GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH/finclude -fintrinsic-modules-path $INSTALL_PATH/include -fintrinsic-modules-path $INSTALL_PATH/include/gfortran -fintrinsic-modules-path $INSTALL_PATH/lib/gfortran/modules -L$INSTALL_PATH/lib -L$INSTALL_PATH/lib64 -fuse-ld=bfd"
	export GALACTICUS_CFLAGS="-fuse-ld=bfd"
	export GALACTICUS_CPPFLAGS="-fuse-ld=bfd"
	export GALACTICUS_EXEC_PATH=/opt/galacticus
	export GALACTICUS_DATA_PATH=/opt/datasets

%runscript
	echo "Building Galacticus container..."

%post
	echo Begin begin: `date`
        export INSTALL_PATH=/usr/local
	export PATH=$INSTALL_PATH/gcc-11/bin:$INSTALL_PATH/bin:$PATH
	export LD_LIBRARY_PATH=$INSTALL_PATH/lib64:$INSTALL_PATH/lib:$INSTALL_PATH/gcc-11/lib64:$INSTALL_PATH/gcc-11/lib:/usr/lib/x86_64-linux-gnu
	export LIBRARY_PATH=/usr/lib/x86_64-linux-gnu
	export GALACTICUS_FCFLAGS="-fintrinsic-modules-path $INSTALL_PATH/finclude -fintrinsic-modules-path $INSTALL_PATH/include -fintrinsic-modules-path $INSTALL_PATH/include/gfortran -fintrinsic-modules-path $INSTALL_PATH/lib/gfortran/modules -L$INSTALL_PATH/lib -L$INSTALL_PATH/lib64 -fuse-ld=bfd"
	export GALACTICUS_CFLAGS="-fuse-ld=bfd"
	export GALACTICUS_CPPFLAGS="-fuse-ld=bfd"
	export GALACTICUS_EXEC_PATH=/opt/galacticus
	export GALACTICUS_DATA_PATH=/opt/datasets

