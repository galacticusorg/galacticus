export GALACTICUS_EXEC_PATH=`pwd`
export FCCOMPILER=gfortran-mp-12
export CCOMPILER=gcc-mp-12
export CPPCOMPILER=g++-mp-12
export GALACTICUS_FCFLAGS="-fintrinsic-modules-path /usr/local/include -fintrinsic-modules-path /usr/local/finclude -L/usr/local/lib -L/opt/local/lib"
if [[ "${ver}" -eq 13 ]]; then
    export GALACTICUS_FCFLAGS="$GALACTICUS_FCFLAGS -Wl,-ld_classic"
fi
export GALACTICUS_CFLAGS="-I/usr/local/include -I/opt/local/include"
export GALACTICUS_CPPFLAGS="-I/usr/local/include -I/opt/local/include"
make -j16 Galacticus.exe

