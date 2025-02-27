#!/bin/bash
#
# needed env. variables: FAUST_VERSION, DURL, DFILE, EXPERIMENTAL_PKG, NUX_PY_VER
# LINUX_HDF5_SLIB_PATH LINUX_ZLIB_SLIB_PATH LINUX_MATIO_SLIB_PATH

for V in FAUST_VERSION DURL DFILE EXPERIMENTAL_PKG NUX_PY_VER LINUX_HDF5_SLIB_PATH LINUX_ZLIB_SLIB_PATH LINUX_MATIO_SLIB_PATH
do
	[[ -z $(env | grep ^$V=) ]] && echo "ERROR: $V variable must be set in the environment." >&2 && exit 1
done

if [[ ! -d 'build' ]]; then mkdir build;fi; cd build
export PYTHON_PATH=$(which python$NUX_PY_VER)
cmake -DBUILD_WRAPPER_PYTHON=ON -DBUILD_WRAPPER_MATLAB=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust -DCPACK_PACKAGE_FILE_NAME=faust-$FAUST_VERSION-static -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DEXCLUDE_FAUST_LIB_INSTALL=ON -DUSE_MATIO_STATIC_LIBS=ON -DMATIO_STATIC_LIB_PATH=$LINUX_MATIO_SLIB_PATH -DZ_STATIC_LIB_PATH=$LINUX_ZLIB_SLIB_PATH -DHDF5_STATIC_LIB_PATH=$LINUX_HDF5_SLIB_PATH -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DBUILD_MULTITHREAD=ON -DNOPY2=ON -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=$PWD/../gpu_mod -DBUILD_FLOAT_PYX=ON -DBUILD_FLOAT_MEX=ON -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG ..
# concise output for make (gitlab output is limited)
make -j8 2>&1 | tee /tmp/log_$(basename $0)_make_$(date +%s) | grep "error:\|Building\|Link\|creating"
make clean
# clean might be not enough, explicitly delete object files (more space for rpm/deb generation)
find ./ -name "*.o" -delete
cpack -G RPM -C CPackConfig.cmake
# remove package temporary files
rm -Rf _CPack_Packages/x86_64/RPM
[[ ! -r $(ls "faust-$FAUST_VERSION"*rpm) ]] && echo "Failed to generate RPM package" && exit 2
cpack -G DEB -C CPackConfig.cmake
rm -Rf _CPack_Packages/x86_64/DEB
[[ ! -r $(ls "faust-$FAUST_VERSION"*deb) ]] && echo "Failed to generate DEB package" && exit 3
[[ -n "$COPY_PKG_TO_HOME" ]] && cp faust-$FAUST_VERSION*rpm faust-$FAUST_VERSION*.deb $HOME || exit 0
