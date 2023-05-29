#!/bin/bash

# needed env. variables: FAUST_VERSION, DURL, DFILE, EXPERIMENTAL_PKG, NUX_PY_VER

for V in FAUST_VERSION DURL DFILE EXPERIMENTAL_PKG NUX_PY_VER
do
	[[ -z $(env | grep ^$V=) ]] && echo "ERROR: $V variable must be set in the environment." >&2 && exit 1
done

export PYTHON_PATH=$(which python$NUX_PY_VER)
if [[ ! -d 'build' ]]; then mkdir build;fi; cd build

# build python and matlab wrappers separately to use clang for python and gcc for matlab
cmake -DBUILD_WRAPPER_PYTHON=OFF -DBUILD_WRAPPER_MATLAB=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust -DCPACK_PACKAGE_FILE_NAME=faust-$FAUST_VERSION -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG -DNOPY2=ON -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=$PWD/../gpu_mod -DBUILD_FLOAT_PYX=ON -DBUILD_FLOAT_MEX=ON ..
# concise output for make (gitlab output is limited)
make 2>&1 | tee /tmp/log_$(basename $0)_make_$(date +%s) | grep "error:"
cmake -DCMAKE_CXX_COMPILER=clang++ .. # it needs to be made separately because it cleans up other variables (for building both faust.a and python wrapper with clang -- necessary to avoid unresolved c++ symbols which happens by mixing up gcc and clang objects)
cmake -DBUILD_WRAPPER_PYTHON=ON -DBUILD_WRAPPER_MATLAB=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust -DCPACK_PACKAGE_FILE_NAME=faust-$FAUST_VERSION -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG -DNOPY2=ON -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=$PWD/../gpu_mod -DBUILD_FLOAT_PYX=ON  -DBUILD_FLOAT_MEX=ON ..
make clean
# clean might be not enough, explicitly delete object files
find ./ -name "*.o" -delete
make faust 2>&1 | tee /tmp/log_$(basename $0)_make_faust_$(date +%s) | grep -i 'error:\|Building\|Link\|creating'
make faust_python 2>&1 | tee /tmp/log_$(basename $0)_make_faust_python_$(date +%s) | grep -i 'error:\|Building\|Link\|creating'
find ./ -name "*.o" -delete # (more space for rpm/deb generation)
# generate package via cpack
cpack -G RPM -C CPackConfig.cmake
# remove package temporary files
rm -Rf _CPack_Packages/x86_64/RPM
[[ ! -r $(ls "faust-$FAUST_VERSION"*rpm) ]] && echo "Failed to generate RPM package" && exit 2
cpack -G DEB -C CPackConfig.cmake
rm -Rf _CPack_Packages/x86_64/DEB
[[ ! -r $(ls "faust-$FAUST_VERSION"*deb) ]] && echo "Failed to generate DEB package" && exit 3
[[ -n "$COPY_PKG_TO_HOME" ]] && cp faust-$FAUST_VERSION*rpm faust-$FAUST_VERSION*deb $HOME || exit 0
