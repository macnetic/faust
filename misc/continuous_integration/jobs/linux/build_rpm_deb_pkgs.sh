#!/bin/bash

# needed env. variables: FAUST_VERSION, DURL, DFILE, EXPERIMENTAL_PKG, NYX_PY_VER

export PYTHON_PATH=$(which python$NUX_PY_VER)
if [[ ! -d 'build' ]]; then  mkdir build;fi; cd build
        # build python and matlab wrappers separately to use clang for python and gcc for matlab
cmake -DBUILD_WRAPPER_PYTHON=OFF -DBUILD_WRAPPER_MATLAB=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust-$FAUST_VERSION -DCPACK_PACKAGE_FILE_NAME=faust-$FAUST_VERSION -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG -DNOPY2=ON -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=$PWD/../gpu_mod -DBUILD_FLOAT_PYX=ON ..
make
cmake -DCMAKE_CXX_COMPILER=clang++ .. # it needs to be made separately because it cleans up other variables (for building both faust.a and python wrapper with clang -- necessary to avoid unresolved c++ symbols which happens by mixing up gcc and clang objects)
cmake -DBUILD_WRAPPER_PYTHON=ON -DBUILD_WRAPPER_MATLAB=ON -DBUILD_DOCUMENTATION=ON -DCMAKE_INSTALL_PREFIX=/opt/local/faust-$FAUST_VERSION -DCPACK_PACKAGE_FILE_NAME=faust-$FAUST_VERSION -DCPACK_PACKAGE_VERSION=$FAUST_VERSION -DBUILD_TESTING=OFF -DREMOTE_DATA_URL="$DURL" -DREMOTE_DATA_FILE="$DFILE" -DEXPERIMENTAL_PKG=$EXPERIMENTAL_PKG -DNOPY2=ON -DUSE_GPU_MOD=ON -DCMAKE_PREFIX_PATH=$PWD/../gpu_mod -DBUILD_FLOAT_PYX=ON ..
make clean
make faust
make faust_python
# generate package via cpack
cpack -G RPM -C CPackConfig.cmake
cpack -G DEB -C CPackConfig.cmake
cp faust-$FAUST_VERSION*rpm faust-$FAUST_VERSION*deb $HOME

