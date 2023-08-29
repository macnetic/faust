#!/usr/bin/env bash

DEFAULT_VERSIONS="11.4 12.1"

if [ $# -lt 1 ]
then
	CUDA_VERSIONS="$DEFAULT_VERSIONS"
else
	CUDA_VERSIONS=$*
fi

for CU_VER in $CUDA_VERSIONS
do
	cd gpu_mod; if [[ ! -d build-$CU_VER ]]; then mkdir build-cu$CU_VER; fi; cd build-cu$CU_VER
	cmake -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-$CU_VER -DCMAKE_CUDA_COMPILER=/usr/local/cuda-$CU_VER/bin/nvcc .. && make -j8
	cd ../..
done
