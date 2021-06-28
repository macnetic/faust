#!/usr/bin/env bash


# This script allows to count lines of code as part of three sets: 1. C++ 2. Matlab 3. Python
# The goal is to keep the track of how the code is counted for the APP deposit

FAUST_PATH=$1

[[ -z "$FAUST_PATH" || ! -d "$FAUST_PATH" ]] && echo "USAGE: $0 <faust_path> [-v]\n-v: verbose mode (all files are listed)." && exit 1

[[ $* = *-v* ]] && VERBOSE=1 || VERBOSE=0

echo "C++ code number of lines:"
wc -l $(find "$FAUST_PATH/src" "$FAUST_PATH/misc/test/src/C++" "$FAUST_PATH/wrapper/matlab/src/" "$FAUST_PATH/wrapper/python/src/" "$FAUST_PATH/gpu_mod/src" | grep -v "matlab.*Cplx.cpp$\|matlab.*Real.cpp" | grep -i "\(.hpp\|.cpp\|.cpp.in\|.hpp.in\)$") | ([[ ! $VERBOSE = 1 ]] && tail -1 || tee /dev/null)
echo "Python code number of lines:"
wc -l $(find "$FAUST_PATH/wrapper/python" "$FAUST_PATH/misc/test/src/Python" -name "*.py" -o -name "*.pyx" -o -name "*.pyx.in" -o -name "*.py") | ( [[ ! $VERBOSE = 1 ]] && tail -1 || tee -a /dev/null)
echo "Matlab code number of lines:"
wc -l $(find "$FAUST_PATH/wrapper/matlab" "$FAUST_PATH/misc/test/src/Matlab" -name "*.m.in" -o -name "*.m") | ([[ ! $VERBOSE = 1 ]] && tail -1 || tee -a /dev/null)




