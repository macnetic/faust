##################################################################
####################### run_nightly ##############################
# This shell-bash SCRIPT performs : 
# 	- svn checkout
#	- compile the selected target of faust project 
#	- run some examples relating to these targets. 
# 
# faust_test directory contains the sources of the faust project. 
# faust_test_build directory contains the build of faust project.
# faust_output directory contains the build of faust project.
#
# WARNING : Please adjust your environment variables, corresponding 
# to your computer configurations (library path) 
###################################################################


# Directory of the local path of the nightly project
export PATH_DIR_RUN_NIGHTLY=/home/aleman/WORK/FAUST/faust_nightly


# Directory of the library used in the FAUST PROJECT 
export OPENBLASDIR=/opt/OpenBLAS
export MATIODIR=/usr/local
export CUDADIR=/usr/local/cuda-7.5
export EIGENDIR=/usr/local/include/eigen3

#export HDF5_ROOT_DIR=/home/aleman/Library/hdf5-1.8.16/src/.libs


cd $PATH_DIR_RUN_NIGHTLY
rm -rf faust_test
rm -rf faust_test_build
rm -rf faust_output

mkdir faust_test
mkdir faust_test_build
mkdir faust_output

# launch ctest file 
ctest -VV -C Release -S faustTest.cmake 


