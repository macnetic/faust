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
#
# svn checkout --username XXX https://scm.gforge.inria.fr/authscm/XXX/svn/faust/trunk/devcpp/test/CTest_nightly
#
###################################################################


# Directory of the local path of the nightly project
export PATH_DIR_RUN_NIGHTLY=/home/aleman/Documents/WORK/FAUST/faust_nightly

if [ ! -d "$PATH_DIR_RUN_NIGHTLY" ];
then
	echo "ERROR : $PATH_DIR_RUN_NIGHTLY directory is not defined or do no exist !!!!!! Please select a valid PATH_DIR_RUN_NIGHTLY "
	exit 1
fi

# Directory of the library used in the FAUST PROJECT 
export OPENBLASDIR=/home/aleman/library/OpenBLAS
export MATIODIR=/home/aleman/library/matio-1.5.3/src/.libs
export MATIODIRINC=/home/aleman/library/matio-1.5.3/src
export CUDADIR=/usr/local/cuda-6.5
export EIGENDIR=/usr/local/include/eigen3
# export version of gcc
export CC=/usr/lib64/ccache/gcc
export CXX=/usr/lib64/ccache/g++

# /usr/local/bin/matlab in the PATH 
export PATH=/usr/local/cuda-7.5/bin:/usr/local/bin:$PATH




#export HDF5_ROOT_DIR=/home/aleman/Library/hdf5-1.8.16/src/.libs

# CTEST OPTION 
export CTEST_CONFIGURE_OPTIONS=



cd $PATH_DIR_RUN_NIGHTLY
rm -rf faust_test
rm -rf faust_test_build
rm -rf faust_output
mkdir faust_test
mkdir faust_test_build
mkdir faust_output

# launch ctest file 
ctest -VV -C Release -S faustTest.cmake 


