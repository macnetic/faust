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
# svn checkout --username XXX https://scm.gforge.inria.fr/authscm/XXX/svn/faust/trunk/devcpp/test/CTest_nightly ./
#
###################################################################

# Directory of the library used in the FAUST PROJECT 
export PATH=$PATH:/usr/local/MATLAB/R2014b/bin;
#export PATH=$PATH:usr/local/include

# Directory of the local path of the nightly project
export PATH_DIR_RUN_NIGHTLY='/home/ci/testAL'

#export EIGENDIR='/home/ci/Library/eigen-eigen-07105f7124f9'
#export OPENBLASDIR='/home/ci/local/OPENBLAS'

#export HDF5_ROOT_DIR='/home/ci/local/HDF5/lib'


#export MATIODIR=/usr/local
#export CUDADIR=/usr/local/cuda-6.5

# export version of gcc
# export CC=/usr/lib64/ccache/gcc
# export CXX=/usr/lib64/ccache/g++

# /usr/local/bin/matlab in the PATH 
#export PATH=/usr/local/cuda-7.5/bin:/usr/local/bin:$PATH
#export PATH=/usr/local/cuda-7.5/bin:/usr/local/bin:/usr/lib64/ccache/:$PATH

if [ ! -d "$PATH_DIR_RUN_NIGHTLY" ];
then
	echo "ERROR : $PATH_DIR_RUN_NIGHTLY directory is not defined or do no exist !!!!!! Please select a valid PATH_DIR_RUN_NIGHTLY "
	exit 1
fi



# CTEST OPTION 
export CTEST_CONFIGURE_OPTIONS="
-DFAUST_USE_OPENBLAS=ON;\
-DFAUST_USE_MATIO=ON;\
-DFAUST_USE_XML=ON;\
-DFAUST_USE_MEX=ON;\
-DFAUST_USE_GPU:BOOL=OFF;\
-DFAUST_USE_PROFILING=OFF;\
-DFAUST_ISVERBOSE=OFF;\
-DFAUST_GEN_DOC:BOOL=OFF"

#-DDASH_TESTING:BOOL=ON;\
#-DBUILD_DEBUG:BOOL=OFF;\
#-DBUILD_COVERAGE:BOOL=ON;\

cd $PATH_DIR_RUN_NIGHTLY

rm -rf ./faust_test
rm -rf ./faust_test_build
rm -rf ./faust_output

mkdir ./faust_test
mkdir ./faust_test_build
mkdir ./faust_output

# launch ctest file 
ctest -VV -C Release -S faustTest.cmake 


