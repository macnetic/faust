REM ##################################################################
REM ####################### run_nightly ##############################
REM # This shell-bash SCRIPT performs : 
REM # 	- svn checkout
REM #	- compile the selected target of faust project 
REM #	- run some examples relating to these targets. 
REM # 
REM # faust_test directory contains the sources of the faust project. 
REM # faust_test_build directory contains the build of faust project.
REM # faust_output directory contains the build of faust project.
REM #
REM # WARNING : Please adjust your environment variables, corresponding 
REM # to your computer configurations (library path) 
REM #
REM # svn checkout --username XXX https://scm.gforge.inria.fr/authscm/XXX/svn/faust/trunk/devcpp/test/CTest_nightly
REM #
REM #####set MATLAB_ROOT_DIR="C:\Program Files\MATLAB\R2015b"

REM # cmake -G "MinGw Makefiles" ..
REM # cmake -G "MinGW Makefiles" -Wno-dev .. 


REM # Directory of the local path of the nightly project
set PATH_DIR_RUN_NIGHTLY=C:\Users\ci\FAUST\CTest_nightly

REM # Directory of the library used in the FAUST PROJECT 
REM #export EIGENDIR='/usr/include/eigen3'
REM #export OPENBLASDIR='/opt/OpenBLAS'
REM #export MATIODIR='/usr/local'
REM #export CUDADIR='/usr/local/cuda-7.5'

REM # export version of gcc
set CC=C:\mingw-w64\mingw64\bin\gcc.exe
set CXX=C:\mingw-w64\mingw64\bin\g++.exe

REM The compiler used is MinGW for windows
set (CTEST_CMAKE_GENERATOR_TMP "MinGW Makefiles")
REM set (CTEST_CMAKE_GENERATOR_TMP "Visual Studio 12 2013")



REM set matlab="C:\Program Files\MATLAB\R2015b\bin\matlab.exe"
REM # /usr/local/bin/matlab in the PATH 
REM #export PATH=/usr/local/bin:$PATH
REM # cuda in the PATH
REM #export PATH=/usr/local/cuda-7.5/bin:/usr/lib64/ccache/:$PATH

if not exist %PATH_DIR_RUN_NIGHTLY% Exit 
REM # ( echo "ERROR : %PATH_DIR_RUN_NIGHTLY% directory is not defined or do no exist. Please select a valid PATH_DIR_RUN_NIGHTLY " && Exit )


REM #export HDF5_ROOT_DIR=/home/aleman/Library/hdf5-1.8.16/src/.libs

REM # CTEST OPTION 
set CTEST_CONFIGURE_OPTIONS=-DBUILD_OPENBLAS=OFF; -DBUILD_READ_MAT_FILE=OFF; -DBUILD_READ_XML_FILE=OFF; -DBUILD_MATLAB_MEX_FILES=ON; -DBUILD_USE_GPU:BOOL=OFF; -DFAUST_USE_PROFILING=OFF; -DBUILD_VERBOSE=OFF; -DBUILD_DOCUMENTATION:BOOL=OFF;

REM #-DDASH_TESTING:BOOL=ON;\
REM #-DBUILD_DEBUG:BOOL=OFF;\
REM #-DBUILD_COVERAGE:BOOL=ON;\				

cd %PATH_DIR_RUN_NIGHTLY%
rmdir /s /q faust_test
rmdir /s /q faust_test_build
rmdir /s /q faust_output

mkdir faust_test
mkdir faust_test_build
mkdir faust_output


REM # look for the path of matlab
set "find_exe=matlab.exe"
(where matlab.exe) > logPath.txt 
(where /R "C:\\Program Files\\MATLAB" matlab.exe) >> logPath.txt 
(where /R "C:\\Program Files (x86)\\MATLAB" matlab.exe) >> logPath.txt
set MATLAB_EXE_DIR_TMP=
set /p MATLAB_EXE_DIR_TMP=<logPath.txt
REM for /f "delims=" %%i in ('type logPath.txt') do (set MATLAB_DIR_TMP=%%i && echo %%i)
echo environment variable is defined for matlab Path : %MATLAB_EXE_DIR_TMP%


REM # launch ctest file 
ctest -VV -C Release -S faustTest.cmake 
REM #ctest -VV -C Release -S faustTestContinuous.cmake
