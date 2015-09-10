## Main CMakeLists.txt of the FAUST project
## Other important files for CMake are
#  - CMakeLists.txt files in the sub-directories
#  - CMake/* files, which define some tools to search for libraries, etc.

cmake_minimum_required(VERSION 2.6)
project(FAUST CXX)
set(FAUST_VERSION_MAJOR 1)
set(FAUST_VERSION_MINOR 0)

#Add the c++11 flag, whatever it is
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++11 COMPILER_SUPPORTS_CXX11)
check_cxx_compiler_flag(-std=c++0x COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has neither C++11 nor c++0x support.")
endif()

# SOURCE DIRECTORIES
set(FAUST_EXCEPTION_SRC_DIR ${PROJECT_SOURCE_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_SRC_DIR ${PROJECT_SOURCE_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core CACHE INTERNAL "")
set(FAUST_MATIO_SRC_DIR ${PROJECT_SOURCE_DIR}/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_SRC_DIR ${PROJECT_SOURCE_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_PALM4MSATEST_SRC_DIR ${PROJECT_SOURCE_DIR}/palm4MSA/test CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_SRC_DIR ${PROJECT_SOURCE_DIR}/interface_matlab CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/mex_functions CACHE INTERNAL "")

# BINARY DIRECTORIES
set(FAUST_EXCEPTION_BIN_DIR ${PROJECT_BINARY_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_BIN_DIR ${PROJECT_BINARY_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_BIN_DIR ${PROJECT_BINARY_DIR}/faust_core CACHE INTERNAL "")
set(FAUST_MATIO_BIN_DIR ${PROJECT_BINARY_DIR}/testing/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_BIN_DIR ${PROJECT_BINARY_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_PALM4MSATEST_BIN_DIR ${PROJECT_BINARY_DIR}/testing CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_BIN_DIR ${PROJECT_BINARY_DIR}/interface_matlab CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "")

# TARGETS
set(FAUST_EXCEPTION_TARGET faust_exception CACHE INTERNAL "")
set(FAUST_MATRIX_TARGET faust_matrix CACHE INTERNAL "")
set(FAUST_PALM4MSA_TARGET faust_palm4MSA CACHE INTERNAL "")
set(FAUST_FAUSTCORE_TARGET faust_core CACHE INTERNAL "")
set(FAUST_MEX_TARGET faust_mex CACHE INTERNAL "")
set(FAUST_MEXTOOLS_TARGET faust_mextools CACHE INTERNAL "")

# Default values for FAUST_USE_SINGLEPRECISION and FAUST_USE_OPENBLAS
set(FAUST_USE_SINGLEPRECISION ON CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_OPENBLAS ON CACHE BOOL "Using openBLAS for matrix and vector computations")

# adding compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -fexceptions -fPIC -fno-omit-frame-pointer -pthread")
if(FAUST_USE_SINGLEPRECISION)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE")
	set(FAUST_COMPILE_SINGLE_FLAG "-DFAUST_SINGLE")
endif()
if(FAUST_USE_OPENBLAS)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D__GEMM_WITH_OPENBLAS__")
endif()



# find external libraries
include(CMake/add_library_path.cmake)
add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "opt" "opt/local")
add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{EIGENDIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "opt" "opt/local")
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")

include(CMake/check_external_libraries.cmake)
check_external_libraries(hdf5 HDF5_LIB_FILE 0)
check_external_libraries(matio MATIO_LIB_FILE 0)
if(FAUST_USE_OPENBLAS)
	check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
endif()
check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)
check_external_includes("matio.h" MATIO_INC_DIR 0)

#find_path(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE})
#find_path(EIGEN_LIB_DIR ${EIGEN_LIB_FILE})




#MEX_EXT
if(UNIX)
	set(MEX_EXT "mexa64" CACHE INTERNAL "mex extension")
else()
	set(MEX_EXT "mexw64" CACHE INTERNAL "mex extension")
endif()


if(UNIX)

	exec_program("matlab_dir_tmp=`which matlab`; if echo \"$matlab_dir_tmp\" | grep -q \"/bin/\"; then echo '1'; else echo '0'; fi" OUTPUT_VARIABLE IS_MATLAB_DIR_FOUND)
	if(IS_MATLAB_DIR_FOUND EQUAL 1)
		exec_program("MATLAB_ROOT_TMP=`which matlab |xargs dirname`; while echo \"$MATLAB_ROOT_TMP\" | grep -q \"/bin\"; do MATLAB_ROOT_TMP=`dirname $MATLAB_ROOT_TMP`; done; echo $MATLAB_ROOT_TMP;" OUTPUT_VARIABLE MATLAB_ROOT)
		set(MATLAB_ROOT ${MATLAB_ROOT} CACHE INTERNAL "Matlab root directory")
	else()
		set(MATLAB_ROOT "" CACHE PATH "Matlab root directory")
	endif()

	exec_program("openblas_dir_tmp=`locate libopenblas.so | head -n 1`; if echo \"$openblas_dir_tmp\" | grep -q \"libopenblas\"; then echo '1'; else echo '0'; fi" OUTPUT_VARIABLE IS_OPENBLAS_DIR_FOUND)
	if(IS_OPENBLAS_DIR_FOUND EQUAL 1)
		exec_program("OPENBLAS_LIB_TMP=`locate libopenblas.so | head -n 1 |xargs dirname`; echo $OPENBLAS_LIB_TMP;" OUTPUT_VARIABLE OPENBLAS_LIB_PATH)
		set(OPENBLAS_LIB_PATH ${OPENBLAS_LIB_PATH} CACHE INTERNAL "Openblas library directory")
	else()
		set(OPENBLAS_LIB_PATH "" CACHE PATH "Openblas library directory")
	endif()
else()
	set(MATLAB_ROOT "" CACHE PATH "Matlab root directory")
	set(OPENBLAS_LIB_PATH "" CACHE PATH "Openblas library directory")
endif()

set(MATLAB_INCLUDE_DIR "${MATLAB_ROOT}/extern/include" CACHE INTERNAL "Matlab include directory")


#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")


add_subdirectory(${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EXCEPTION_BIN_DIR})
add_subdirectory(${FAUST_MATRIX_SRC_DIR} ${FAUST_MATRIX_BIN_DIR})
add_subdirectory(${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_FAUSTCORE_BIN_DIR})
add_subdirectory(${FAUST_MATIO_SRC_DIR} ${FAUST_MATIO_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSA_SRC_DIR} ${FAUST_PALM4MSA_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSATEST_SRC_DIR} ${FAUST_PALM4MSATEST_BIN_DIR})
add_subdirectory(${FAUST_MEXINTERFACE_SRC_DIR} ${FAUST_MEXINTERFACE_BIN_DIR})


#add_dependencies(${FAUST_FAUSTCORE_TARGET} ${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})
#add_dependencies(${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})




#set(FAUST_INCLUDE_DIRS "")
#file(GLOB FAUST_CPP "${PROJECT_SOURCE_DIR}/*.cpp")
#file(GLOB FAUST_H "${PROJECT_SOURCE_DIR}/*.h")
#foreach(header ${FAUST_H})
#	get_filename_component(tmp_dir ${header} DIRECTORY)
#	list(APPEND FAUST_INCLUDE_DIRS ${tmp_dir})
#endforeach()
#list(REMOVE_DUPLICATES FAUST_INCLUDE_DIRS)
#include_directories(${FAUST_INCLUDE_DIRS})



#include_directories (BEFORE ${FAUST_MATRIX_DIR} ${FAUST_FAUSTCORE_DIR} 
#	${FAUST_MATIO_DIR} ${FAUST_PALM4MSA_DIR} ${FAUST_MEXINTERFACE_DIR}
#	${FAUST_MEXINTERFACETOOLS_DIR})








