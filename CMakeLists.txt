## Main CMakeLists.txt of the FAUST project
## Other important files for CMake are
#  - CMakeLists.txt files in the sub-directories
#  - CMake/* files, which define some tools to search for libraries, etc.

cmake_minimum_required(VERSION 2.8.8)
project(FAUST CXX)
set(FAUST_VERSION_MAJOR 1)
set(FAUST_VERSION_MINOR 0)

#Add the c++11 flag, whatever it is
#include(CheckCXXCompilerFlag)
#check_cxx_compiler_flag(-std=c++11 COMPILER_SUPPORTS_CXX11)
#check_cxx_compiler_flag(-std=c++0x COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" CACHE STRING "compile flags" FORCE)
#elseif(COMPILER_SUPPORTS_CXX0X)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" CACHE STRING "compile flags" FORCE)
#else()
#	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has neither C++11 nor c++0x support.")
#endif()



IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/output" CACHE STRING "default install path" FORCE )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


# tmp directory where temporary objects will be located
set(FAUST_TMP_BUILD_DIR "${PROJECT_BINARY_DIR}/tmp")
file(MAKE_DIRECTORY ${FAUST_TMP_BUILD_DIR})
set(FAUST_TMPMEX_DIR "${FAUST_TMP_BUILD_DIR}/mex_obj")
file(MAKE_DIRECTORY ${FAUST_TMPMEX_DIR})

if(UNIX)
	set(MEXOBJ_EXT "o")
elseif(WIN32)
	set(MEXOBJ_EXT ".obj")
endif()


#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin" CACHE STRING "" FORCE)

# SOURCE DIRECTORIES
set(FAUST_EXCEPTION_SRC_DIR ${PROJECT_SOURCE_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_SRC_DIR ${PROJECT_SOURCE_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core CACHE INTERNAL "")
set(FAUST_MATIO_SRC_DIR ${PROJECT_SOURCE_DIR}/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_SRC_DIR ${PROJECT_SOURCE_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_SRC_DIR ${PROJECT_SOURCE_DIR}/test/runtime_comparison CACHE INTERNAL "")
set(FAUST_TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test CACHE INTERNAL "")
set(FAUST_TESTSRC_SRC_DIR ${PROJECT_SOURCE_DIR}/test/src CACHE INTERNAL "")
set(FAUST_TESTDATA_SRC_DIR ${PROJECT_SOURCE_DIR}/test/data CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_SRC_DIR ${PROJECT_SOURCE_DIR}/interface_matlab CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/mex_functions CACHE INTERNAL "")
set(FAUST_EIGEN_SRC_DIR ${PROJECT_SOURCE_DIR}/eigen CACHE INTERNAL "")

# BINARY DIRECTORIES
set(FAUST_EXCEPTION_BIN_DIR ${PROJECT_BINARY_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_BIN_DIR ${PROJECT_BINARY_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_BIN_DIR ${PROJECT_BINARY_DIR}/faust_core CACHE INTERNAL "")
#set(FAUST_MATIO_BIN_DIR ${PROJECT_BINARY_DIR}/testing/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_BIN_DIR ${PROJECT_BINARY_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_TEST_BIN_DIR ${PROJECT_BINARY_DIR}/testing/bin CACHE INTERNAL "")
set(FAUST_TESTSRC_BIN_DIR ${PROJECT_BINARY_DIR}/testing/src CACHE INTERNAL "")
set(FAUST_TESTOUTPUT_BIN_DIR ${PROJECT_BINARY_DIR}/testing/output CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_BIN_DIR ${PROJECT_BINARY_DIR}/testing/runtime_comparison CACHE INTERNAL "")
set(FAUST_TESTMATLABSCRIPT_BIN_DIR ${PROJECT_BINARY_DIR}/testing/matlab_scripts CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_BIN_DIR ${PROJECT_BINARY_DIR}/mex CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "")
file(MAKE_DIRECTORY ${FAUST_TESTOUTPUT_BIN_DIR})

# INSTALL DIRECTORIES (RELATIVES PATHS FROM CMAKE_INSTALL_PREFIX)
set(FAUST_INSTALL_MEX ${CMAKE_INSTALL_PREFIX}/mex)
set(FAUST_INSTALL_LIB ${CMAKE_INSTALL_PREFIX}/lib)
set(FAUST_INSTALL_TESTING_BIN ${CMAKE_INSTALL_PREFIX}/testing/bin)
set(FAUST_INSTALL_TESTING_DATA ${CMAKE_INSTALL_PREFIX}/testing/data)
set(FAUST_INSTALL_TESTRUNCOMP ${CMAKE_INSTALL_PREFIX}/testing/runtime_comparison CACHE INTERNAL "")

# TARGETS
set(FAUST_EXCEPTION_TARGET faust_exception CACHE INTERNAL "")
set(FAUST_MATRIX_TARGET faust_matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_TARGET faust_core CACHE INTERNAL "")
set(FAUST_PALM4MSA_TARGET faust_palm4MSA CACHE INTERNAL "")
set(FAUST_MEX_TARGET faust_mex CACHE INTERNAL "")
set(FAUST_MEXTOOLS_TARGET faust_mextools CACHE INTERNAL "")
set(FAUST_MATIO_TARGET faust_matio CACHE INTERNAL "")
set(FAUST_TARGET faust CACHE INTERNAL "")



if(WIN32)
	if(CMAKE_BUILD_TYPE MATCHES "Debug")
		set(FAUST_DEBUG_RELEASE "debug")
		set(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/Debug)
		set(FAUST_MEXINTERFACE_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR}/Debug)
	else()
		set(FAUST_DEBUG_RELEASE "release")
		set(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/Release)
		set(FAUST_MEXINTERFACE_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR}/Release)
	endif()
else()
	set(FAUST_BIN_DIR ${PROJECT_BINARY_DIR})
endif()



# Default values for FAUST_USE_SINGLEPRECISION and FAUST_USE_OPENBLAS
set(FAUST_USE_SINGLEPRECISION ON CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_OPENBLAS ON CACHE BOOL "Using openBLAS for matrix and vector computations")
set(FAUST_USE_MKL OFF CACHE BOOL "Using MKL for sparse matrix and vector computations")



# find external libraries

if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif()

include(CMake/add_library_path.cmake)
add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{OPENBLASDIR}" "$ENV{MKLDIR}/lib/intel64" "$ENV{MKL_COMPILER_DIR}/lib/intel64" "/usr" "/usr/local" "/opt" "/opt/local")
add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{OPENBLASDIR}" "$ENV{MKLDIR}/include" "/usr" "/usr/local" "/opt" "/opt/local")
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")
message(STATUS "LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")

include(CMake/check_external_libraries.cmake)
check_external_libraries(hdf5 HDF5_LIB_FILE 0)
check_external_libraries(matio MATIO_LIB_FILE 0)
if(FAUST_USE_OPENBLAS)
	check_external_libraries(openblas OPENBLAS_LIB_FILE 1)
	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
endif()
if(FAUST_USE_MKL)
	check_external_libraries(mkl_core MKL_LIB_CORE 1)
	check_external_libraries(mkl_intel_thread MKL_LIB_INTEL_THREAD 1)
	check_external_libraries(iomp5 MKL_LIB_IOMP5 1)
	check_external_libraries(mkl_intel_lp64 MKL_LIB_LP64 1)
	check_external_libraries(m MKL_LIB_M 1)
	check_external_includes("mkl_spblas.h" MKL_INC_DIR 0)
endif()



#check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)
check_external_includes("matio.h" MATIO_INC_DIR 0)

#find_path(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE})
#find_path(EIGEN_LIB_DIR ${EIGEN_LIB_FILE})

# adding compiler flags
if(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fexceptions -fno-omit-frame-pointer -pthread -fPIC " CACHE STRING "compile flags" FORCE)
elseif(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /EHs /openmp /MD" CACHE STRING "compile flags" FORCE)
endif()

if(FAUST_USE_SINGLEPRECISION)
	set(CXX_MEX_FLAGS "${CXX_MEX_FLAGS} -DFAUST_SINGLE")
	if(UNIX)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	elseif(WIN32)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	endif()
endif()

# test if executable matlab is in the path
if(UNIX)
	exec_program("which matlab | xargs echo" OUTPUT_VARIABLE MATLAB_DIR_TMP)
elseif(WIN32)
	exec_program("where matlab.exe" OUTPUT_VARIABLE MATLAB_DIR_TMP)
endif()
if( ${MATLAB_DIR_TMP} MATCHES "matlab")
	if(UNIX)
		string(REGEX REPLACE "([a-zA-Z0-9_/:]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	elseif(WIN32)
		string(REGEX REPLACE "([a-zA-Z0-9_\\:]+)\\\\bin\\\\matlab.exe" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	endif()
	set(MATLAB_ROOT ${MATLAB_ROOT} CACHE PATH "Matlab root directory")
	message(STATUS "MATLAB_ROOT has been found : ${MATLAB_ROOT}")

	set(MATLAB_INCLUDE_DIR "${MATLAB_ROOT}/extern/include" CACHE INTERNAL "Matlab include directory")

	set(MATLAB_ARCH_FILE "${FAUST_TMP_BUILD_DIR}/matlab_arch.txt")
	if(UNIX) 
		# METHODE 1 (without using matlab)
		exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i glnx" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		if("${MEX_SUBDIR_LIB}" STREQUAL "glnxa64")
			set(MEX_EXT  "mexa64")
		elseif("${MEX_SUBDIR_LIB}" STREQUAL "glnx86")
			set(MEX_EXT  "mexa32")
		endif()	

		# METHODE 2 (using matlab)
		#exec_program("matlab -wait -nodesktop -nojvm -nodisplay -r \"fid=fopen('${MATLAB_ARCH_FILE}','w');fprintf(fid,'%s\\n%s\\n',computer('arch'),mexext);fclose(fid);exit\" > ${FAUST_TMP_BUILD_DIR}/matlab_output.log")
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 1" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 2 | tail -n 1" OUTPUT_VARIABLE MEX_EXT)

	elseif(WIN32)
		# METHODE 1 (without using matlab)
		execute_process(COMMAND ${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat "2" "${MATLAB_ROOT}" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		string(REGEX REPLACE "\n" "" MEX_SUBDIR_LIB ${MEX_SUBDIR_LIB})
		if("${MEX_SUBDIR_LIB}" STREQUAL "win64")
			set(MEX_EXT  "mexw64")
		elseif("${MEX_SUBDIR_LIB}" STREQUAL "win32")
			set(MEX_EXT  "mexw32")
		endif()	

		# METHODE 2 (using matlab)
		#exec_program("${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat 0 \"${MATLAB_ARCH_FILE}\"" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		#exec_program("${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat 1 \"${MATLAB_ARCH_FILE}\"" OUTPUT_VARIABLE MEX_EXT)
	else()
		message(WARNING "Unknown type of plateform for matlab")
	endif()
else()
	set(MATLAB_ROOT "" CACHE PATH "Matlab root directory")
	message(WARNING "Matlab executable seems not to be in the path. So \"MATLAB_ROOT\" and \"MATLAB_INCLUDE_DIR\" wont'be defined and mex files won't be compiled. if matlab is installed on your computer, please add matlabroot/bin tu the PATH and try again.")	
endif()

#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")



# Methode 1
add_subdirectory(${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EXCEPTION_BIN_DIR})
add_subdirectory(${FAUST_MATRIX_SRC_DIR} ${FAUST_MATRIX_BIN_DIR})
add_subdirectory(${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_FAUSTCORE_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSA_SRC_DIR} ${FAUST_PALM4MSA_BIN_DIR})
#add_subdirectory(${FAUST_MATIO_SRC_DIR} ${FAUST_MATIO_BIN_DIR})
add_subdirectory(${FAUST_TEST_SRC_DIR} ${FAUST_TEST_BIN_DIR})
add_subdirectory(${FAUST_MEXINTERFACE_SRC_DIR} ${FAUST_MEXINTERFACE_BIN_DIR})
add_library(${FAUST_TARGET} SHARED $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}>)
#target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/glnxa64/libmex.so )


if(FAUST_USE_OPENBLAS)
	if(FAUST_USE_MKL)
		target_link_libraries(${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
	else()
		target_link_libraries(${FAUST_TARGET} ${OPENBLAS_LIB_FILE} ${MKL_LIB_CORE} ${MKL_LIB_INTEL_THREAD} ${MKL_LIB_IOMP5} ${MKL_LIB_LP64} ${MKL_LIB_M})
	endif()	
else()
	if(FAUST_USE_MKL)
		target_link_libraries(${FAUST_TARGET} ${MKL_LIB_CORE} ${MKL_LIB_INTEL_THREAD} ${MKL_LIB_IOMP5} ${MKL_LIB_LP64} ${MKL_LIB_M})
	else()
		target_link_libraries(${FAUST_TARGET} )
	endif()	
endif()

install(TARGETS ${FAUST_TARGET} DESTINATION ${FAUST_INSTALL_LIB})

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









