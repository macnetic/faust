###### Main CMakeLists.txt of the FAUST project ######
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
#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin" CACHE STRING "" FORCE)
##################################################################

IF (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/output" CACHE STRING "default install path" FORCE )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

###### tmp directory where temporary objects will be located ######
set(FAUST_TMP_BUILD_DIR "${PROJECT_BINARY_DIR}/tmp")
file(MAKE_DIRECTORY ${FAUST_TMP_BUILD_DIR})
set(FAUST_TMPMEX_DIR "${FAUST_TMP_BUILD_DIR}/mex_obj")
file(MAKE_DIRECTORY ${FAUST_TMPMEX_DIR})

if(UNIX)
	set(MEXOBJ_EXT "o")
elseif(WIN32)
	set(MEXOBJ_EXT ".obj")
endif()

###### Different option of configuration ######
set(FAUST_USE_SINGLEPRECISION OFF CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_OPENBLAS ON CACHE BOOL "Using openBLAS for matrix and vector computations")
set(FAUST_USE_MEX ON CACHE BOOL "Generate Mexfiles")
set(FAUST_USE_PROFILING OFF CACHE BOOL "Profiling the code")
set(FAUST_USE_MATIO ON CACHE BOOL "Using matio library to read matfile")
set(FAUST_USE_XML OFF CACHE BOOL "Using xml configuration file")
set(FAUST_ISVERBOSE OFF CACHE BOOL "Display useful message for debugging")
set(FAUST_USE_GPU OFF CACHE BOOL "Using both CPU and GPU process")
set(FAUST_GEN_DOC OFF CACHE BOOL "Generate html documentation with doxygen")

if (FAUST_USE_SINGLEPRECISION)
	message(STATUS "**********SINGLEPRECISION ACTIF************")
endif()
if (FAUST_USE_OPENBLAS)
	message(STATUS "**********OPENBLAS ACTIF************")
endif()
if (FAUST_USE_MEX)
	message(STATUS "**********MEX ACTIF************")
endif()
if (FAUST_USE_PROFILING)
	message(STATUS "**********PROFILING ACTIF************")
endif()
if (FAUST_USE_MATIO)
	message(STATUS "**********MATIO ACTIF************")
endif()
if (FAUST_USE_XML)
	message(STATUS "**********XML ACTIF************")
endif()
if (FAUST_ISVERBOSE)
	message(STATUS "**********VERBOSE ACTIF************")
endif()
if (FAUST_USE_GPU)
	message(STATUS "**********GPU ACTIF************")
endif()
if (FAUST_GEN_DOC)
	message(STATUS "**********DOC ACTIF************")
endif()
##################################################################

###### Enable testing: ######
OPTION(BUILD_TESTING "Build Tests" ON) 
include(CTest)
enable_testing()
##################################################################

###### original source DIRECTORIES ######
set(FAUST_EXCEPTION_SRC_DIR ${PROJECT_SOURCE_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_SRC_DIR ${PROJECT_SOURCE_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core CACHE INTERNAL "")
set(FAUST_MATIO_SRC_DIR ${PROJECT_SOURCE_DIR}/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_SRC_DIR ${PROJECT_SOURCE_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_SRC_DIR ${FAUST_TEST_SRC_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_TEST_SRC_TOOL_DIR ${FAUST_TEST_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_TESTSRC_SRC_DIR ${FAUST_TEST_SRC_DIR}/src CACHE INTERNAL "")
set(FAUST_TESTDATA_SRC_DIR ${FAUST_TEST_SRC_DIR}/data CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_SRC_DIR ${PROJECT_SOURCE_DIR}/interface_matlab CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/mex_functions CACHE INTERNAL "")
#set(FAUST_EIGEN_SRC_DIR ${PROJECT_SOURCE_DIR}/eigen CACHE INTERNAL "")
if (FAUST_USE_GPU)
	set(FAUST_MATRIX_CU_SRC_DIR ${PROJECT_SOURCE_DIR}/matrix/matrix_cu CACHE INTERNAL "")
	set(FAUST_FAUSTCORE_CU_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core/faust_core_cu CACHE INTERNAL "")
	set(FAUST_PALM4MSA_CU_SRC_DIR ${PROJECT_SOURCE_DIR}/palm4MSA/palm4MSA_cu CACHE INTERNAL "")
	set(FAUST_TEST_GPU_SRC_DIR ${PROJECT_SOURCE_DIR}/test/GPU_test CACHE INTERNAL "")
endif(FAUST_USE_GPU)
##################################################################


###### BINARY DIRECTORIES ######
# MEX function
set(FAUST_WRAPPER_BIN_DIR ${PROJECT_BINARY_DIR}/wrapper)
set(FAUST_MEXINTERFACE_BIN_DIR ${FAUST_WRAPPER_BIN_DIR}/mex CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "") #  a MODIFIER ET METTRE LES MEXFONCTIONS.mexa DEDANS

## examples directories including some examples of the use of the library FAuST
set(FAUST_TEST_DIR ${FAUST_WRAPPER_BIN_DIR}/run_test CACHE INTERNAL "")
set(FAUST_TEST_BIN_DIR ${FAUST_TEST_DIR}/bin CACHE INTERNAL "")
set(FAUST_TESTSRC_BIN_DIR ${FAUST_TEST_DIR}/src CACHE INTERNAL "")
set(FAUST_TESTTMP_BIN_DIR ${FAUST_TEST_DIR}/tmp CACHE INTERNAL "")
set(FAUST_TESTOUTPUT_BIN_DIR ${FAUST_TEST_DIR}/output CACHE INTERNAL "")
set(FAUST_TESTFIG_BIN_DIR ${FAUST_TEST_DIR}/fig CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_BIN_DIR ${FAUST_TEST_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_TESTTOOLS_BIN_DIR ${FAUST_TEST_DIR}/tools CACHE INTERNAL "")


# sources directorie
set(FAUST_EXCEPTION_BIN_DIR ${PROJECT_BINARY_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_BIN_DIR ${PROJECT_BINARY_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_BIN_DIR ${PROJECT_BINARY_DIR}/faust_core CACHE INTERNAL "")
set(FAUST_PALM4MSA_BIN_DIR ${PROJECT_BINARY_DIR}/palm4MSA CACHE INTERNAL "")
#set(FAUST_MATIO_BIN_DIR ${PROJECT_BINARY_DIR}/testing/init_from_matlab CACHE INTERNAL "")
if (FAUST_USE_GPU)
	set(FAUST_MATRIX_CU_BIN_DIR ${FAUST_MATRIX_BIN_DIR}/matrix_cu CACHE INTERNAL "")
	set(FAUST_FAUSTCORE_CU_BIN_DIR ${FAUST_FAUSTCORE_BIN_DIR}/faust_core_cu CACHE INTERNAL "")
	set(FAUST_PALM4MSA_CU_BIN_DIR ${FAUST_PALM4MSA_BIN_DIR}/palm4MSA_cu CACHE INTERNAL "")
	set(FAUST_TEST_GPU_BIN_DIR ${FAUST_TEST_DIR}/binGPU CACHE INTERNAL "") 
endif(FAUST_USE_GPU)






###### DEMO DIRECTORIES ######
set(FAUST_DEMO_SRC_DIR ${PROJECT_SOURCE_DIR}/demo CACHE INTERNAL "")
set(FAUST_DEMODATA_SRC_DIR ${FAUST_DEMO_SRC_DIR}/data CACHE INTERNAL "")
set(FAUST_DEMOTOOLS_SRC_DIR ${FAUST_DEMO_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_DEMO_BIN_DIR ${FAUST_DEMO_SRC_DIR} CACHE INTERNAL "")
set(FAUST_DEMOTOOLS_BIN_DIR ${PROJECT_BINARY_DIR}/tools CACHE INTERNAL "")
##################################################################

###### CMDLINE DIRECTORIES ######
set(FAUST_CMDLINE_SRC_DIR ${PROJECT_SOURCE_DIR}/cmd_line CACHE INTERNAL "")
set(FAUST_CMDLINE_SRC_DATA_DIR ${FAUST_CMDLINE_SRC_DIR}/data CACHE INTERNAL "")
set(FAUST_CMDLINE_SRC_CONFIGFILE_DIR ${FAUST_CMDLINE_SRC_DIR}/config_file CACHE INTERNAL "")
set(FAUST_CMDLINE_TOOLS_SRC_DIR ${FAUST_CMDLINE_SRC_DIR}/xml_tools CACHE INTERNAL "")
set(FAUST_CMDLINE_FUNC_SRC_DIR ${FAUST_CMDLINE_SRC_DIR}/cmdline_function CACHE INTERNAL "")
set(FAUST_CMDLINE_BIN_SRC_DIR ${FAUST_WRAPPER_BIN_DIR}/cmd_line/src CACHE INTERNAL "")
set(FAUST_CMDLINE_BIN_DIR ${FAUST_WRAPPER_BIN_DIR}/cmd_line/bin CACHE INTERNAL "")
##################################################################

###### DOCUMENTATION DIRECTORIES ######
set(FAUST_DOC_SRC_DIR ${PROJECT_SOURCE_DIR}/gen_doc CACHE INTERNAL "")
set(FAUST_DOC_BIN_DIR ${PROJECT_BINARY_DIR}/doc CACHE INTERNAL "")
##################################################################


file(MAKE_DIRECTORY ${FAUST_TESTOUTPUT_BIN_DIR})
file(MAKE_DIRECTORY ${FAUST_TESTFIG_BIN_DIR})
file(MAKE_DIRECTORY ${FAUST_TESTTMP_BIN_DIR})

###### INSTALL DIRECTORIES (RELATIVES PATHS FROM CMAKE_INSTALL_PREFIX) ######
set(FAUST_INSTALL_MEX ${CMAKE_INSTALL_PREFIX}/mex)
set(FAUST_INSTALL_LIB ${CMAKE_INSTALL_PREFIX}/lib)
set(FAUST_INSTALL_CMDLINE ${CMAKE_INSTALL_PREFIX}/CmdLine)
set(FAUST_INSTALL_TESTING_BIN ${CMAKE_INSTALL_PREFIX}/testing/bin)
set(FAUST_INSTALL_TESTING_DATA ${CMAKE_INSTALL_PREFIX}/testing/data)
set(FAUST_INSTALL_TESTRUNCOMP ${CMAKE_INSTALL_PREFIX}/testing/runtime_comparison CACHE INTERNAL "")
##################################################################

###### TARGETS ######
set(FAUST_EXCEPTION_TARGET faust_exception CACHE INTERNAL "")
set(FAUST_MATRIX_TARGET faust_matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_TARGET faust_core CACHE INTERNAL "")
set(FAUST_PALM4MSA_TARGET faust_palm4MSA CACHE INTERNAL "")
set(FAUST_MEX_TARGET faust_mex CACHE INTERNAL "")
set(FAUST_MEXTOOLS_TARGET faust_mextools CACHE INTERNAL "")
set(FAUST_MATIO_TARGET faust_matio CACHE INTERNAL "")
set(FAUST_XML_TARGET faust_xml CACHE INTERNAL "")
set(FAUST_TARGET faust CACHE INTERNAL "")
if (FAUST_USE_GPU)
	set(FAUST_MATRIX_CU_TARGET faust_matrix_cu CACHE INTERNAL "")
	set(FAUST_FAUSTCORE_CU_TARGET faust_core_cu CACHE INTERNAL "")
	set(FAUST_PALM4MSA_CU_TARGET faust_palm4MSA_cu CACHE INTERNAL "")
endif(FAUST_USE_GPU)
##################################################################


###### DEBUG vs RELEASE TARGETS ######
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
endif(WIN32)

# adding compiler flags 
set(CMAKE_CXX_FLAGS "")
if( (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
	if (UNIX)	
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} -O1")
	elseif(WIN32)
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} /O0")
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")
	endif() 
	message(STATUS "**********DEBUG mode************")
else()
	if (UNIX)	
		set(CMAKE_CXX_FLAGS " -O2 ")
	elseif(WIN32)
		set(CMAKE_CXX_FLAGS " /O2 ")
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
	endif()
	message(STATUS "**********RELEASE mode************")
endif()

if(UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fexceptions -fno-omit-frame-pointer -pthread -fPIC " CACHE STRING "compile flags" FORCE)
elseif(WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /EHs /openmp /MD" CACHE STRING "compile flags" FORCE)
else()
	message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
endif(UNIX)

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
##################################################################



###### find external libraries ######
if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif(WIN32)

# looking for library and include files in the environment variables
include(CMake/add_library_path.cmake)
add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/local/lib" "/opt" "/opt/local" "/usr/lib/x86_64-linux-gnu/" )
add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "$ENV{EIGENDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" )

set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")

message(STATUS "LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")

include(CMake/check_external_libraries.cmake)

check_external_libraries(matio MATIO_LIB_FILE 0)
check_external_libraries(xml2 XML2_LIB_FILE 0)
check_external_libraries(hdf5 HDF5_LIB_FILE 1)

#LDFLAGS = -L$(CUDADIR)/lib64 -L$(MATIODIR)/lib -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio -lhdf5
check_external_includes("libxml/parser.h" XML_INC_DIR 0)
check_external_includes("matio.h" MATIO_INC_DIR 0)

check_external_includes("Eigen/Dense" EIGEN_INC_DIR 1)


if (FAUST_USE_GPU)
	check_external_libraries(cublas CUBLAS_LIB_FILE 1)
	check_external_libraries(cudart CUDART_LIB_FILE 1)
	check_external_libraries(cusparse CUSPARSE_LIB_FILE 1)

	check_external_includes("cublas_v2.h" CUBLAS_V2_INC_DIR 1)
	check_external_includes("cusparse.h" CUSPARSE_INC_DIR 1)
	check_external_includes("cuda.h" CUDA_INC_DIR 1)
	check_external_includes("cuda_runtime_api.h" CUDA_RUNTIME_API_INC_DIR 1)
endif(FAUST_USE_GPU)

if(FAUST_USE_OPENBLAS)
	check_external_libraries(openblas OPENBLAS_LIB_FILE 1)
	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
endif(FAUST_USE_OPENBLAS)
##################################################################

#find_path(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE})
#find_path(EIGEN_LIB_DIR ${EIGEN_LIB_FILE})
# if(FAUST_USE_SINGLEPRECISION)
	# set(CXX_MEX_FLAGS "${CXX_MEX_FLAGS} -DFAUST_SINGLE")
	# if(UNIX)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# elseif(WIN32)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# endif()
	# message(STATUS "**** SINGLE PRECISION USE *****")
# endif()


###### test if executable matlab is in the path ######
if(UNIX)
	#message(STATUS "MATLAB_DIR_TMP 1 = ${MATLAB_DIR_TMP}")
	exec_program("which matlab | xargs echo" OUTPUT_VARIABLE MATLAB_DIR_TMP)
	#message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")
   	exec_program("readlink ${MATLAB_DIR_TMP}" OUTPUT_VARIABLE READLINK_TMP)
	#message(STATUS "READLINK_TMP = ${READLINK_TMP}")
	if(${READLINK_TMP} MATCHES matlab)
		set(MATLAB_DIR_TMP ${READLINK_TMP})
		#message(STATUS "MATLAB_DIR_TMP 3 = ${MATLAB_DIR_TMP}")
   	endif()
	#message(STATUS "MATLAB_DIR_TMP 4 = ${MATLAB_DIR_TMP}")

elseif(WIN32)
	exec_program("where matlab.exe" OUTPUT_VARIABLE MATLAB_DIR_TMP)
else()
	message(WARNING "Unknown type of plateform for matlab")		
endif()
if( ${MATLAB_DIR_TMP} MATCHES "matlab")
	if(UNIX)
		string(REGEX REPLACE "([a-zA-Z0-9_/:]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	elseif(WIN32)
		string(REGEX REPLACE "([a-zA-Z0-9_\\:]+)\\\\bin\\\\matlab.exe" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	else()
		message(WARNING "Unknown type of plateform for matlab")	
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
##################################################################



###### TARGET GENERATION ######
add_subdirectory(${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EXCEPTION_BIN_DIR})
add_subdirectory(${FAUST_MATRIX_SRC_DIR} ${FAUST_MATRIX_BIN_DIR})
add_subdirectory(${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_FAUSTCORE_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSA_SRC_DIR} ${FAUST_PALM4MSA_BIN_DIR})
#add_subdirectory(${FAUST_MATIO_SRC_DIR} ${FAUST_MATIO_BIN_DIR})
add_subdirectory(${FAUST_TEST_SRC_DIR} ${FAUST_TEST_BIN_DIR})
add_subdirectory(${FAUST_DOC_SRC_DIR} ${FAUST_DOC_BIN_DIR})
add_subdirectory(${FAUST_CMDLINE_SRC_DIR} ${FAUST_CMDLINE_BIN_DIR})
if (FAUST_USE_GPU)
	add_subdirectory(${FAUST_MATRIX_CU_SRC_DIR} ${FAUST_MATRIX_CU_BIN_DIR})
	add_subdirectory(${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_FAUSTCORE_CU_BIN_DIR})
	add_subdirectory(${FAUST_PALM4MSA_CU_SRC_DIR} ${FAUST_PALM4MSA_CU_BIN_DIR})
	add_subdirectory(${FAUST_TEST_GPU_SRC_DIR} ${FAUST_TEST_GPU_BIN_DIR})
endif (FAUST_USE_GPU)

if (FAUST_USE_MEX)
	add_subdirectory(${FAUST_MEXINTERFACE_SRC_DIR} ${FAUST_MEXINTERFACE_BIN_DIR})
	add_subdirectory(${FAUST_DEMO_SRC_DIR} ${FAUST_DEMO_BIN_DIR})
endif()
##################################################################
file(GLOB MATRIX_HPP_FILES "matrix/*.hpp")


###### TARGET GENERATION ######
if (FAUST_USE_GPU)
	if (FAUST_USE_MEX)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_CU_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_CU_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (FAUST_USE_MEX)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_CU_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_CU_TARGET}>)
	endif(FAUST_USE_MEX)
else (FAUST_USE_GPU)
	if (FAUST_USE_MEX)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (FAUST_USE_MEX)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}>)
	endif(FAUST_USE_MEX)
endif (FAUST_USE_GPU)

target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/glnxa64/libmex.so )

# In case of GPU, the "CUDA" library libfaust_matrix_cu.a is linked with faust (libfaust_matrix_cu.a is not a OBJECT library) 
if (FAUST_USE_GPU)
	target_link_libraries(${FAUST_TARGET} ${FAUST_MATRIX_CU_TARGET} )	
endif(FAUST_USE_GPU)


if(FAUST_USE_OPENBLAS)
	target_link_libraries(${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
else()
	target_link_libraries(${FAUST_TARGET} )	
endif(FAUST_USE_OPENBLAS)

install(TARGETS ${FAUST_TARGET} DESTINATION ${FAUST_INSTALL_LIB})
##################################################################

#add_dependencies(${FAUST_FAUSTCORE_TARGET} ${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})
#add_dependencies(${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})







