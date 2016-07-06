###### Main CMakeLists.txt of the FAUST project ######
## Other important files for CMake are
#  - CMakeLists.txt files in the sub-directories
#  - CMake/* files, which define some tools to search for libraries, etc.

cmake_minimum_required(VERSION 2.8.8)
project(FAUST CXX)

#set(FAUST_VERSION_MAJOR 1)
#set(FAUST_VERSION_MINOR 0)
set(FAUST_VERSION 2.0)


if(UNIX)
	#message(STATUS "UNIX OPERATING SYSTEM")
	#if(LINUX) DON'T WORK SO AN UNIX OS WHICH IS NOT APPLE IS ASSUMED TO BE LINUX OS 	
	if(APPLE)
		message(STATUS "APPLE OPERATING SYSTEM")
	else(APPLE)
		message(STATUS "LINUX OPERATING SYSTEM")
	endif(APPLE)
else(UNIX)
	message(STATUS "WINDOWS OPERATING SYSTEM")
endif(UNIX)

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


#cmake .. -DCMAKE_INSTALL_PREFIX="/home/aleman/Faust_install" 
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
	SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/build/FAuST_OUTPUT" CACHE STRING "default install path" FORCE )
else()
	SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
	message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


##############################################################################
# Setting compilation option
##############################################################################
# Allows building executable:
#option(BUILD_EXECUTABLE "Build executable." ON)
# -->  matio  

# Allows building of shared libraries:
#option(BUILD_EXTERNALS_LIBS "Build externals libraries if they are not installed." ON)
# Enable testing:
option(BUILD_TESTING "Build Tests" ON) 
#Allows generating the doc:
option(BUILD_DOCUMENTATION "Generating the doxygen documentation." OFF)
# Allows multithreading:
option(BUILD_MULTITHREAD "Enable multithread with OpenMP Multithreading." OFF)
# Allows verbose:
option(BUILD_VERBOSE "Enable verbose option when compile (-v)." OFF)
# Allows MPTK Debug:
option(BUILD_DEBUG "Enable MPTK Debug mode." OFF)
# Allows Matlab MEX files :
option(BUILD_USE_GPU  "Using both CPU and GPU process" OFF)

option(BUILD_MATLAB_MEX_FILES "Enable building Matlab MEX files." ON)
# Using openBLAS for matrix and vector computations
option(BUILD_OPENBLAS "Using openBLAS for matrix and vector computations" OFF)
option(BUILD_READ_XML_FILE "Using xml configuration to read xml files" OFF)
option(BUILD_READ_MAT_FILE "Using matio library to read mat files" OFF)


#CMAKE_DEPENDENT_OPTION

###### Different option of configuration ######
#option(BUILD_TESTING "Build Tests" ON) 
set(FAUST_USE_SINGLEPRECISION OFF CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_PROFILING OFF CACHE BOOL "Profiling the code")
#set(BUILD_MULTITHREAD OFF CACHE BOOL "OpenMP Multithreading")
#set(BUILD_VERBOSE OFF CACHE BOOL "Display useful message for debugging")
#set(BUILD_USE_GPU OFF CACHE BOOL "Using both CPU and GPU process")
#set(BUILD_DOCUMENTATION OFF CACHE BOOL "Generate html documentation with doxygen")

#set(BUILD_READ_XML_FILE OFF CACHE BOOL "Using xml configuration to read xml files")
#set(BUILD_MATLAB_MEX_FILES ON CACHE BOOL "Generate Mexfiles")
#set(BUILD_OPENBLAS OFF CACHE BOOL "Using openBLAS for matrix and vector computations")
#set(BUILD_READ_MAT_FILE OFF CACHE BOOL "Using matio library to read mat files")

if (FAUST_USE_SINGLEPRECISION)
	message(STATUS "**********SINGLEPRECISION ACTIF************")
endif()
if (BUILD_OPENBLAS)
	message(STATUS "**********OPENBLAS ACTIF************")
endif()
if (BUILD_MATLAB_MEX_FILES)
	message(STATUS "**********MEX ACTIF************")
endif()
if (FAUST_USE_PROFILING)
	message(STATUS "**********PROFILING ACTIF************")
endif()
if (BUILD_READ_MAT_FILE)
	message(STATUS "**********MATIO ACTIF************")
endif()
if (BUILD_READ_XML_FILE)
	message(STATUS "**********XML ACTIF************")
endif()
if (BUILD_VERBOSE)
	message(STATUS "**********VERBOSE ACTIF************")
endif()
if (BUILD_USE_GPU)
	message(STATUS "**********GPU ACTIF************")
endif()
if(BUILD_MULTITHREAD)
	message(STATUS "**********OPENMP ACTIF*********")
endif()
if (BUILD_DOCUMENTATION)
	message(STATUS "**********DOC ACTIF************")
endif()
##################################################################




###### Enable testing: ######
include(CTest)
enable_testing()
##################################################################



###### SOURCES DIRECTORIES ######
set(FAUST_SRC_DIR ${PROJECT_SOURCE_DIR}/src CACHE INTERNAL "")
set(FAUST_SRC_LINEAR_OPERATOR_DIR ${FAUST_SRC_DIR}/faust_linear_operator CACHE INTERNAL "")
set(FAUST_LINEAR_OPERATOR_CPU_SRC_DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}/CPU CACHE INTERNAL "")

set(FAUST_BIN_LINEAR_OPERATOR_DIR ${PROJECT_BINARY_DIR}/faust_linear_operator CACHE INTERNAL "")
set(FAUST_LINEAR_OPERATOR_CPU_BIN_DIR ${FAUST_BIN_LINEAR_OPERATOR_DIR}/CPU CACHE INTERNAL "")

set(FAUST_ALGORITHM_SRC_DIR ${FAUST_SRC_DIR}/algorithm CACHE INTERNAL "")
set(FAUST_ALGORITHM_CONSTRAINT_SRC_DIR ${FAUST_ALGORITHM_SRC_DIR}/constraint CACHE INTERNAL "")
set(FAUST_ALGORITHM_FACTORIZATION_SRC_DIR ${FAUST_ALGORITHM_SRC_DIR}/factorization CACHE INTERNAL "")

set(FAUST_ALGORITHM_BIN_DIR ${PROJECT_BINARY_DIR}/algorithm CACHE INTERNAL "")
set(FAUST_ALGORITHM_CONSTRAINT_BIN_DIR ${FAUST_ALGORITHM_BIN_DIR}/constraint CACHE INTERNAL "")
set(FAUST_ALGORITHM_FACTORIZATION_BIN_DIR ${FAUST_ALGORITHM_BIN_DIR}/factorization CACHE INTERNAL "")

#set(FAUST_EXCEPTION_SRC_DIR ${FAUST_SRC_DIR}/error_display CACHE INTERNAL "")
#set(FAUST_FAUSTCORE_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core CACHE INTERNAL "")


##############  INSTALL DIRECTORIES   #############################
set(FAUST_INSTALL_LIB 			${CMAKE_INSTALL_PREFIX}/lib)
set(FAUST_INSTALL_BIN 			${CMAKE_INSTALL_PREFIX}/bin)
set(FAUST_INSTALL_BIN_TIME_COMP ${CMAKE_INSTALL_PREFIX}/bin/runtime_comparison)
set(FAUST_INSTALL_DOC 			${CMAKE_INSTALL_PREFIX}/doc)

set(FAUST_INSTALL_CMDLINE 		${CMAKE_INSTALL_PREFIX}/CmdLine)
set(FAUST_INSTALL_TESTING_DATA 	${CMAKE_INSTALL_PREFIX}/testing/data)
##################################################################


###### WRAPPER DIRECTORIES ######
set(FAUST_SRC_WRAPPER_DIR ${PROJECT_SOURCE_DIR}/wrapper CACHE INTERNAL "")
set(FAUST_BIN_WRAPPER_DIR ${PROJECT_BINARY_DIR}/wrapper CACHE INTERNAL "")
###### WRAPPER CMDLINE SRC DIRECTORIES ######
set(FAUST_SRC_CMDLINE_DIR 			${FAUST_SRC_WRAPPER_DIR}/cmd_line CACHE INTERNAL "")
set(FAUST_CMDLINE_FUNC_SRC_DIR 		${FAUST_SRC_WRAPPER_DIR}/cmd_line/src CACHE INTERNAL "")
set(FAUST_CMDLINE_TYPE_FORMAT_DIR 	${FAUST_SRC_WRAPPER_DIR}/cmd_line/type_format CACHE INTERNAL "")
set(FAUST_CMDLINE_TYPE_FORMAT_XML_SRC_DIR ${FAUST_CMDLINE_TYPE_FORMAT_DIR}/xml CACHE INTERNAL "")
set(FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR ${FAUST_CMDLINE_TYPE_FORMAT_DIR}/mat CACHE INTERNAL "")


###### WRAPPER CMDLINE BUILD DIRECTORIES ######
set(FAUST_BIN_CMDLINE_DIR  		${FAUST_BIN_WRAPPER_DIR}/cmd_line CACHE INTERNAL "")
set(FAUST_BIN_CMDLINE_SRC_DIR 	${FAUST_BIN_WRAPPER_DIR}/cmd_line/src CACHE INTERNAL "")
set(FAUST_BIN_CMDLINE_BIN_DIR 	${FAUST_BIN_WRAPPER_DIR}/cmd_line/bin CACHE INTERNAL "")
##################################################################

## fichier d'entrée doivent être renseigné dans l'executable
###### WRAPPER CMDLINE DATA DIRECTORIES ######
set(FAUST_MISC_DIR ${PROJECT_SOURCE_DIR}/misc CACHE INTERNAL "")
set(FAUST_DATA_DIR 			${FAUST_MISC_DIR}/data CACHE INTERNAL "")
set(FAUST_DATA_TXT_DIR 		${FAUST_MISC_DIR}/data/txt CACHE INTERNAL "")
set(FAUST_CONFIG_DIR 		${FAUST_MISC_DIR}/config_file CACHE INTERNAL "")
set(FAUST_CONFIG_XML_DIR 	${FAUST_MISC_DIR}/config_file/xml CACHE INTERNAL "")


###### MISC DIRECTORIES ######
set(FAUST_SRC_TEST_DIR 	${FAUST_MISC_DIR}/test CACHE INTERNAL "")
set(FAUST_DATA_MAT_DIR 	${FAUST_MISC_DIR}/data/mat CACHE INTERNAL "")

set(FAUST_SRC_TEST_TIME_COMP_DIR 	${FAUST_SRC_TEST_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_SRC_TEST_TOOL_DIR 		${FAUST_SRC_TEST_DIR}/tools CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_DIR 			${FAUST_SRC_TEST_DIR}/src CACHE INTERNAL "")


###### WRAPPER MATLAB DIRECTORIES ######
if (BUILD_MATLAB_MEX_FILES)
	###### SRC DIRECTORY #####
	set(FAUST_MATLAB_SRC_DIR ${FAUST_SRC_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/tools CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/src CACHE INTERNAL "")
	set(FAUST_MATLAB_DOC_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/doc CACHE INTERNAL "")
	###### BIN DIRECTORY ######
	set(FAUST_MATLAB_BIN_DIR ${FAUST_BIN_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/mex CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/tools CACHE INTERNAL "")
	####  INSTALL DIRECTORY ####
	set(FAUST_MATLAB_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/matlab)
	set(FAUST_MATLAB_MEX_INSTALL_DIR 	${FAUST_MATLAB_INSTALL_DIR}/mex)
	set(FAUST_MATLAB_TOOLS_INSTALL_DIR 	${FAUST_MATLAB_INSTALL_DIR}/tools CACHE INTERNAL "")	


	########## DEMO DIRECTORIES #############
	###### SRC DIRECTORY #####
	set(FAUST_DEMO_SRC_DIR ${FAUST_MISC_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_SRC_DIR 		${FAUST_DEMO_SRC_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Brain_source_localization/data/ CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_SRC_DIR 	${FAUST_DEMO_SRC_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_TIMECOMP_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Runtime_comparison CACHE INTERNAL "")

	###### BIN DIRECTORY ######
	set(FAUST_DEMO_BIN_DIR 		${FAUST_MATLAB_BIN_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_BIN_DIR 		${FAUST_DEMO_BIN_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_BIN_DIR 	${FAUST_DEMO_BIN_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_BIN_DIR ${FAUST_DEMO_BSL_BIN_DIR}/data/ CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_OUTPUT_BIN_DIR 	${FAUST_DEMO_BSL_BIN_DIR}/output/ CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_OUTPUT_BIN_DIR 	${FAUST_DEMO_HADAMARD_BIN_DIR}/output/ CACHE INTERNAL "")
	set(FAUST_DEMO_TIMECOMP_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Runtime_comparison CACHE INTERNAL "")
	

	####  INSTALL DIRECTORY ####
	set(FAUST_DEMO_INSTALL_DIR ${FAUST_MATLAB_INSTALL_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_INSTALL_DIR 		${FAUST_DEMO_INSTALL_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_INSTALL_DIR 	${FAUST_DEMO_INSTALL_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_INSTALL_DIR ${FAUST_DEMO_INSTALL_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/data CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/data CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_OUTPUT_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/output CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_OUTPUT_INSTALL_DIR 	${FAUST_DEMO_HADAMARD_INSTALL_DIR}/output/ CACHE INTERNAL "")
	set(FAUST_DEMO_TIMECOMP_INSTALL_DIR ${FAUST_DEMO_INSTALL_DIR}/Runtime_comparison CACHE INTERNAL "")
endif(BUILD_MATLAB_MEX_FILES)


## examples directories including some examples of the use of the library FAuST
set(FAUST_BIN_TEST_DIR ${PROJECT_BINARY_DIR}/run_test CACHE INTERNAL "")
set(FAUST_BIN_TEST_BIN_DIR 		${FAUST_BIN_TEST_DIR}/bin CACHE INTERNAL "")
set(FAUST_BIN_TEST_SRC_DIR 		${FAUST_BIN_TEST_DIR}/src CACHE INTERNAL "")
set(FAUST_BIN_TEST_TMP_DIR 		${FAUST_BIN_TEST_DIR}/tmp CACHE INTERNAL "")
set(FAUST_BIN_TEST_OUTPUT_DIR 	${FAUST_BIN_TEST_DIR}/output CACHE INTERNAL "")
set(FAUST_BIN_TEST_FIG_DIR 		${FAUST_BIN_TEST_DIR}/fig CACHE INTERNAL "")
set(FAUST_BIN_TEST_TIME_COMPARE_DIR 	${FAUST_BIN_TEST_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_BIN_TEST_TOOLS_DIR 	${FAUST_BIN_TEST_DIR}/tools CACHE INTERNAL "")
##################################################################



### pour le GPU 
if (BUILD_USE_GPU)
	set(FAUST_SRC_LINEAR_OPERATOR_GPU_DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")
	set(FAUST_BIN_LINEAR_OPERATOR_GPU_DIR ${FAUST_BIN_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")

	set(FAUST_SRC_TEST_GPU_DIR ${FAUST_MISC_DIR}/test/GPU_test CACHE INTERNAL "")
	set(FAUST_BIN_TEST_GPU_DIR ${FAUST_BIN_TEST_DIR}/binGPU CACHE INTERNAL "") 
endif(BUILD_USE_GPU)


###### DOCUMENTATION DIRECTORIES ######
if (BUILD_DOCUMENTATION)
	set(FAUST_DOC_SRC_DIR ${PROJECT_SOURCE_DIR}/gen_doc CACHE INTERNAL "")
	set(FAUST_DOC_BIN_DIR ${PROJECT_BINARY_DIR}/doc CACHE INTERNAL "")
endif()
##################################################################

#message(STATUS "FAUST_BIN_TEST_OUTPUT_DIR : ${FAUST_BIN_TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_OUTPUT_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_FIG_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_TMP_DIR})

if (BUILD_MATLAB_MEX_FILES)
	file(MAKE_DIRECTORY ${FAUST_DEMO_BSL_OUTPUT_BIN_DIR})
	install(DIRECTORY ${FAUST_DEMO_BSL_OUTPUT_BIN_DIR} DESTINATION ${FAUST_DEMO_BSL_OUTPUT_INSTALL_DIR})
	file(MAKE_DIRECTORY ${FAUST_DEMO_HADAMARD_OUTPUT_BIN_DIR})
	install(DIRECTORY ${FAUST_DEMO_HADAMARD_OUTPUT_BIN_DIR} DESTINATION ${FAUST_DEMO_HADAMARD_OUTPUT_INSTALL_DIR})
endif (BUILD_MATLAB_MEX_FILES)

# CANNOT CREATE EMPTY DIRECTORY IN INSTALL_DIR  !!!!! MODIFY !!!!!
#INSTALL(CODE "FILE(MAKE_DIRECTORY /logs)")

###### INSTALL DIRECTORIES (RELATIVES PATHS FROM CMAKE_INSTALL_PREFIX) ######
#set(FAUST_INSTALL_MEX ${CMAKE_INSTALL_PREFIX}/mex)
#set(FAUST_INSTALL_LIB ${CMAKE_INSTALL_PREFIX}/lib)
#set(FAUST_INSTALL_CMDLINE ${CMAKE_INSTALL_PREFIX}/CmdLine)
#set(FAUST_INSTALL_BIN ${CMAKE_INSTALL_PREFIX}/testing/bin)
#set(FAUST_INSTALL_TESTING_DATA ${CMAKE_INSTALL_PREFIX}/testing/data)
#set(FAUST_INSTALL_BIN_TIME_COMP ${CMAKE_INSTALL_PREFIX}/testing/runtime_comparison CACHE INTERNAL "")

#message(STATUS "ALALALLLLLLLLLLL CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}")
#message(STATUS "ALALALLLLLLLLLLL FAUST_INSTALL_BIN=${FAUST_INSTALL_BIN}")
#message(STATUS "ALALALLLLLLLLLLL FAUST_INSTALL_LIB=${FAUST_INSTALL_LIB}")


###### TARGETS ######
set(FAUST_LINEAR_OPERATOR_TARGET faust_linear_operator CACHE INTERNAL "")
set(FAUST_ALGORITHM_TARGET faust_algorithm CACHE INTERNAL "")
set(FAUST_MEX_TARGET faust_mex CACHE INTERNAL "")
set(FAUST_MEXTOOLS_TARGET faust_mextools CACHE INTERNAL "")
set(FAUST_MATIO_TARGET faust_matio CACHE INTERNAL "")
set(FAUST_XML_TARGET faust_xml CACHE INTERNAL "")
set(FAUST_TARGET faust CACHE INTERNAL "")
if (BUILD_USE_GPU)
	set(FAUST_LINEAR_OPERATOR_GPU_TARGET faust_linear_operator_gpu CACHE INTERNAL "")
	#set(FAUST_MATRIX_CU_TARGET faust_gpu_matrix CACHE INTERNAL "")
	#set(FAUST_FAUSTCORE_CU_TARGET faust_gpu_core CACHE INTERNAL "")
	set(FAUST_ALGORITHM_CU_TARGET faust_algorithm_gpu CACHE INTERNAL "")
endif(BUILD_USE_GPU)
##################################################################


###### FLAG for DEBUG or RELEASE TARGETS ######
if(WIN32)
	if(CMAKE_BUILD_TYPE MATCHES "Debug")
		set(FAUST_DEBUG_RELEASE "debug")
		set(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/Debug)
		set(FAUST_MATLAB_BIN_DIR ${FAUST_MATLAB_BIN_DIR}/Debug)
	else()
		set(FAUST_DEBUG_RELEASE "release")
		set(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/Release)
		set(FAUST_MATLAB_BIN_DIR ${FAUST_MATLAB_BIN_DIR}/Release)
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
		if (BUILD_MULTITHREAD)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
		endif(BUILD_MULTITHREAD)
					
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -pthread -fPIC -O2 " CACHE STRING "compile flags" FORCE)
		
		
	elseif(WIN32)
		if (BUILD_MULTITHREAD)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp " CACHE STRING "compile flags" FORCE)
		endif(BUILD_MULTITHREAD)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /EHs /MD /O2" CACHE STRING "compile flags" FORCE)
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
	endif()
	message(STATUS "**********RELEASE mode************")
endif()

#if(UNIX)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp -fexceptions -fno-omit-frame-pointer -pthread -fPIC " CACHE STRING "compile flags" FORCE)
#elseif(WIN32)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /EHs /openmp /MD" CACHE STRING "compile flags" FORCE)
#else()
#	message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
#endif(UNIX)

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
##################################################################


#########################################
###### FIND EXTERNALS LIBRARY  ##########
include(CMake/findExternalLibs.cmake)
#########################################
#########################################

########### add_subdirectory ###############
add_subdirectory(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_BIN_LINEAR_OPERATOR_DIR})
add_subdirectory(${FAUST_ALGORITHM_SRC_DIR} ${FAUST_ALGORITHM_BIN_DIR})
add_subdirectory(${FAUST_SRC_TEST_DIR} ${FAUST_BIN_TEST_BIN_DIR})
if (BUILD_DOCUMENTATION)
	add_subdirectory(${FAUST_DOC_SRC_DIR} ${FAUST_DOC_BIN_DIR})
endif ()
add_subdirectory(${FAUST_SRC_CMDLINE_DIR} ${FAUST_BIN_CMDLINE_BIN_DIR})
if (BUILD_USE_GPU)
	add_subdirectory(${FAUST_SRC_LINEAR_OPERATOR_GPU_DIR} ${FAUST_BIN_LINEAR_OPERATOR_GPU_DIR})
	#add_subdirectory(${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_FAUSTCORE_CU_BIN_DIR})
	#add_subdirectory(${FAUST_ALGORITHM_CU_SRC_DIR} ${FAUST_ALGORITHM_CU_BIN_DIR})
	add_subdirectory(${FAUST_SRC_TEST_GPU_DIR} ${FAUST_BIN_TEST_GPU_DIR})
endif (BUILD_USE_GPU)

if (BUILD_MATLAB_MEX_FILES)
	add_subdirectory(${FAUST_MATLAB_SRC_DIR} ${FAUST_MATLAB_MEX_BIN_DIR})
	add_subdirectory(${FAUST_DEMO_SRC_DIR} ${FAUST_DEMO_BIN_DIR})
endif()
##################################################################


###### TARGET GENERATION ######
if (BUILD_USE_GPU)
	if (BUILD_MATLAB_MEX_FILES)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_CU_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (BUILD_MATLAB_MEX_FILES)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_CU_TARGET}>)
	endif(BUILD_MATLAB_MEX_FILES)
else (BUILD_USE_GPU)
	if (BUILD_MATLAB_MEX_FILES)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (BUILD_MATLAB_MEX_FILES)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_TARGET}>)
	endif(BUILD_MATLAB_MEX_FILES)
endif (BUILD_USE_GPU)


if (BUILD_MATLAB_MEX_FILES)
	if (UNIX) 
		if(APPLE)
			target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}/libmex.dylib)
		else(APPLE)
			target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}/libmex.so)
		endif(APPLE)
		
	elseif(WIN32)
		message(FATAL_ERROR "MEXFILE : CANNOT SEARCH LIBMEX IN WINDOWS")
	endif()
endif(BUILD_MATLAB_MEX_FILES)

# In case of GPU, the "CUDA" library libfaust_matrix_cu.a is linked with faust (libfaust_matrix_cu.a is not a OBJECT library) 
if (BUILD_USE_GPU)
	#target_link_libraries(${FAUST_TARGET} ${FAUST_MATRIX_CU_TARGET} )	
	target_link_libraries(${FAUST_TARGET} ${FAUST_LINEAR_OPERATOR_GPU_TARGET} )	
endif(BUILD_USE_GPU)


if(BUILD_OPENBLAS)
	target_link_libraries(${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
else()
	target_link_libraries(${FAUST_TARGET} )	
endif(BUILD_OPENBLAS)

install(TARGETS ${FAUST_TARGET} DESTINATION ${FAUST_INSTALL_LIB})
##################################################################

#add_dependencies(${FAUST_FAUSTCORE_TARGET} ${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})
#add_dependencies(${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})



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



