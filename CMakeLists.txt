##############################################################################
##                              Description:                                ##
##  Main cmakelist file for the faust project                               ##
##  Other important files for CMake are                                     ##
##      - CMakeLists.txt files in the sub-directories                       ##
##      - CMake/* files, which define some tools to search for libraries.   ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##############################################################################

#cmake_minimum_required(VERSION 3.2.2) #compatibilite UTF-8 ASCII (Windows Ubunt)u)	
cmake_minimum_required(VERSION 3.0.2) #compatibilite UTF-8 ASCII (Windows Ubunt)u)
#cmake_minimum_required(VERSION 2.8.8)
project(FAUST CXX)

##############    FAUST_VERSION 2.1)  ############################
set(FAUST_VERSION_MAJOR 2)
set(FAUST_VERSION_MINOR 1)

set(FAUST_VERSION_PATCH 0)
########################################################

##############  INSTALL DIRECTORIES   #############################
if(UNIX)
	if(APPLE)
		message(STATUS "APPLE OPERATING SYSTEM")
	else(APPLE)
		message(STATUS "LINUX OPERATING SYSTEM")
	endif(APPLE)

#cmake .. -DCMAKE_INSTALL_PREFIX="home/project/Faust_install" 
#	if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#		SET(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/build/FAuST_OUTPUT" CACHE STRING "default install path" FORCE )
#	else()
#		SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
#		message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
#	endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

message (STATUS "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}") 
	if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

		# default location of library is in the path  		
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}/lib" CACHE STRING "default install path for the library faust" FORCE )

#		set(FAUST_INSTALL_BIN 			"$ENV{HOME}/Documents/faust/bin" CACHE STRING "default install path for the binary faust" FORCE )
#		set(FAUST_INSTALL_BIN_TIME_COMP "$ENV{HOME}/Documents/faust/bin/runtime_comparison" CACHE STRING "default install path for the binary faust" FORCE )
#		set(FAUST_INSTALL_DOC 			"$ENV{HOME}/Documents/faust/doc" CACHE STRING "default install path for the documentation faust" FORCE )
#		set(FAUST_INSTALL_TESTING_DATA 	"$ENV{HOME}/Documents/faust/testing/data" CACHE STRING "default install path for the test faust" FORCE )

#		set(FAUST_INSTALL_BIN 			"${PROJECT_BINARY_DIR}/output/bin" CACHE STRING "default install path for the binary faust" FORCE )
#		set(FAUST_INSTALL_BIN_TIME_COMP "${PROJECT_BINARY_DIR}/output/bin/runtime_comparison" CACHE STRING "default install path for the binary faust" FORCE )

		#set(FAUST_INSTALL_DOC 			"${PROJECT_BINARY_DIR}/output/doc" CACHE STRING "default install path for the documentation faust" FORCE )
		#set(FAUST_INSTALL_TESTING_DATA 	"${PROJECT_BINARY_DIR}/output/testing/data" CACHE STRING "default install path for the test faust" FORCE )
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "$ENV{HOME}/Documents/MATLAB/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		SET(CMAKE_INSTALL_PYTHON_PREFIX "$ENV{HOME}/Documents/PYTHON/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		# CMAKE_INSTALL_PREFIX path is deleted for next cmake .. configure
		SET(CMAKE_INSTALL_PREFIX "" CACHE STRING "default install path" FORCE )
		#message (STATUS "default install") 
	# Install in the CMAKE_INSTALL_PREFIX path defined by the user (local install) 	
	elseif(CMAKE_INSTALL_PREFIX)
		SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
		message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}/lib")

#		set(FAUST_INSTALL_BIN 			"${CMAKE_INSTALL_PREFIX}/bin")
#		set(FAUST_INSTALL_BIN_TIME_COMP "${CMAKE_INSTALL_PREFIX}/bin/runtime_comparison")
		#set(FAUST_INSTALL_DOC 			"${CMAKE_INSTALL_PREFIX}/doc")
		#set(FAUST_INSTALL_TESTING_DATA 	"${CMAKE_INSTALL_PREFIX}/testing/data")
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "${CMAKE_INSTALL_PREFIX}/matlab")
		SET(CMAKE_INSTALL_PYTHON_PREFIX "${CMAKE_INSTALL_PREFIX}/python")
	else()
		message(STATUS "Install FAuST path already defined to : ${CMAKE_INSTALL_MATLAB_PREFIX}")

	endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


else(UNIX)
	message(STATUS "WINDOWS OPERATING SYSTEM")

message (STATUS "CMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}") 
	if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

		# default location of library is in the path  		
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}" CACHE STRING "default install path for the library faust" FORCE )

#		set(FAUST_INSTALL_BIN 			"$ENV{HOME}/Documents/faust/bin" CACHE STRING "default install path for the binary faust" FORCE )
#		set(FAUST_INSTALL_BIN 			"${PROJECT_BINARY_DIR}/output/bin" CACHE STRING "default install path for the binary faust" FORCE )
#		set(FAUST_INSTALL_BIN_TIME_COMP "${PROJECT_BINARY_DIR}/output/bin/runtime_comparison" CACHE STRING "default install path for the binary faust" FORCE )
		#set(FAUST_INSTALL_DOC 			"${PROJECT_BINARY_DIR}/output/doc" CACHE STRING "default install path for the documentation faust" FORCE )
		#set(FAUST_INSTALL_TESTING_DATA 	"${PROJECT_BINARY_DIR}/output/testing/data" CACHE STRING "default install path for the test faust" FORCE )
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "C:/Users/$ENV{USERNAME}/Documents/MATLAB/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		SET(CMAKE_INSTALL_MATLAB_PREFIX "C:/Users/$ENV{USERNAME}/Documents/PYTHON/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		# CMAKE_INSTALL_PREFIX path is deleted for next cmake .. configure
		SET(CMAKE_INSTALL_PREFIX "" CACHE STRING "default install path" FORCE )
		#message (STATUS "default install") 
	# Install in the CMAKE_INSTALL_PREFIX path defined by the user (local install) 	
	elseif(CMAKE_INSTALL_PREFIX)
		SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
		message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
		
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}/lib")
#		set(FAUST_INSTALL_BIN 			"${CMAKE_INSTALL_PREFIX}/bin")
#		set(FAUST_INSTALL_BIN_TIME_COMP "${CMAKE_INSTALL_PREFIX}/bin/runtime_comparison")
		#set(FAUST_INSTALL_DOC 			"${CMAKE_INSTALL_PREFIX}/doc")
		#set(FAUST_INSTALL_TESTING_DATA 	"${CMAKE_INSTALL_PREFIX}/testing/data")
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "${CMAKE_INSTALL_PREFIX}/matlab")
		SET(CMAKE_INSTALL_PYTHON_PREFIX "${CMAKE_INSTALL_PREFIX}/python")
	else()
		message(STATUS "Install FAuST path already defined to : ${CMAKE_INSTALL_MATLAB_PREFIX}")

	endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

endif(UNIX)

##################################################################
## For windows platform :
# cmake -G "MinGW Makefiles" ..
# message(STATUS "CXX=$ENV{CXX}")
# sous windows on utilise MinGW pour le compilateur GCC G++ au lieu de visual studio plous lourd à installer. MinGW version 4.9.2 pour les mex functions..
# set(CTEST_CMAKE_GENERATOR "MinGW Makefiles") 


###### tmp directory where temporary objects will be located ######
set(FAUST_TMP_BUILD_DIR "${PROJECT_BINARY_DIR}/tmp")
file(MAKE_DIRECTORY ${FAUST_TMP_BUILD_DIR})
set(FAUST_TMPMEX_DIR "${FAUST_TMP_BUILD_DIR}/mex_obj")
file(MAKE_DIRECTORY ${FAUST_TMPMEX_DIR})

if (MATLAB_EXE_DIR)
# nothing 
else(MATLAB_EXE_DIR)
	SET(MATLAB_EXE_DIR " " CACHE STRING "Force the directory of your expected matlab binary (default Value: None)" FORCE )
endif(MATLAB_EXE_DIR)


if(UNIX)
	set(MEXOBJ_EXT "o")
elseif(WIN32)
	set(MEXOBJ_EXT "obj")
endif()



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
option(BUILD_DEBUG "Enable FAUST Debug mode." OFF)
# Allows Matlab MEX files :

option(BUILD_WRAPPER_MATLAB "Enable building Matlab MEX files." OFF)

# different type of files 
option(BUILD_READ_XML_FILE "Using xml configuration to read xml files" OFF)
option(BUILD_READ_MAT_FILE "Using matio library to read mat files" ON)
#set(BUILD_READ_XML_FILE OFF)
#set(BUILD_READ_MAT_FILE OFF)

option(BUILD_WRAPPER_PYTHON "Enable building PyFaust (Python librairie)." OFF)
#set(BUILD_WRAPPER_PYTHON OFF)
option(USE_MATIO_STATIC_LIBS OFF) # if ON the script waits several variables to
# be set: MATIO_STATIC_LIB_PATH, Z_STATIC_LIB_PATH, HDF5_STATIC_LIB_PATH
#CMAKE_DEPENDENT_OPTION


###### Different option of configuration ######
#option(BUILD_TESTING "Build Tests" ON) 
set(BUILD_USE_SINGLEPRECISION OFF) # CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
#set(BUILD_USE_PROFILING OFF) # CACHE BOOL "Profiling the code")
option(BUILD_USE_PROFILING "Profiling the code" OFF) # CACHE BOOL )

option(SLOW_TESTS "Enabling tests which take much time." ON)
option(EXCLUDE_FAUST_LIB_INSTALL "Excluding faust static core library from installation files." OFF)
option(INCLUDE_ALL_SYS_INSTALL_INSTRUCS "In generated doc/html/README.html put all systems instructions for installation, not just the building system ones." OFF)
option(API_DOC_BASE_URL "Base URL for the online or local documentation.") # see gen_doc/CMakeLists.txt
option(REMOTE_DATA_URL "Base URL for getting FAµST data archive (filename excluded from URL, see variable REMOTE_DATA_FILE)." "")
option(REMOTE_DATA_FILE "Filename of the FAµST data zip archive." "")
option(NOCPPTESTS "Enabling/Disabling CPP tests (only testing python or matlab wrappers if built)." OFF)
option(EXPERIMENTAL_PKG "Enabling/Disabling experimental code." OFF)
# experimental code is delimited for example in matlab scripts by comments like "% experimental block start / end"
option(NOPY2 "Enabling/Disabling python2 support." ON)
option(PYTHON_ENCODING "The *.py charset for pyfaust/__init__.py (utf-8, windows-1252)." utf-8)
#set(BUILD_MULTITHREAD OFF CACHE BOOL "OpenMP Multithreading")
#set(BUILD_VERBOSE OFF CACHE BOOL "Display useful message for debugging")
#set(BUILD_DOCUMENTATION OFF CACHE BOOL "Generate html documentation with doxygen")
#set(BUILD_READ_XML_FILE OFF CACHE BOOL "Using xml configuration to read xml files")
#set(BUILD_WRAPPER_MATLAB ON CACHE BOOL "Generate Mexfiles")
#set(BUILD_READ_MAT_FILE OFF CACHE BOOL "Using matio library to read mat files")

option(FAUST_TORCH "Faust torch backend for the Faust-Matrix/Vector product." OFF)
option(PYTORCH_EXP "pytorch Faust class experimental code." OFF)
option(USE_GPU_MOD "Faust cuda backend for the Faust-Matrix/Vector product." OFF)
option(USE_OPENBLASO "Enable the use of OpenBLAS-OpenMP eigen backend in pyfaust (library embedded in the pip package)." OFF)
option(BUILD_COMPLEX_PYX "Building the python wrapper in order to support complex matrices or not." ON)
option(BUILD_FLOAT_PYX "Enabling the support of float matrices in the python wrapper." OFF)
option(BUILD_FLOAT_MEX "Enabling the support of float matrices in the matlab wrappers." OFF)
option(DISABLE_ALGO_PYX "If OFF the algorithms (pyfaust.fact, pyfaust.poly etc.) python wrappers are built otherwise they're not (useful to avoid large compilation time when updating only the Faust class)." OFF)
option(BUILD_EIGTJ_SVDTJ "if OFF the eigjt and svdtj C++ implementations are not built." ON)
option(MATLAB_2017_API "if ON the deprecated Matlab API will be used to build -- including separate complex API instead of new interleaved complex (R2018a)" OFF)
option(VCOMPLIB_PATH "Set the library path of OpenMP which comes from the Visual Studio environment, e.g. C:/Program Files (x86)/Microsoft Visual Studio/*/vcomp140.dll")

if (BUILD_USE_SINGLEPRECISION)
	message(STATUS "**********SINGLEPRECISION ACTIF************")
endif()
if (BUILD_WRAPPER_MATLAB)
	message(STATUS "**********MEX ACTIF************")
	message(STATUS "CMAKE_INSTALL_MATLAB_PREFIX = ${CMAKE_INSTALL_MATLAB_PREFIX}")
endif()
if (BUILD_USE_PROFILING)
	message(STATUS "**********PROFILING ACTIF************")
endif()
if (BUILD_READ_MAT_FILE)
	message(STATUS "**********MATIO ACTIF************")
	#message(STATUS "FAUST_INSTALL_BIN = ${FAUST_INSTALL_BIN}")
endif()
if (BUILD_READ_XML_FILE)
	message(STATUS "**********XML ACTIF************")
endif()
if (BUILD_VERBOSE)
	message(STATUS "**********VERBOSE ACTIF************")
endif()
if (USE_GPU_MOD)
	message(STATUS "**********GPU MODULE ENABLED************")
endif()
if(BUILD_MULTITHREAD)
	message(STATUS "**********OPENMP ACTIF*********")
endif()
if (BUILD_DOCUMENTATION)
	message(STATUS "**********DOC ACTIF************")
endif()
if(BUILD_WRAPPER_PYTHON)
	message(STATUS "**********PYTHON ACTIF************")
endif(BUILD_WRAPPER_PYTHON)
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

set(FAUST_UTILS_SRC_DIR ${FAUST_SRC_DIR}/utils CACHE INTERNAL "")
set(FAUST_UTILS_BIN_DIR ${PROJECT_BINARY_DIR}/utils CACHE INTERNAL "")

#set(FAUST_EXCEPTION_SRC_DIR ${FAUST_SRC_DIR}/error_display CACHE INTERNAL "")
#set(FAUST_FAUSTCORE_SRC_DIR ${PROJECT_SOURCE_DIR}/faust_core CACHE INTERNAL "")


###### WRAPPER DIRECTORIES ######
set(FAUST_SRC_WRAPPER_DIR ${PROJECT_SOURCE_DIR}/wrapper CACHE INTERNAL "")
set(FAUST_BIN_WRAPPER_DIR ${PROJECT_BINARY_DIR}/wrapper CACHE INTERNAL "")
##################################################################

## fichier d'entrée doivent être renseigné dans l'executable
###### WRAPPER DATA DIRECTORIES ######
set(FAUST_MISC_DIR ${PROJECT_SOURCE_DIR}/misc CACHE INTERNAL "")
set(FAUST_DATA_DIR 			${FAUST_MISC_DIR}/data CACHE INTERNAL "")
set(FAUST_DATA_TXT_DIR 		${FAUST_MISC_DIR}/data/txt CACHE INTERNAL "")
set(FAUST_CONFIG_DIR 		${FAUST_MISC_DIR}/config_file CACHE INTERNAL "")
set(FAUST_CONFIG_XML_DIR 	${FAUST_MISC_DIR}/config_file/xml CACHE INTERNAL "")
set(FAUST_CONFIG_MAT_DIR 	${FAUST_MISC_DIR}/config_file/mat CACHE INTERNAL "")


###### MISC DIRECTORIES ######
set(FAUST_SRC_TEST_DIR 	${FAUST_MISC_DIR}/test CACHE INTERNAL "")
set(FAUST_DATA_MAT_DIR 	${FAUST_MISC_DIR}/data/mat CACHE INTERNAL "")

set(FAUST_SRC_TEST_TIME_COMP_DIR 	${FAUST_SRC_TEST_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_SRC_TEST_TOOL_DIR 		${FAUST_SRC_TEST_DIR}/tools CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_DIR 			${FAUST_SRC_TEST_DIR}/src CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_CPP_DIR ${FAUST_SRC_TEST_SRC_DIR}/C++/ CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_UNIT_DIR ${FAUST_SRC_TEST_SRC_CPP_DIR}/unit/ CACHE INTERNAL "")	
set(FAUST_SRC_TEST_SRC_MATLAB_DIR ${FAUST_SRC_TEST_SRC_DIR}/Matlab/ CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_PYTHON_DIR ${FAUST_SRC_TEST_SRC_DIR}/Python/ CACHE INTERNAL "")


###### WRAPPER MATLAB DIRECTORIES ######
if (BUILD_WRAPPER_MATLAB)
	#version MATLAB du code (entierement MATLAB, pas de mexfunction)
	#file(MAKE_DIRECTORY ${CMAKE_INSTALL_MATLAB_PREFIX})
	###### SRC DIRECTORY #####
	set(FAUST_MATLAB_SRC_DIR ${FAUST_SRC_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/tools CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/src CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_MFILES_SRC_DIR ${FAUST_MATLAB_SRC_DIR}/mex CACHE INTERNAL "")
	set(FAUST_MATLAB_DOC_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/doc CACHE INTERNAL "")
	###### BIN DIRECTORY ######
	set(FAUST_MATLAB_BIN_DIR ${FAUST_BIN_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/mex CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/tools CACHE INTERNAL "")
	set(FAUST_MATLAB_NAMESPACE_REL_DIR "+matfaust")
	set(FAUST_MATLAB_NAMESPACE_SRC_DIR	${FAUST_MATLAB_SRC_DIR}/${FAUST_MATLAB_NAMESPACE_REL_DIR})
	set(FAUST_MATLAB_NAMESPACE_BIN_DIR	${FAUST_MATLAB_BIN_DIR}/${FAUST_MATLAB_NAMESPACE_REL_DIR})
	file(GLOB_RECURSE FAUST_MATLAB_MATFAUST_FILE_LIST ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/*.m)
	####  INSTALL DIRECTORY ####
	#set(FAUST_MATLAB_INSTALL_DIR ${CMAKE_INSTALL_MATLAB_PREFIX})

	#set(FAUST_MATLAB_MEX_INSTALL_DIR 	${CMAKE_INSTALL_MATLAB_PREFIX}/mex)
	#set(FAUST_MATLAB_TOOLS_INSTALL_DIR 	${CMAKE_INSTALL_MATLAB_PREFIX}/tools CACHE INTERNAL "")	


	###### DEMO DATA BIN DIRECTORY ######
	set(FAUST_MATFAUST_DEMO_DATA_BIN_DIR ${FAUST_MATLAB_BIN_DIR}/data CACHE INTERNAL "")



endif(BUILD_WRAPPER_MATLAB)
##################################################################





if(BUILD_WRAPPER_PYTHON)
	set(FAUST_PYTHON_SRC_DIR ${FAUST_SRC_WRAPPER_DIR}/python CACHE INTERNAL "")
	set(FAUST_PYTHON_SRC_SRC_DIR ${FAUST_PYTHON_SRC_DIR}/src/ CACHE INTERNAL "")
	set(FAUST_PYTHON_BIN_DIR ${FAUST_BIN_WRAPPER_DIR}/python CACHE INTERNAL "")
	set(FAUST_PYTHON_TARGET faust_python CACHE INTERNAL "")
	set(PYFAUST_PYTHON_REQUIREMENTS "pygsp scipy numpy>=1.24 matplotlib>=2.0.0") # pathlib") #ENOTE: order matters in setup.py, the last package is installed first
endif(BUILD_WRAPPER_PYTHON)


set(FAUST_BIN_TEST_DIR ${PROJECT_BINARY_DIR}/run_test CACHE INTERNAL "")
set(FAUST_BIN_TEST_BIN_DIR 		${FAUST_BIN_TEST_DIR}/bin CACHE INTERNAL "")
set(FAUST_BIN_TEST_SRC_DIR 		${FAUST_BIN_TEST_DIR}/src CACHE INTERNAL "")
set(FAUST_BIN_TEST_TMP_DIR 		${FAUST_BIN_TEST_DIR}/tmp CACHE INTERNAL "")
set(FAUST_BIN_TEST_OUTPUT_DIR 	${FAUST_BIN_TEST_DIR}/output CACHE INTERNAL "")
set(FAUST_BIN_TEST_FIG_DIR 		${FAUST_BIN_TEST_DIR}/fig CACHE INTERNAL "")
set(FAUST_BIN_TEST_TOOLS_DIR 	${FAUST_BIN_TEST_DIR}/tools CACHE INTERNAL "")


##################################################################

if (USE_GPU_MOD)
	set(FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU2 CACHE INTERNAL "")
	set(FAUST_BIN_LINEAR_OPERATOR_GPU2_DIR ${FAUST_BIN_LINEAR_OPERATOR_DIR}/GPU2 CACHE INTERNAL "")
endif(USE_GPU_MOD)
##################################################################

###### DOCUMENTATION DIRECTORIES ######
if (BUILD_DOCUMENTATION)
	#file(MAKE_DIRECTORY ${FAUST_INSTALL_DOC})
	set(FAUST_DOC_SRC_DIR ${PROJECT_SOURCE_DIR}/gen_doc CACHE INTERNAL "")
	set(FAUST_DOC_BIN_DIR ${PROJECT_BINARY_DIR}/doc CACHE INTERNAL "")
endif()
##################################################################
#message(STATUS "FAUST_BIN_TEST_OUTPUT_DIR : ${FAUST_BIN_TEST_OUTPUT_DIR}")

file(MAKE_DIRECTORY ${FAUST_BIN_TEST_OUTPUT_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_FIG_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_TMP_DIR})



#       [FILE_PERMISSIONS permissions...]
#       [DIRECTORY_PERMISSIONS permissions...]
set(INSTALL_DIR_PERMISSION OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
set(INSTALL_FILE_PERMISSION OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)


###### TARGETS ######
set(FAUST_LINEAR_OPERATOR_TARGET faust_linear_operator CACHE INTERNAL "")
set(FAUST_ALGORITHM_TARGET faust_algorithm CACHE INTERNAL "")
set(FAUST_MEX_TARGET faust_mex CACHE INTERNAL "")
set(FAUST_MEXTOOLS_TARGET faust_mextools CACHE INTERNAL "")
set(FAUST_MATIO_TARGET faust_matio CACHE INTERNAL "")
set(FAUST_XML_TARGET faust_xml CACHE INTERNAL "")
set(FAUST_TARGET faust CACHE INTERNAL "")
##################################################################


######################################"
# Case of Visual Studio 
IF(WIN32)
    IF (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
		# Case of debug mode
		#SET(FAUST_DEBUG_RELEASE "debug")
		SET(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/lib/debug)
    ELSE (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
    	# Case of release mode
		#SET(FAUST_DEBUG_RELEASE "release")
		SET(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/lib/release)
    ENDIF (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )

ELSE(WIN32)
	# Case of Xcode
	IF( CMAKE_GENERATOR STREQUAL "Xcode")
    	# We would like to select a single build mode, would we ?
	    # or to select it dynamically there is the CMAKE_CFG_INTDIR stuff
	    IF (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
    		SET(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/lib/Debug)
	    ELSE (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
    		SET(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/lib/Release)
#file(MAKE_DIRECTORY ${FAUST_BIN_DIR})

	    ENDIF (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
	ELSE( CMAKE_GENERATOR STREQUAL "Xcode")
	    # Generic case
	    SET(FAUST_BIN_DIR ${PROJECT_BINARY_DIR}/lib)
	ENDIF(CMAKE_GENERATOR STREQUAL "Xcode") 
    #SET(MPTK_LIB_DIR_INSTALL ${CMAKE_INSTALL_PREFIX}/mptk/plugins)
ENDIF(WIN32)

#message(STATUS "FAUST_BIN_DIR : ${FAUST_BIN_DIR}")
#message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY : ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${FAUST_BIN_DIR} CACHE STRING "" FORCE)
#message(STATUS "CMAKE_ARCHIVE_OUTPUT_DIRECTORY : ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${FAUST_BIN_DIR} CACHE STRING "" FORCE)
############################################################################



## adding compiler flags ###################################################
set(CMAKE_CXX_FLAGS "")
#Add the c++11 flag, whatever it is
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-std=c++11 COMPILER_SUPPORTS_CXX11)
#set(CMAKE_CXX_FLAGS)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" CACHE STRING "compile flags" FORCE)
	message(STATUS "C++ compiler supprt c++11")
else()
	message(STATUS "Failed to determine the compiler ${CMAKE_CXX_COMPILER} C++11 flag.")
endif()



if(BUILD_MULTITHREAD)
	find_package(OpenMP) #not REQUIRED, allow to fail the first time
	if(NOT OpenMP_FOUND) # this section is particularly useful for clang/macOS configuration (OpenMP is not enabled by default, installed it via macports)
		#  message(STATUS "OpenMP not found automatically, let's help it.")
		if(DEFINED ENV{OpenMP_INC_DIR})
			set(OpenMP_INC_DIR $ENV{OpenMP_INC_DIR})
		endif()
		if(DEFINED ENV{OpenMP_gomp_LIBRARY})
			set(OpenMP_gomp_LIBRARY $ENV{OpenMP_gomp_LIBRARY})
		endif()
		if(DEFINED OpenMP_INC_DIR AND DEFINED OpenMP_gomp_LIBRARY)
			set(OpenMP_CXX_FLAGS "-I${OpenMP_INC_DIR}")
			set(OpenMP_CXX_LIB_NAMES "gomp")
		else()
			message(FATAL_ERROR "OpenMP not found, please set cmake variables OpenMP_INC_DIR to omp.h parent directory and OpenMP_gomp_LIBRARY to libgomp.* library path")
		endif()
		# now retry
		find_package(OpenMP) # last chance, don't use REQUIRED option here too, in order to print the fatal error with lead to fix
	endif()
	message(STATUS "OpenMP is supported OpenMP_CXX_FOUND:  ${OpenMP_CXX_FOUND}, OpenMP_CXX_FLAGS: ${OpenMP_CXX_FLAGS}, OpenMP_CXX_LIB_NAMES: ${OpenMP_CXX_LIB_NAMES}, OpenMP_CXX_LIBRARY: ${OpenMP_CXX_LIBRARY}, OpenMP_CXX_LIBRARIES: ${OpenMP_CXX_LIBRARIES}")
	foreach(OMP_FLAG IN LISTS OpenMP_CXX_FLAGS)
		set(OMP_CXX_FLAGS "${OMP_CXX_FLAGS} ${OMP_FLAG}")
	endforeach()
	if(UNIX OR (WIN32 AND ${CMAKE_GENERATOR} MATCHES "MinGW Makefiles"))
		foreach(OMP_LFLAG IN LISTS OpenMP_CXX_LIB_NAMES)
			set(OMP_CXX_LFLAGS "${OMP_CXX_LFLAGS} -l${OMP_LFLAG}")
		endforeach()
		foreach(OMP_LFLAG IN LISTS OpenMP_CXX_LIBRARIES)
			string(REGEX REPLACE "(.*)(/|\\\\).*" "\\1" OMP_LIB_PATH ${OMP_LFLAG})
			set(OMP_CXX_LFLAGS " ${OMP_CXX_LFLAGS} -L${OMP_LIB_PATH}")
		endforeach()
	elseif(WIN32 AND ${CMAKE_GENERATOR} MATCHES "Visual Studio.*")
		# strangely I noticed that findOpenMP gives -openmp instead of /openmp when compiler is msvc cl.exe
		string(REGEX REPLACE "-" "/" OMP_CXX_FLAGS ${OMP_CXX_FLAGS})
	endif()
	string(REGEX REPLACE "(^[ \t\r\n]+)|([ \t\r\n]+$)" "" OMP_CXX_FLAGS ${OMP_CXX_FLAGS})
	if(DEFINED OMP_CXX_LFLAGS)
		string(REGEX REPLACE "(^[ \t\r\n]+)|([ \t\r\n]+$)" "" OMP_CXX_LFLAGS ${OMP_CXX_LFLAGS})
	endif()
	message(STATUS OMP_CXX_FLAGS=${OMP_CXX_FLAGS} OMP_CXX_LFLAGS=${OMP_CXX_LFLAGS})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OMP_CXX_FLAGS}")
	set(OMP_ENABLED ON) # (CPP preprocessor compilation constant) used in cpp code to detect if OMP is usable
endif(BUILD_MULTITHREAD)
configure_file(${FAUST_SRC_DIR}/utils/faust_openmp.h.in ${FAUST_SRC_DIR}/utils/faust_openmp.h)

#if( (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
if (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
	if (UNIX)	
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} -g -pg -std=c++11")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -fPIC " CACHE STRING "compile flags" FORCE)
	elseif(WIN32)
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} /O0")
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")
	endif() 
	message(STATUS "**********DEBUG mode************")
else() # RELEASE CONFIG
	if (UNIX)	
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -fPIC -O3" CACHE STRING "compile flags" FORCE)
	elseif(WIN32)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHs /MD /O2" CACHE STRING "compile flags" FORCE)
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")
	endif()
	message(STATUS "**********RELEASE mode************")
endif(BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
############################################################################

###### MATIO inc and lib env. vars #########
if(NOT "$ENV{MATIO_INC_DIR}" MATCHES "^$")
    message(STATUS "MATIO_INC_DIR=$ENV{MATIO_INC_DIR}")
    set(MATIO_INC_DIR "$ENV{MATIO_INC_DIR}")
    set(MATIO_LIB_FILE "$ENV{MATIO_LIB_FILE}")
endif()
############################################


#########################################
###### FIND EXTERNALS LIBRARY  ##########
include(CMake/findExternalLibs.cmake)
#########################################
#########################################

if(USE_GPU_MOD)
	# GPU2
	find_package(GPU_MOD) #REQUIRED)
	if(NOT GPU_MOD_INCLUDE_DIR)
		set(GPU_MOD_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/gpu_mod/src) # fallback to submodule directory in case find_package didn't work
	endif()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GPU_MOD_CFLAGS} -I${GPU_MOD_INCLUDE_DIR} -DUSE_GPU_MOD")
	# TODO: set a function or a foreach to factorize the code configuration

	set(FAUST_SCALAR_FOR_GM float)
	set(GM_SCALAR float)
	set(GM_REINTERPRET_CAST_SCALAR float)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu_float.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu_float.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu_float.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu_float.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu_float.cpp)

	set(FAUST_SCALAR_FOR_GM double)
	set(GM_SCALAR double)
	set(GM_REINTERPRET_CAST_SCALAR double)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu_double.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu_double.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu_double.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu_double.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu_double.cpp)
	#target_sources(${FAUST_TARGET} PUBLIC ${GPU_MOD_CPP_FILES})
	#list(APPEND CPP_FILES )
	set(FAUST_SCALAR_FOR_GM std::complex<double>)
	set(GM_SCALAR cuDoubleComplex)
	set(GM_REINTERPRET_CAST_SCALAR double)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatDense_gpu_complexdouble.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatSparse_gpu_complexdouble.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_MatBSR_gpu_complexdouble.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Transform_gpu_complexdouble.cpp)
	configure_file(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu.cpp.in ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}/faust_Vect_gpu_complexdouble.cpp)
	# set install relative path of the lib (for each OS)
	if(APPLE)
		set(GPU_MOD_LIB libgm.dylib)
	elseif(UNIX)
		set(GPU_MOD_LIB libgm.so)
	elseif(WIN32)
		set(GPU_MOD_LIB gm.dll)
	endif()
	file(TO_NATIVE_PATH lib/${GPU_MOD_LIB} GPU_MOD_LIB_PATH) # used in setup.py
	file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${GPU_MOD_LIB} GPU_MOD_LIB_INSTALL_PATH) # only for system packages (not for pip package)
endif()


include_directories(BEFORE ${EIGEN_INC_DIR} ${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}  ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${FAUST_UTILS_SRC_DIR} ${FAUST_UTILS_SRC_DIR}/matio ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR})

if(FAUST_TORCH)
	find_package(Torch REQUIRED)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS} -DFAUST_TORCH") #-D_GLIBCXX_USE_CXX11_ABI=0")
	message(STATUS TORCH_LIBRARY=${TORCH_LIBRARY})
endif()


########### add_subdirectory ###############
if (BUILD_WRAPPER_PYTHON)
	add_subdirectory(${FAUST_PYTHON_SRC_DIR} ${FAUST_PYTHON_BIN_DIR})
endif(BUILD_WRAPPER_PYTHON)
if(MATIO_INC_DIR)
	include_directories(${MATIO_INC_DIR})
endif()
add_subdirectory(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_BIN_LINEAR_OPERATOR_DIR})
add_subdirectory(${FAUST_ALGORITHM_SRC_DIR} ${FAUST_ALGORITHM_BIN_DIR})
if(BUILD_TESTING)
	add_subdirectory(${FAUST_SRC_TEST_DIR} ${FAUST_BIN_TEST_BIN_DIR})
endif()
if (BUILD_DOCUMENTATION)
	add_subdirectory(${FAUST_DOC_SRC_DIR} ${FAUST_DOC_BIN_DIR})
endif ()


if (BUILD_WRAPPER_MATLAB)
	add_subdirectory(${FAUST_MATLAB_SRC_DIR} ${FAUST_MATLAB_MEX_BIN_DIR})
endif()
##################################################################




##### Add subdirectory ALGORITHM CPU et GPU version dans le main CmakeLists pour add library avec les cpp et pas que des lib objet 
#####    --> bug avec ide xcode ################################################################
#set(SRC_FILES "")
#file(GLOB CPP_FILES "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*.cpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp" "${FAUST_UTILS_SRC_DIR}/*.cpp")
file(GLOB_RECURSE CPP_FILES ${FAUST_SRC_DIR}/*.cpp)
file(GLOB_RECURSE GPU_CPP_FILES ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU/*.cpp)
file(GLOB_RECURSE GPU_MOD_CPP_FILES ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU2/*.cpp)
file(GLOB_RECURSE EIGTJ_SVDTJ_CPP_FILES ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*Givens* ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*SVDTJ*)
if(NOT USE_GPU_MOD)
	foreach(GPU_CPP_FILE IN LISTS GPU_MOD_CPP_FILES)
		list(REMOVE_ITEM CPP_FILES ${GPU_CPP_FILE})
	endforeach()
endif()
if(NOT BUILD_EIGTJ_SVDTJ)
	foreach(FILE IN LISTS EIGTJ_SVDTJ_CPP_FILES)
		message(STATUS "remove ${FILE} from build list of cpp files")
		list(REMOVE_ITEM CPP_FILES ${FILE})
	endforeach()
endif()
if(NOT FAUST_TORCH)
	list(REMOVE_ITEM CPP_FILES ${FAUST_SRC_LINEAR_OPERATOR_DIR}/CPU/faust_torch.cpp)
endif()
message(STATUS "FAµST cpp files: ${CPP_FILES}")
### include directory
#include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${EIGEN_INC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${FAUST_UTILS_SRC_DIR})
###############################################################################################################


if (BUILD_WRAPPER_MATLAB)
	add_library(${FAUST_TARGET} STATIC ${CPP_FILES} $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
else (BUILD_WRAPPER_MATLAB)
	add_library(${FAUST_TARGET} STATIC ${CPP_FILES} $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}>)
endif(BUILD_WRAPPER_MATLAB)

if(FAUST_TORCH)
	target_compile_definitions(${FAUST_TARGET} PUBLIC FAUST_TORCH -D_GLIBCXX_USE_CXX11_ABI=1)
	#target_compile_options(${FAUST_TARGET} PUBLIC ${TORCH_CXX_FLAGS})
	#message(STATUS "TORCH_CXX_FLAGS=" ${TORCH_CXX_FLAGS)
	get_filename_component(Torch_DIR ${TORCH_LIBRARY} DIRECTORY)
	message(STATUS TORCH_DIR=${Torch_DIR})
	include_directories(${FAUST_TARGET} ${Torch_DIR}/../include ${Torch_DIR}/../include/torch/csrc/api/include)
	#	target_link_libraries(${FAUST_TARGET} torch c10)
	target_include_directories(${FAUST_TARGET} PUBLIC  ${Torch_DIR}/../include ${Torch_DIR}/../include/torch/csrc/api/include)
	target_include_directories(${FAUST_LINEAR_OPERATOR_TARGET} PUBLIC  ${Torch_DIR}/../include ${Torch_DIR}/../include/torch/csrc/api/include)
endif()

if(BUILD_MULTITHREAD)
	target_compile_definitions(${FAUST_TARGET} PUBLIC _MUL_OMP_) #cf. ./src/faust_linear_operator/CPU/faust_prod_opt.hpp
endif()

if(USE_GPU_MOD)
	# TODO: test existence of build-cu9.2 # TODO: it must be a cmake variable
	# the headers for gpu_mod are the same for any version of cuda
	include_directories(${FAUST_TARGET} ${GPU_MOD_INCLUDE_DIR} ${GPU_MOD_INCLUDE_DIR}/../build-cu11.4)
	include_directories(${FAUST_TARGET} ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR})
	target_include_directories(${FAUST_TARGET} PUBLIC  ${GPU_MOD_INCLUDE_DIR} ${GPU_MOD_INCLUDE_DIR}/../build-cu11.4)
	target_link_libraries(${FAUST_TARGET} ${GPU_MOD_LIBS})
	if(UNIX)
		target_link_libraries(${FAUST_TARGET} dl)
	endif(UNIX)
	# TODO: used versions of CUDA and include dir should be CMake variables
	if(WIN32)
		# do not support CUDA 9.2 on Windows (require too old version VS <= 2017)
		set(CUDA_VERSIONS cu11.4)
	else()
		set(CUDA_VERSIONS cu9.2;cu11.4)
	endif()
	foreach(CUDA_CONF IN LISTS CUDA_VERSIONS)
		string(REGEX REPLACE "\\." "-${CUDA_CONF}." GPU_MOD_LIB_CUDA_CONF ${GPU_MOD_LIB})
		# if ${GPU_MOD_LIB} is already a symlink, then we've already gone through it, skip
		if(NOT IS_SYMLINK ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB})
			file(COPY ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB} DESTINATION ${PROJECT_BINARY_DIR})
			file(RENAME ${PROJECT_BINARY_DIR}/${GPU_MOD_LIB} ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB_CUDA_CONF})
			if (CUDA_CONF MATCHES "cu11.4")
				file(CREATE_LINK ${GPU_MOD_LIB_CUDA_CONF} ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB} RESULT SYMLINK_OK SYMBOLIC)
				install(FILES  ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_EXECUTE WORLD_READ)
			endif()
		endif()
		install(FILES  ${GPU_MOD_INCLUDE_DIR}/../build-${CUDA_CONF}/${GPU_MOD_LIB_CUDA_CONF} DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/ PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_EXECUTE WORLD_READ)
	endforeach()
	include_directories(${FAUST_TARGET} ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR})
endif()

target_link_libraries(${FAUST_TARGET})

if(CMAKE_BUILD_TYPE MATCHES "Release")
	install(CODE "file(GLOB TO_STRIP_FILES lib/*.a wrapper/python/*.so wrapper/matlab/mex/*.mexa64 wrapper/matlab/mex/*.mexa32 wrapper/matlab/mex/*.maci wrapper/matlab/mex/*.maci64)
	message(STATUS \"TO_STRIP_FILES=\${TO_STRIP_FILES}\")
	execute_process(COMMAND ${CMAKE_STRIP} \${TO_STRIP_FILES})
	")
endif(CMAKE_BUILD_TYPE MATCHES "Release")

if(NOT ${EXCLUDE_FAUST_LIB_INSTALL})
	install(TARGETS ${FAUST_TARGET} DESTINATION ${FAUST_INSTALL_LIB})
endif()

if(${BUILD_DOCUMENTATION})
	if(${CMAKE_INSTALL_PREFIX} MATCHES "faust")
		set(DOC_INSTALL_PREFIX  "${CMAKE_INSTALL_PREFIX}/doc")
	else()
		# CMAKE_INSTALL_PREFIX is not set
		string(REGEX REPLACE "(.*)/lib/?.*" "\\1/share/faust/doc" DOC_INSTALL_PREFIX ${FAUST_INSTALL_LIB})
	endif()
	install(FILES ${PROJECT_SOURCE_DIR}/license.txt ${PROJECT_SOURCE_DIR}/license_matio.txt DESTINATION ${DOC_INSTALL_PREFIX})
	install(DIRECTORY "${PROJECT_BINARY_DIR}/doc/html" DESTINATION "${DOC_INSTALL_PREFIX}")
endif()


##################################################################



##############################################################################
# Packaging using CPack:
##############################################################################
# USAGE 
# Use "package" target for a binary distribution. For example "make package" or "mingw32-make package"
# Use "package_source" target for a source distribution. For example "make package_source" or "mingw32-make package_source"
#
# For Unix Makefiles / MinGW Makefiles: Generates standard UNIX makefiles:
# A hierarchy of UNIX makefiles is generated into the build tree.
# Any standard UNIX-style make program can build the project through the default make target.
# A "install" target is provided. For example "make install" or "mingw32-make install"

# For Visual Studio and XCode project files
# A hierarchy of target is generated into the build tree and displayed in the IDE 
# Project provides "INSTALL" and "PACKAGE" target 

#In order to make a release:
#Set the release number in the file "version/buildversion.cmake"
#Make a tag with the number of the release
#Use Unix make to make a source release with command "make package_source"
#Upload the archive on the INRIA GForge FAUST project using the files admin page
#Add release in the corresponding categories FAUST with Add files page
#Filling the various required informations and validate the release
#Download the release and test if build correctly
#Make a binary relesase using Visual Studio using target "PACKAGE"
#Upload the archive on the INRIA GForge FAUST project using the files admin page
#Add release in the corresponding categories FAUST_win with Add files page
#Filling the various required informations and validate the release
#Download the release and test if build correctly

##############################################################################
# CPACK RULES
##############################################################################
#INCLUDE(InstallRequiredSystemLibraries)
if(WIN32)
	# TODO: move in NSIS specific CMakeLists.txt file
	if(BUILD_WRAPPER_PYTHON)
		if(DEFINED ENV{WIN_PY_VER})
			set(WIN_PY_VER $ENV{WIN_PY_VER})
		else()
			message(FATAL_ERROR "The environment variable WIN_PY_VER must be defined to indicate against which version of python the shared library (of pyfaust) must be compiled on Windows , e.g. 3.9)")
		endif()
	endif()
	configure_file(${PROJECT_SOURCE_DIR}/misc/nsis/faust.nsi ${PROJECT_BINARY_DIR}/faust.nsi @ONLY)
endif()
SET(CPACK_SET_DESTDIR "ON")
SET(CPACK_PACKAGE_RELOCATABLE "FALSE")
SET(CPACK_CMAKE_GENERATOR ${CMAKE_GENERATOR})
SET(CPACK_PACKAGE_NAME "FAUST")
SET(CPACK_PACKAGE_VENDOR "PANAMA Project-team INRIA")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FAUST_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FAUST_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FAUST_VERSION_PATCH}")
if(NOT DEFINED CPACK_PACKAGE_VERSION)
	SET(CPACK_PACKAGE_VERSION "${FAUST_VERSION_MAJOR}.${FAUST_VERSION_MINOR}.${FAUST_VERSION_PATCH}")
endif()
SET(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "FAUST ${FAUST_VERSION_MAJOR}.${FAUST_VERSION_MINOR}")

SET(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/gen_doc/images/logo.png")

SET(CPACK_OUTPUT_CONFIG_FILE "${PROJECT_BINARY_DIR}/CPackConfig.cmake")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flexible Approximate Multi-Layer Sparse Transform")
SET(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/license.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_URL "http://faust.inria.fr/" )
#INSTALL(FILES "${PROJECT_SOURCE_DIR}/README.txt" DESTINATION ${PROJECT_BINARY_DIR})
#INSTALL(FILES "${PROJECT_SOURCE_DIR}/license.txt" DESTINATION ${PROJECT_BINARY_DIR})


##############################################################################
# SOURCE DISTRIBUTION:
##############################################################################
set(CPACK_SOURCE_PACKAGE_FILE_NAME "FAUST-Source-${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.bzr/;/CMakeFiles/;/_CPack_Packages/;/.svn/;/.git/;/headerCPP.txt;/headerMatlab.txt;/headerCMakeList.txt;/gen_doc/LaTeX/images/;/gen_doc/installation_README.odt;/gen_doc/installation_Virtual_Machine.txt;/misc/data/txt/;/misc/data/mat/visage/;/misc/CTest_nightly/;/CTestConfig.cmake;/gen_doc/LaTeX/;/README.developer;")

set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/externals/unix/tarLibs/libxml2-2.9.4.tar.gz;/externals/unix/tarLibs/matio-1.5.7.7z;")
#if(WIN32 AND NOT UNIX)
#	set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/externals/unix/;/externals/win/eigen*/;/externals/win/OpenBLAS-v0.2.14-Win64-int64/;/externals/win/eigen-eigen-07105f7124f9/;")
#else()
#	set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/externals/win/;/externals/unix/eigen/;/externals/unix/OpenBLAS/;/externals/unix/sdk_OpenBLAS/;/externals/unix/matio/;")
#endif()
SET(CPACK_SOURCE_GENERATOR "TGZ")

##############################################################################
# BINARY DISTRIBUTION:
##############################################################################
SET(CPACK_PACKAGE_CONTACT "remi.gribonval@inria.fr")
SET(CPACK_INSTALL_CMAKE_PROJECTS "${PROJECT_BINARY_DIR};FAUST;ALL;/")

#SET(CPACK_IGNORE_FILES "/CMake/;/build/;/CMakeFiles/;/_CPack_Packages/;/src/;/bin/make_regression_constants;/install_manifest_/;.*~;/www/;/CVS//;/.svn/;/.git/;.cdtprojects;.project;/.settings/")
SET(CPACK_STRIP_FILES "")
#SET(CPACK_SOURCE_IGNORE_FILES "/build/;/CMakeFiles/;/_CPack_Packages/;/bin/make_regression_constants;/install_manifest_/;.*~;/www/;/CVS//;/.svn/;/.git/;.cdtprojects;.project;/.settings/")

IF(WIN32 AND NOT UNIX)
	# NSIS is used independently (see misc/nsis/faust.nsi)
	SET(CPACK_GENERATOR "NSIS")
	IF(CMAKE_CL_64)
		SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Windows")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Windows")
	ELSE(CMAKE_CL_64)
		SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Windows")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Windows")
  	ENDIF(CMAKE_CL_64)
	# NSIS not used yet, cause of a bug with suitable ZIP programme
	# There is a bug in NSI that does not handle full unix paths properly. Make
	# sure there is at least one set of four (4) backslashes.
	# UPDATE: NSIS is used externally, see misc/nsis/faust.nis
	# SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
	SET(CPACK_NSIS_DISPLAY_NAME "FAUST")
#	SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.irisa.fr/metiss/home_html-fr")
#	SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\http://www.irisa.fr/metiss/home_html-fr")
#	SET(CPACK_NSIS_CONTACT "xxxxxxxxxx@irisa.fr")
	SET(CPACK_NSIS_MODIFY_PATH OFF)
ELSE(WIN32 AND NOT UNIX)
	IF(APPLE)
		SET(CPACK_GENERATOR "PackageMaker")
		IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
			SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Mac")
			SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Mac")
		ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
			SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Mac")
			SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Mac")
		ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
		file(REMOVE ${PROJECT_BINARY_DIR}/postinstall) # needed to remove possible postinstall script left here 
							       # really needed only when BUILD_WRAPPER_MATLAB==OFF and BUILD_WRAPPER_PYTHON==ON
		if(BUILD_WRAPPER_MATLAB)
		    # prepare script to set matlab path for faust automatically after package installation
		    configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/macosx/postinst_matlab.sh.in ${PROJECT_BINARY_DIR}/postinstall @ONLY)
		endif()
		if(BUILD_WRAPPER_PYTHON)
		    # prepare the script allowing the user to avoid PYTHONPATH setting
		    configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/macosx/postinst_py.sh ${PROJECT_BINARY_DIR}/postinstall2 @ONLY)
		    file(READ ${PROJECT_BINARY_DIR}/postinstall2 POSTINST2_LINES)
		    file(APPEND ${PROJECT_BINARY_DIR}/postinstall "${POSTINST2_LINES}")
		    file(REMOVE ${PROJECT_BINARY_DIR}/postinstall2)
		endif()
		# post-install script to auto open readme
		configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/macosx/postinst_open_readme.sh ${PROJECT_BINARY_DIR}/postinstall2 @ONLY)
		file(READ ${PROJECT_BINARY_DIR}/postinstall2 POSTINST2_LINES)
		file(APPEND ${PROJECT_BINARY_DIR}/postinstall "${POSTINST2_LINES}")
		file(REMOVE ${PROJECT_BINARY_DIR}/postinstall2)
		# post-install script to download FAµST data from remote site into matlab wrapper data folder
		if(REMOTE_DATA_URL AND REMOTE_DATA_FILE)
			file(APPEND ${PROJECT_BINARY_DIR}/postinstall "mkdir -p @CMAKE_INSTALL_MATLAB_PREFIX@/data/;python  @CMAKE_INSTALL_PYTHON_PREFIX@/pyfaust/datadl.py @CMAKE_INSTALL_MATLAB_PREFIX@/data/ 2>&1 | tee -a /tmp/log_faust_install;")
		endif()
	ELSE(APPLE)
		SET(CPACK_GENERATOR "DEB;RPM")
		#SET(CPACK_DEBIAN_PACKAGE_DEPENDS "python, python-numpy, libfftw3-3, libsndfile1")
		#SET(CPACK_DEBIAN_PACKAGE_SECTION "math")
		SET(PACKAGE_DESCRIPTION "FAuST is a toolbox for decomposing a given dense matrix into a product of sparse matrices in order to reduce its computational complexity (both for storage and manipulation).")
		SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
		SET(CPACK_RPM_PACKAGE_DESCRIPTION "${PACKAGE_DESCRIPTION}")
		file(REMOVE ${PROJECT_BINARY_DIR}/postinst) # needed to remove possible postinstall script left here 
		# really needed only when BUILD_WRAPPER_MATLAB==OFF and BUILD_WRAPPER_PYTHON==ON
		if(BUILD_WRAPPER_MATLAB)
			# prepare script to set matlab path for faust automatically after package installation
			configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/debrpm/postinst_matlab.sh.in ${PROJECT_BINARY_DIR}/postinst @ONLY)
		endif()
		if(BUILD_WRAPPER_PYTHON)
			# prepare the script allowing the user to avoid PYTHONPATH setting
			configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/debrpm/postinst_py.sh ${PROJECT_BINARY_DIR}/postinst2 @ONLY)
			file(READ ${PROJECT_BINARY_DIR}/postinst2 POSTINST2_LINES)
			file(APPEND ${PROJECT_BINARY_DIR}/postinst "${POSTINST2_LINES}")
			file(REMOVE ${PROJECT_BINARY_DIR}/postinst2)
			#set(CPACK_RPM_PACKAGE_REQUIRES "${CPACK_RPM_PACKAGE_REQUIRES}, python3 >= 3.6.5, python3-numpy >= 1.13.3, python3-scipy >= 1.0.0")
			#set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, python3, python3-numpy, python3-scipy)
		endif()
		# postscript to auto open/indicate readme
		configure_file(${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/debrpm/postinst_open_readme.sh ${PROJECT_BINARY_DIR}/postinstall2 @ONLY)
		file(READ ${PROJECT_BINARY_DIR}/postinstall2 POSTINST2_LINES)
		file(APPEND ${PROJECT_BINARY_DIR}/postinst "${POSTINST2_LINES}")
		file(REMOVE ${PROJECT_BINARY_DIR}/postinstall2)

		#SET(CPACK_RPM_PACKAGE_DEPENDS "libmat, libmex, libmx ")
		#CPACK_RPM_PACKAGE_GROUP
		#SET (CPACK_RPM_PACKAGE_REQUIRES )
		#set(CPACK_RPM_COMPONENT_INSTALL ON)

		set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST "/opt") # to not conflict with filesystem package on fedora

		if(NOT ${USE_MATIO_STATIC_LIBS})
			set(CPACK_RPM_PACKAGE_REQUIRES "matio >= 1.5.7")
			set(CPACK_DEBIAN_PACKAGE_DEPENDS "libmatio-dev (>= 1.5.7)")
			# add postinst script to symlink matio lib if needed
			file(READ ${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/debrpm/postinst_matio_symlink.sh POSTINST2_LINES)
			file(APPEND ${PROJECT_BINARY_DIR}/postinst "${POSTINST2_LINES}")
			file(REMOVE ${PROJECT_BINARY_DIR}/postinst2)
		endif()
		# post-install script to download FAµST data from remote site into matlab wrapper data folder
		if(REMOTE_DATA_URL AND REMOTE_DATA_FILE)
			file(APPEND ${PROJECT_BINARY_DIR}/postinst "mkdir -p @CMAKE_INSTALL_MATLAB_PREFIX@/data/;python @CMAKE_INSTALL_PYTHON_PREFIX@/pyfaust/datadl.py  @CMAKE_INSTALL_MATLAB_PREFIX@/data/;")
		endif()
		#no need to exclude library auto-deps for debian package/ubuntu
		#it's only matters for rpm package (for example gpu_mod lib is linked to cusparse and cublas but we don't want to block install with auto-deps to libcusparse/cublas)
		# so let's ignore these library as deps
		set(CPACK_RPM_SPEC_MORE_DEFINE "%define __requires_exclude ^libmex[.]so.*|libmat[.]so.*|libmx[.]so.*|libpython.*[.]so.*|libcublas[.]so[.].*|libcusparse[.]so[.].*$ ")
		set(CPACK_RPM_PACKAGE_REQUIRES "python >= 3")
		set(CPACK_RPM_PACKAGE_RELEASE "1")
		SET(CPACK_RPM_PACKAGE_LICENSE "INRIA License (faust.inria.fr)")
		set(CPACK_RPM_PACKAGE_GROUP "unknown")
		set(CPACK_RPM_PACKAGE_VENDOR ${CPACK_PACKAGE_VENDOR})
		set(CPACK_RPM_PACKAGE_URL ${CPACK_PACKAGE_URL} )
		#set(CPACK_RPM_PACKAGE_DESCRIPTION ${CPACK_PACKAGE_DESCRIPTION_SUMMARY})
		IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
			if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
				SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Linux")
			else()
				SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-i386")
			endif()
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
			SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
		ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
			if(NOT DEFINED CPACK_PACKAGE_FILE_NAME)
				SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Linux")
			else()
				SET(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}-x86_64")
			endif()
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
			SET(CPACK_RPM_PACKAGE_ARCHITECTURE "x86_64")
		ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
		SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${CPACK_PACKAGE_FILE_NAME}")
		# Set post installation scripts for RPM and DEB packages
		SET(CPACK_RPM_POST_INSTALL_SCRIPT_FILE "${PROJECT_BINARY_DIR}/postinst")
		SET(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${PROJECT_BINARY_DIR}/postinst")
	ENDIF(APPLE)
ENDIF(WIN32 AND NOT UNIX)
#include(CPackConfig.cmake)
INCLUDE(CPack)
##############################################################################
##############################################################################





#add_dependencies(${FAUST_FAUSTCORE_TARGET} ${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})
#add_dependencies(${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET})




#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib" CACHE STRING "" FORCE)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin" CACHE STRING "" FORCE)
##################################################################


