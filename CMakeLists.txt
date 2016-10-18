##############################################################################
##                              Description:                                ##
##  Main cmakelist file for the faust project                               ##
##  Other important files for CMake are                                     ##
##      - CMakeLists.txt files in the sub-directories                       ##
##      - CMake/* files, which define some tools to search for libraries.   ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.gforge.inria.fr>                         ##
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

set(FAUST_VERSION_MAJOR 2)
set(FAUST_VERSION_MINOR 0)
set(FAUST_VERSION_PATCH 0)
#set(FAUST_VERSION 2.0)

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
#		set(FAUST_INSTALL_CMDLINE 		"$ENV{HOME}/Documents/faust/CmdLine" CACHE STRING "default install path for the cmd faust" FORCE )
#		set(FAUST_INSTALL_TESTING_DATA 	"$ENV{HOME}/Documents/faust/testing/data" CACHE STRING "default install path for the test faust" FORCE )
		set(FAUST_INSTALL_BIN 			"${PROJECT_BINARY_DIR}/output/bin" CACHE STRING "default install path for the binary faust" FORCE )
		set(FAUST_INSTALL_BIN_TIME_COMP "${PROJECT_BINARY_DIR}/output/bin/runtime_comparison" CACHE STRING "default install path for the binary faust" FORCE )
		set(FAUST_INSTALL_DOC 			"${PROJECT_BINARY_DIR}/output/doc" CACHE STRING "default install path for the documentation faust" FORCE )
		#set(FAUST_INSTALL_CMDLINE 		"${PROJECT_BINARY_DIR}/output/CmdLine" CACHE STRING "default install path for the cmd faust" FORCE )
		#set(FAUST_INSTALL_TESTING_DATA 	"${PROJECT_BINARY_DIR}/output/testing/data" CACHE STRING "default install path for the test faust" FORCE )
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "$ENV{HOME}/Documents/MATLAB/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		# CMAKE_INSTALL_PREFIX path is deleted for next cmake .. configure
		SET(CMAKE_INSTALL_PREFIX "" CACHE STRING "default install path" FORCE )
		#message (STATUS "default install") 
	# Install in the CMAKE_INSTALL_PREFIX path defined by the user (local install) 	
	elseif(CMAKE_INSTALL_PREFIX)
		SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
		message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
		
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}/lib")
		set(FAUST_INSTALL_BIN 			"${CMAKE_INSTALL_PREFIX}/bin")
		set(FAUST_INSTALL_BIN_TIME_COMP "${CMAKE_INSTALL_PREFIX}/bin/runtime_comparison")
		set(FAUST_INSTALL_DOC 			"${CMAKE_INSTALL_PREFIX}/doc")
		#set(FAUST_INSTALL_CMDLINE 		"${CMAKE_INSTALL_PREFIX}/CmdLine")
		#set(FAUST_INSTALL_TESTING_DATA 	"${CMAKE_INSTALL_PREFIX}/testing/data")
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "${CMAKE_INSTALL_PREFIX}/matlab")
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
		set(FAUST_INSTALL_BIN 			"${PROJECT_BINARY_DIR}/output/bin" CACHE STRING "default install path for the binary faust" FORCE )
		set(FAUST_INSTALL_BIN_TIME_COMP "${PROJECT_BINARY_DIR}/output/bin/runtime_comparison" CACHE STRING "default install path for the binary faust" FORCE )
		set(FAUST_INSTALL_DOC 			"${PROJECT_BINARY_DIR}/output/doc" CACHE STRING "default install path for the documentation faust" FORCE )
		#set(FAUST_INSTALL_CMDLINE 		"${PROJECT_BINARY_DIR}/output/CmdLine" CACHE STRING "default install path for the cmd faust" FORCE )
		#set(FAUST_INSTALL_TESTING_DATA 	"${PROJECT_BINARY_DIR}/output/testing/data" CACHE STRING "default install path for the test faust" FORCE )
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "C:/Users/$ENV{USERNAME}/Documents/MATLAB/faust" CACHE STRING "default install path for the matlab wrapper" FORCE )
		# CMAKE_INSTALL_PREFIX path is deleted for next cmake .. configure
		SET(CMAKE_INSTALL_PREFIX "" CACHE STRING "default install path" FORCE )
		#message (STATUS "default install") 
	# Install in the CMAKE_INSTALL_PREFIX path defined by the user (local install) 	
	elseif(CMAKE_INSTALL_PREFIX)
		SET(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" CACHE STRING "new install path" FORCE )
		message(STATUS "Install FAuST path : ${CMAKE_INSTALL_PREFIX}")
		
		set(FAUST_INSTALL_LIB 			"${CMAKE_INSTALL_PREFIX}/lib")
		set(FAUST_INSTALL_BIN 			"${CMAKE_INSTALL_PREFIX}/bin")
		set(FAUST_INSTALL_BIN_TIME_COMP "${CMAKE_INSTALL_PREFIX}/bin/runtime_comparison")
		set(FAUST_INSTALL_DOC 			"${CMAKE_INSTALL_PREFIX}/doc")
		#set(FAUST_INSTALL_CMDLINE 		"${CMAKE_INSTALL_PREFIX}/CmdLine")
		#set(FAUST_INSTALL_TESTING_DATA 	"${CMAKE_INSTALL_PREFIX}/testing/data")
		# default location of install dir for matlab component
		SET(CMAKE_INSTALL_MATLAB_PREFIX "${CMAKE_INSTALL_PREFIX}/matlab")
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
option(BUILD_USE_GPU  "Using both CPU and GPU process" OFF)

option(BUILD_MATLAB_MEX_FILES "Enable building Matlab MEX files." ON)
# Using openBLAS for matrix and vector computations
option(BUILD_OPENBLAS "Using openBLAS for matrix and vector computations" OFF)
# different type of files 
option(BUILD_READ_XML_FILE "Using xml configuration to read xml files" OFF)
option(BUILD_READ_MAT_FILE "Using matio library to read mat files" OFF)
#CMAKE_DEPENDENT_OPTION

###### Different option of configuration ######
#option(BUILD_TESTING "Build Tests" ON) 
set(BUILD_USE_SINGLEPRECISION OFF) # CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(BUILD_USE_PROFILING OFF) # CACHE BOOL "Profiling the code")

#set(BUILD_MULTITHREAD OFF CACHE BOOL "OpenMP Multithreading")
#set(BUILD_VERBOSE OFF CACHE BOOL "Display useful message for debugging")
#set(BUILD_USE_GPU OFF CACHE BOOL "Using both CPU and GPU process")
#set(BUILD_DOCUMENTATION OFF CACHE BOOL "Generate html documentation with doxygen")
#set(BUILD_READ_XML_FILE OFF CACHE BOOL "Using xml configuration to read xml files")
#set(BUILD_MATLAB_MEX_FILES ON CACHE BOOL "Generate Mexfiles")
#set(BUILD_OPENBLAS OFF CACHE BOOL "Using openBLAS for matrix and vector computations")
#set(BUILD_READ_MAT_FILE OFF CACHE BOOL "Using matio library to read mat files")

if (BUILD_USE_SINGLEPRECISION)
	message(STATUS "**********SINGLEPRECISION ACTIF************")
endif()
if (BUILD_OPENBLAS)
	message(STATUS "**********OPENBLAS ACTIF************")
endif()
if (BUILD_MATLAB_MEX_FILES)
	message(STATUS "**********MEX ACTIF************")
	message(STATUS "CMAKE_INSTALL_MATLAB_PREFIX = ${CMAKE_INSTALL_MATLAB_PREFIX}")
endif()
if (BUILD_USE_PROFILING)
	message(STATUS "**********PROFILING ACTIF************")
endif()
if (BUILD_READ_MAT_FILE)
	message(STATUS "**********MATIO ACTIF************")
	message(STATUS "FAUST_INSTALL_BIN = ${FAUST_INSTALL_BIN}")
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
set(FAUST_CONFIG_MAT_DIR 	${FAUST_MISC_DIR}/config_file/mat CACHE INTERNAL "")


###### MISC DIRECTORIES ######
set(FAUST_SRC_TEST_DIR 	${FAUST_MISC_DIR}/test CACHE INTERNAL "")
set(FAUST_DATA_MAT_DIR 	${FAUST_MISC_DIR}/data/mat CACHE INTERNAL "")

set(FAUST_SRC_TEST_TIME_COMP_DIR 	${FAUST_SRC_TEST_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_SRC_TEST_TOOL_DIR 		${FAUST_SRC_TEST_DIR}/tools CACHE INTERNAL "")
set(FAUST_SRC_TEST_SRC_DIR 			${FAUST_SRC_TEST_DIR}/src CACHE INTERNAL "")



###### WRAPPER MATLAB DIRECTORIES ######
if (BUILD_MATLAB_MEX_FILES)
	#version MATLAB du code (entierement MATLAB, pas de mexfunction)
	#file(MAKE_DIRECTORY ${CMAKE_INSTALL_MATLAB_PREFIX})
	###### SRC DIRECTORY #####
	set(FAUST_MATLAB_SRC_DIR ${FAUST_SRC_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/tools CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/src CACHE INTERNAL "")
	set(FAUST_MATLAB_DOC_SRC_DIR 	${FAUST_MATLAB_SRC_DIR}/doc CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_SRC_DIR ${FAUST_MATLAB_SRC_DIR}/old_matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_PROX_SRC_DIR ${FAUST_MATLAB_OLD_SRC_DIR}/proxs CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_UTILS_SRC_DIR ${FAUST_MATLAB_OLD_SRC_DIR}/utils CACHE INTERNAL "")
	###### BIN DIRECTORY ######
	set(FAUST_MATLAB_BIN_DIR ${FAUST_BIN_WRAPPER_DIR}/matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_MEX_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/mex CACHE INTERNAL "")
	set(FAUST_MATLAB_TOOLS_BIN_DIR 	${FAUST_MATLAB_BIN_DIR}/tools CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_BIN_DIR ${FAUST_MATLAB_BIN_DIR}/old_matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_PROX_BIN_DIR ${FAUST_MATLAB_OLD_BIN_DIR}/proxs CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_UTILS_BIN_DIR ${FAUST_MATLAB_OLD_BIN_DIR}/utils CACHE INTERNAL "")
	####  INSTALL DIRECTORY ####
	#set(FAUST_MATLAB_INSTALL_DIR ${CMAKE_INSTALL_MATLAB_PREFIX})
	
	#set(FAUST_MATLAB_MEX_INSTALL_DIR 	${CMAKE_INSTALL_MATLAB_PREFIX}/mex)
	#set(FAUST_MATLAB_TOOLS_INSTALL_DIR 	${CMAKE_INSTALL_MATLAB_PREFIX}/tools CACHE INTERNAL "")	


	########## DEMO DIRECTORIES #############
	###### SRC DIRECTORY #####
	set(FAUST_DEMO_SRC_DIR ${FAUST_MISC_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_SRC_DIR 		${FAUST_DEMO_SRC_DIR}/Brain_source_localization CACHE INTERNAL "")
	#set(FAUST_DEMO_BSL_DATA_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Brain_source_localization/data/ CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_SRC_DIR 	${FAUST_DEMO_SRC_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_TIMECOMP_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Runtime_comparison CACHE INTERNAL "")
	set(FAUST_DEMO_QUICKSTART_SRC_DIR ${FAUST_DEMO_SRC_DIR}/Quick_start CACHE INTERNAL "")

	###### BIN DIRECTORY ######
	set(FAUST_DEMO_BIN_DIR 		${FAUST_MATLAB_BIN_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_BIN_DIR 		${FAUST_DEMO_BIN_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_BIN_DIR 	${FAUST_DEMO_BIN_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_BIN_DIR ${FAUST_DEMO_BSL_BIN_DIR}/data/ CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_OUTPUT_BIN_DIR 	${FAUST_DEMO_BSL_BIN_DIR}/output/ CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_OUTPUT_BIN_DIR 	${FAUST_DEMO_HADAMARD_BIN_DIR}/output/ CACHE INTERNAL "")
	set(FAUST_DEMO_TIMECOMP_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Runtime_comparison CACHE INTERNAL "")
	set(FAUST_DEMO_FIG_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Figures CACHE INTERNAL "")
	set(FAUST_DEMO_QUICKSTART_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Quick_start CACHE INTERNAL "")
	

	####  INSTALL DIRECTORY ####
	# No need to specify installation demo  directory because 
	# it has the same architecture as FAUST_DEMO_BIN_DIR
	# so FAUST_DEMO_BIN_DIR is directly installed


endif(BUILD_MATLAB_MEX_FILES)
##################################################################

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

### In the case of GPU preocess 
if (BUILD_USE_GPU)
	set(FAUST_SRC_LINEAR_OPERATOR_GPU_DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")
	set(FAUST_BIN_LINEAR_OPERATOR_GPU_DIR ${FAUST_BIN_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")

	set(FAUST_SRC_TEST_GPU_DIR ${FAUST_MISC_DIR}/test/GPU_test CACHE INTERNAL "")
	set(FAUST_BIN_TEST_GPU_DIR ${FAUST_BIN_TEST_DIR}/binGPU CACHE INTERNAL "") 
endif(BUILD_USE_GPU)
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


if (BUILD_MATLAB_MEX_FILES)
	file(MAKE_DIRECTORY ${FAUST_DEMO_BSL_OUTPUT_BIN_DIR})
	file(MAKE_DIRECTORY ${FAUST_DEMO_HADAMARD_OUTPUT_BIN_DIR})
	file(MAKE_DIRECTORY ${FAUST_DEMO_FIG_BIN_DIR})
endif (BUILD_MATLAB_MEX_FILES)



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
#if( (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
if (BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
	if (UNIX)	
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} -O1")

		if(APPLE)
			#message(STATUS "APPLE OPERATING SYSTEM")
			if (BUILD_MULTITHREAD)
				message (WARNING "-fopen is not activated for MAC OS.")
			endif(BUILD_MULTITHREAD)
			# -pthread is not take into account with MAC platform. 
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -fPIC " CACHE STRING "compile flags" FORCE)
		else(APPLE)
			#message(STATUS "UNIX OPERATING SYSTEM") 
			if (BUILD_MULTITHREAD)
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
			endif(BUILD_MULTITHREAD)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -pthread -fPIC " CACHE STRING "compile flags" FORCE)
		endif(APPLE)
	elseif(WIN32)
		set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_DEBUG} /O0")
		if (BUILD_MULTITHREAD)
			if (${CMAKE_GENERATOR} MATCHES "MinGW Makefiles" )
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
			else()
				#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp " CACHE STRING "compile flags" FORCE)
				message (WARNING "openmp is not activated for visual studio compiler.")
			endif()		
		endif(BUILD_MULTITHREAD)
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")
	endif() 
	message(STATUS "**********DEBUG mode************")
else(BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
	if (UNIX)	

		if(APPLE)
			#message(STATUS "APPLE OPERATING SYSTEM")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -fPIC -O3" CACHE STRING "compile flags" FORCE)
			# -pthread is not take into account with MAC platform. 
			if (BUILD_MULTITHREAD)
				message (WARNING "-fopen is not activated for MAC OS.")
			endif(BUILD_MULTITHREAD)
		else(APPLE)
			#message(STATUS "UNIX OPERATING SYSTEM") 
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -pthread -fPIC -O3" CACHE STRING "compile flags" FORCE)
			if (BUILD_MULTITHREAD)
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
			endif(BUILD_MULTITHREAD)
		endif(APPLE)
					
	elseif(WIN32)
		if (BUILD_MULTITHREAD)
			#message (STATUS "CMAKE_GENERATOR=${CMAKE_GENERATOR}")
			if (${CMAKE_GENERATOR} MATCHES "MinGW Makefiles" )
				set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
			else()
				#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp " CACHE STRING "compile flags" FORCE)
				message (WARNING "openmp is not activated for visual studio compiler.")
			endif()	
		endif(BUILD_MULTITHREAD)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHs /MD /O2" CACHE STRING "compile flags" FORCE)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHs /MD /O2" CACHE STRING "compile flags" FORCE)
	else()
		message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
	endif()
	message(STATUS "**********RELEASE mode************")
endif(BUILD_DEBUG OR (CMAKE_BUILD_TYPE MATCHES "Debug") OR  (CMAKE_BUILD_TYPE MATCHES "debug") OR  (CMAKE_BUILD_TYPE MATCHES "DEBUG") )
############################################################################

#if(UNIX)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp -fexceptions -fno-omit-frame-pointer -pthread -fPIC " CACHE STRING "compile flags" FORCE)
#elseif(WIN32)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /EHs /openmp /MD" CACHE STRING "compile flags" FORCE)
#else()
#	message(WARNING "Unknown type of plateform for CMAKE_CXX_FLAGS")		
#endif(UNIX)
#message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
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
		target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}/libmex.dll)
		#message(FATAL_ERROR "MEXFILE : CANNOT SEARCH LIBMEX IN WINDOWS")
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
SET(CPACK_SET_DESTDIR "ON")
SET(CPACK_PACKAGE_RELOCATABLE "FALSE")
SET(CPACK_CMAKE_GENERATOR ${CMAKE_GENERATOR})
SET(CPACK_PACKAGE_NAME "FAUST")
SET(CPACK_PACKAGE_VENDOR "PANAMA Project-team INRIA")
SET(CPACK_PACKAGE_VERSION_MAJOR "${FAUST_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${FAUST_VERSION_MINOR}")
SET(CPACK_PACKAGE_VERSION_PATCH "${FAUST_VERSION_PATCH}")
SET(CPACK_PACKAGE_VERSION "${FAUST_VERSION_MAJOR}.${FAUST_VERSION_MINOR}.${FAUST_VERSION_PATCH}")
SET(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_PROCESSOR}")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "FAUST ${FAUST_VERSION_MAJOR}.${FAUST_VERSION_MINOR}")

SET(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/gen_doc/images/logo.png")

SET(CPACK_OUTPUT_CONFIG_FILE "${PROJECT_BINARY_DIR}/CPackConfig.cmake")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.txt")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flexible Approximate Multi-Layer Sparse Transform")
SET(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/license.txt")
SET(CPACK_RESOURCE_FILE_WELCOME "${PROJECT_SOURCE_DIR}/README.txt")
#INSTALL(FILES "${PROJECT_SOURCE_DIR}/README.txt" DESTINATION ${PROJECT_BINARY_DIR})
#INSTALL(FILES "${PROJECT_SOURCE_DIR}/license.txt" DESTINATION ${PROJECT_BINARY_DIR})


##############################################################################
# SOURCE DISTRIBUTION:
##############################################################################
set(CPACK_SOURCE_PACKAGE_FILE_NAME "FAUST-Source-${FAUST_VERSION_MAJOR}.${FAUST_VERSION_MINOR}.${FAUST_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.bzr/;/CMakeFiles/;/_CPack_Packages/;/.svn/;/.git/;/headerCPP.txt;/headerMatlab.txt;/headerCMakeList.txt;/gen_doc/LaTeX/images/;/gen_doc/installation_README.odt;/gen_doc/installation_Virtual_Machine.txt;/misc/data/txt/;/misc/data/mat/visage/")

#if(WIN32 AND NOT UNIX)
#	set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/externals/unix/;/externals/win/eigen*/;/externals/win/OpenBLAS-v0.2.14-Win64-int64/;/externals/win/eigen-eigen-07105f7124f9/;")
#else()
#	set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};/externals/win/;/externals/unix/eigen/;/externals/unix/OpenBLAS/;/externals/unix/sdk_OpenBLAS/;/externals/unix/matio/;")
#endif()


SET(CPACK_SOURCE_GENERATOR "TGZ")

##############################################################################
# BINARY DISTRIBUTION:
##############################################################################
SET(CPACK_PACKAGE_CONTACT "name@inria.fr")
SET(CPACK_INSTALL_CMAKE_PROJECTS "${PROJECT_BINARY_DIR};FAUST;ALL;/")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Flexible Approximate Multi-Layer Sparse Transform")
#SET(CPACK_IGNORE_FILES "/CMake/;/build/;/CMakeFiles/;/_CPack_Packages/;/src/;/bin/make_regression_constants;/install_manifest_/;.*~;/www/;/CVS//;/.svn/;/.git/;.cdtprojects;.project;/.settings/")
SET(CPACK_STRIP_FILES "")

#SET(CPACK_SOURCE_IGNORE_FILES "/build/;/CMakeFiles/;/_CPack_Packages/;/bin/make_regression_constants;/install_manifest_/;.*~;/www/;/CVS//;/.svn/;/.git/;.cdtprojects;.project;/.settings/")

IF(WIN32 AND NOT UNIX)
	SET(CPACK_GENERATOR "NSIS")
	IF(CMAKE_CL_64)
		SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Windows")
    	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Windows")
	ELSE(CMAKE_CL_64)
    	SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Windows")
    	SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Windows")
  	ENDIF(CMAKE_CL_64)
	SET(CPACK_SET_DESTDIR "OFF")
	# NSIS not used yet, cause of a bug with suitable ZIP programme
	# There is a bug in NSI that does not handle full unix paths properly. Make
	# sure there is at least one set of four (4) backlasshes.
	# SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/Utilities/Release\\\\InstallIcon.bmp")
	SET(CPACK_NSIS_DISPLAY_NAME "FAUST")
#	SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.irisa.fr/metiss/home_html-fr")
#	SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\http://www.irisa.fr/metiss/home_html-fr")
#	SET(CPACK_NSIS_CONTACT "matchingpursuite@irisa.fr")
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
	ELSE(APPLE)
		SET(CPACK_GENERATOR "DEB;RPM")
		#SET(CPACK_DEBIAN_PACKAGE_DEPENDS "python, python-numpy, libfftw3-3, libsndfile1")
		#SET(CPACK_DEBIAN_PACKAGE_SECTION "math")
		SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "FAuST is a C++ toolbox, useful to decompose a given dense matrix into a product of sparse matrices in order to reduce its computational complexity (both for storage and manipulation). ")
		IF(CMAKE_SIZEOF_VOID_P MATCHES "4")
			SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Linux")
			SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-i386-Linux")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")
		ELSE(CMAKE_SIZEOF_VOID_P MATCHES "4")
			SET(CPACK_PACKAGE_FILE_NAME "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Linux")
			SET(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "FAUST-binary-${CPACK_PACKAGE_VERSION}-x86_64-Linux")
			SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
		ENDIF(CMAKE_SIZEOF_VOID_P MATCHES "4")
	ENDIF(APPLE)
ENDIF(WIN32 AND NOT UNIX)
INCLUDE(CPack)
##############################################################################
##############################################################################


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



