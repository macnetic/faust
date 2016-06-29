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
	message(STATUS "********** install path : ${CMAKE_INSTALL_PREFIX} ************")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

###### Different option of configuration ######
option(BUILD_TESTING "Build Tests" ON) 
set(FAUST_USE_SINGLEPRECISION OFF CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_OPENBLAS OFF CACHE BOOL "Using openBLAS for matrix and vector computations")
set(FAUST_USE_MEX ON CACHE BOOL "Generate Mexfiles")
set(FAUST_USE_PROFILING OFF CACHE BOOL "Profiling the code")
set(FAUST_USE_OPENMP OFF CACHE BOOL "OpenMP Multithreading")
set(FAUST_USE_MATIO OFF CACHE BOOL "Using matio library to read mat files")
set(FAUST_USE_XML OFF CACHE BOOL "Using xml configuration to read xml files")
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
if(FAUST_USE_OPENMP)
	message(STATUS "**********OPENMP ACTIF*********")
endif()
if (FAUST_GEN_DOC)
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
if (FAUST_USE_MEX)
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
	###### BIN DIRECTORY ######
	set(FAUST_DEMO_BIN_DIR 		${FAUST_MATLAB_BIN_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_BIN_DIR 		${FAUST_DEMO_BIN_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_BIN_DIR 	${FAUST_DEMO_BIN_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_BIN_DIR ${FAUST_DEMO_BIN_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_BIN_DIR 	${FAUST_DEMO_BSL_BIN_DIR}/data/ CACHE INTERNAL "")
	set(FAUST_BIN_DEMO_BSL_OUTPUT_DIR 	${FAUST_DEMO_BSL_BIN_DIR}/output/ CACHE INTERNAL "")

	####  INSTALL DIRECTORY ####
	set(FAUST_DEMO_INSTALL_DIR ${FAUST_MATLAB_INSTALL_DIR}/demo CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_INSTALL_DIR 		${FAUST_DEMO_INSTALL_DIR}/Brain_source_localization CACHE INTERNAL "")
	set(FAUST_DEMO_TOOLS_INSTALL_DIR 	${FAUST_DEMO_INSTALL_DIR}/tools CACHE INTERNAL "")
	set(FAUST_DEMO_HADAMARD_INSTALL_DIR ${FAUST_DEMO_INSTALL_DIR}/Hadamard_factorization CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/data CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_DATA_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/data CACHE INTERNAL "")
	set(FAUST_DEMO_BSL_OUTPUT_INSTALL_DIR 	${FAUST_DEMO_BSL_INSTALL_DIR}/output CACHE INTERNAL "")

endif(FAUST_USE_MEX)


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
if (FAUST_USE_GPU)
	set(FAUST_SRC_LINEAR_OPERATOR_GPU_DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")
	set(FAUST_BIN_LINEAR_OPERATOR_GPU_DIR ${FAUST_BIN_LINEAR_OPERATOR_DIR}/GPU CACHE INTERNAL "")

	set(FAUST_SRC_TEST_GPU_DIR ${FAUST_MISC_DIR}/test/GPU_test CACHE INTERNAL "")
	set(FAUST_BIN_TEST_GPU_DIR ${FAUST_BIN_TEST_DIR}/binGPU CACHE INTERNAL "") 
endif(FAUST_USE_GPU)


###### DOCUMENTATION DIRECTORIES ######
if (FAUST_GEN_DOC)
	set(FAUST_DOC_SRC_DIR ${PROJECT_SOURCE_DIR}/gen_doc CACHE INTERNAL "")
	set(FAUST_DOC_BIN_DIR ${PROJECT_BINARY_DIR}/doc CACHE INTERNAL "")
endif()
##################################################################

#message(STATUS "FAUST_BIN_TEST_OUTPUT_DIR : ${FAUST_BIN_TEST_OUTPUT_DIR}")
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_OUTPUT_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_FIG_DIR})
file(MAKE_DIRECTORY ${FAUST_BIN_TEST_TMP_DIR})

if (FAUST_USE_MEX)
	file(MAKE_DIRECTORY ${FAUST_BIN_DEMO_BSL_OUTPUT_DIR})
	install(DIRECTORY ${FAUST_BIN_DEMO_BSL_OUTPUT_DIR} DESTINATION ${FAUST_DEMO_BSL_OUTPUT_INSTALL_DIR})
endif (FAUST_USE_MEX)

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
if (FAUST_USE_GPU)
	set(FAUST_LINEAR_OPERATOR_GPU_TARGET faust_linear_operator_gpu CACHE INTERNAL "")
	#set(FAUST_MATRIX_CU_TARGET faust_gpu_matrix CACHE INTERNAL "")
	#set(FAUST_FAUSTCORE_CU_TARGET faust_gpu_core CACHE INTERNAL "")
	set(FAUST_ALGORITHM_CU_TARGET faust_algorithm_gpu CACHE INTERNAL "")
endif(FAUST_USE_GPU)
##################################################################


###### DEBUG vs RELEASE TARGETS ######
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
		if (FAUST_USE_OPENMP)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp " CACHE STRING "compile flags" FORCE)
		endif(FAUST_USE_OPENMP)
					
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fexceptions -fno-omit-frame-pointer -pthread -fPIC -O2 " CACHE STRING "compile flags" FORCE)
		
		
	elseif(WIN32)
		if (FAUST_USE_OPENMP)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp " CACHE STRING "compile flags" FORCE)
		endif(FAUST_USE_OPENMP)
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



###### find external libraries ######
if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif(WIN32)



#find_package(eigen)


# looking for library and include files in the environment variables
include(CMake/add_library_path.cmake)
include(CMake/check_external_libraries.cmake)

# /sw /opt/local
add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/local/lib" "/opt" "/opt/local" "/usr/lib/x86_64-linux-gnu/")
add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" "${PROJECT_SOURCE_DIR}/externals/eigen")


#add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/local/lib" "/opt" "/opt/local" "/usr/lib/x86_64-linux-gnu/" "${PROJECT_SOURCE_DIR}/externals" )
#add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "$ENV{EIGENDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" "${PROJECT_SOURCE_DIR}/externals" )

set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")
message(STATUS "LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")


check_external_libraries(matio MATIO_LIB_FILE 0)
check_external_libraries(xml2 XML2_LIB_FILE 0)
check_external_libraries(hdf5 HDF5_LIB_FILE 0)


#LDFLAGS = -L$(CUDADIR)/lib64 -L$(MATIODIR)/lib -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio -lhdf5
check_external_includes("libxml/parser.h" XML_INC_DIR 0)
check_external_includes("matio.h" MATIO_INC_DIR 0)

### find and/or install eigen lib (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)
check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)

if (EIGEN_INC_DIR)
	message(STATUS "Eigen library is available here : EIGEN_INC_DIR= ${EIGEN_INC_DIR}")
else (EIGEN_INC_DIR)
	if(UNIX)
		#exec_program("which hg | xargs echo" OUTPUT_VARIABLE MERCURIAL_HG_PATH_DIR)
	   	#exec_program("readlink ${MERCURIAL_HG_PATH_DIR}" OUTPUT_VARIABLE READLINK_TMP)
		if (MERCURIAL_HG_PATH_DIR)
			exec_program("hg clone https://bitbucket.org/eigen/eigen/ ${CMAKE_SOURCE_DIR}/externals/eigen")
		else (MERCURIAL_HG_PATH_DIR)
			exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/tmp https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2")
			exec_program("tar jxvf ${CMAKE_SOURCE_DIR}/externals/tmp/*.tar.bz2 -C ${CMAKE_SOURCE_DIR}/externals")
			exec_program("mv ${CMAKE_SOURCE_DIR}/externals/eigen-* ${CMAKE_SOURCE_DIR}/externals/eigen")
		endif()
	elseif(WIN32)
		exec_program("where hg.exe" MERCURIAL_HG_PATH_DIR)
		if (MERCURIAL_HG_PATH_DIR)
			exec_program("hg clone https://bitbucket.org/eigen/eigen/ ${CMAKE_SOURCE_DIR}/externals/eigen")
		else (MERCURIAL_HG_PATH_DIR)	
			exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/eigen https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2")
			exec_program("tar jxvf ${CMAKE_SOURCE_DIR}/externals/eigen/*.tar.bz2 ${CMAKE_SOURCE_DIR}/externals/eigen/")
		endif(MERCURIAL_HG_PATH_DIR)	
	
	else(UNIX)
		message(WARNING "Unknown type of plateform for library Eigen")	
	endif(UNIX)

	add_include_path(INCLUDE_PATH_LIST_TMP2 "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" "${PROJECT_SOURCE_DIR}/externals/eigen")

	#add_include_path(INCLUDE_PATH_LIST_TMP2 "${PROJECT_SOURCE_DIR}/externals/eigen")
	#set(INCLUDE_PATH_LIST2 ${INCLUDE_PATH_LIST_TMP2} CACHE PATH "List of include paths used as PATH parameter in find_path")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP2})
	message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")
	check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)
	message(STATUS "EIGEN_INC_DIR= ${EIGEN_INC_DIR}")
endif (EIGEN_INC_DIR)
#############################################################




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
if (FAUST_USE_MEX)
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
	# LINUX AND APPLE METHOD ARE VERY SIMILAR, CODE COULD BE FACTORIZED
	if(UNIX) 
		if(APPLE)
			exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i mac" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
                	if("${MEX_SUBDIR_LIB}" STREQUAL "maci64")
                        	set(MEX_EXT  "mexmaci64")
                	elseif("${MEX_SUBDIR_LIB}" STREQUAL "maci")
                        	set(MEX_EXT  "mexmaci")
                	endif()		
		# METHODE 1 (without using matlab)
		else(APPLE)		
			exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i glnx" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
			if("${MEX_SUBDIR_LIB}" STREQUAL "glnxa64")
				set(MEX_EXT  "mexa64")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "glnx86")
				set(MEX_EXT  "mexa32")
			endif()
		endif(APPLE)	
	

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
endif(FAUST_USE_MEX)
##################################################################



###### TARGET GENERATION ######
add_subdirectory(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_BIN_LINEAR_OPERATOR_DIR})
add_subdirectory(${FAUST_ALGORITHM_SRC_DIR} ${FAUST_ALGORITHM_BIN_DIR})
add_subdirectory(${FAUST_SRC_TEST_DIR} ${FAUST_BIN_TEST_BIN_DIR})
if (FAUST_GEN_DOC)
	add_subdirectory(${FAUST_DOC_SRC_DIR} ${FAUST_DOC_BIN_DIR})
endif ()
add_subdirectory(${FAUST_SRC_CMDLINE_DIR} ${FAUST_BIN_CMDLINE_BIN_DIR})
if (FAUST_USE_GPU)
	add_subdirectory(${FAUST_SRC_LINEAR_OPERATOR_GPU_DIR} ${FAUST_BIN_LINEAR_OPERATOR_GPU_DIR})
	#add_subdirectory(${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_FAUSTCORE_CU_BIN_DIR})
	#add_subdirectory(${FAUST_ALGORITHM_CU_SRC_DIR} ${FAUST_ALGORITHM_CU_BIN_DIR})
	add_subdirectory(${FAUST_SRC_TEST_GPU_DIR} ${FAUST_BIN_TEST_GPU_DIR})
endif (FAUST_USE_GPU)

if (FAUST_USE_MEX)
	add_subdirectory(${FAUST_MATLAB_SRC_DIR} ${FAUST_MATLAB_MEX_BIN_DIR})
	add_subdirectory(${FAUST_DEMO_SRC_DIR} ${FAUST_DEMO_BIN_DIR})
endif()
##################################################################


###### TARGET GENERATION ######
#$<TARGET_OBJECTS:${FAUST_FAUSTCORE_CU_TARGET}> 
if (FAUST_USE_GPU)
	if (FAUST_USE_MEX)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_CU_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (FAUST_USE_MEX)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_CU_TARGET}>)
	endif(FAUST_USE_MEX)
else (FAUST_USE_GPU)
	if (FAUST_USE_MEX)
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
	else (FAUST_USE_MEX)	
		add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_LINEAR_OPERATOR_TARGET}> $<TARGET_OBJECTS:${FAUST_ALGORITHM_TARGET}>)
	endif(FAUST_USE_MEX)
endif (FAUST_USE_GPU)


if (FAUST_USE_MEX)
	if (UNIX) 
		if(APPLE)
			target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}/libmex.dylib)
		else(APPLE)
			target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}/libmex.so)
		endif(APPLE)
		
	elseif(WIN32)
		message(FATAL_ERROR "MEXFILE : CANNOT SEARCH LIBMEX IN WINDOWS")
	endif()
endif(FAUST_USE_MEX)

# In case of GPU, the "CUDA" library libfaust_matrix_cu.a is linked with faust (libfaust_matrix_cu.a is not a OBJECT library) 
if (FAUST_USE_GPU)
	#target_link_libraries(${FAUST_TARGET} ${FAUST_MATRIX_CU_TARGET} )	
	target_link_libraries(${FAUST_TARGET} ${FAUST_LINEAR_OPERATOR_GPU_TARGET} )	
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



