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
set(FAUST_TEST_SRC_DIR ${PROJECT_SOURCE_DIR}/test CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_SRC_DIR ${FAUST_TEST_SRC_DIR}/runtime_comparison CACHE INTERNAL "")
set(FAUST_TEST_SRC_TOOL_DIR ${FAUST_TEST_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_TESTSRC_SRC_DIR ${FAUST_TEST_SRC_DIR}/src CACHE INTERNAL "")
set(FAUST_TESTDATA_SRC_DIR ${FAUST_TEST_SRC_DIR}/data CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_SRC_DIR ${PROJECT_SOURCE_DIR}/interface_matlab CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_SRC_DIR ${FAUST_MEXINTERFACE_SRC_DIR}/mex_functions CACHE INTERNAL "")
set(FAUST_EIGEN_SRC_DIR ${PROJECT_SOURCE_DIR}/eigen CACHE INTERNAL "")

# BINARY DIRECTORIES
set(FAUST_WRAPPER_BIN_DIR ${PROJECT_BINARY_DIR}/wrapper)
set(FAUST_EXCEPTION_BIN_DIR ${PROJECT_BINARY_DIR}/exception CACHE INTERNAL "")
set(FAUST_MATRIX_BIN_DIR ${PROJECT_BINARY_DIR}/matrix CACHE INTERNAL "")
set(FAUST_FAUSTCORE_BIN_DIR ${PROJECT_BINARY_DIR}/faust_core CACHE INTERNAL "")
#set(FAUST_MATIO_BIN_DIR ${PROJECT_BINARY_DIR}/testing/init_from_matlab CACHE INTERNAL "")
set(FAUST_PALM4MSA_BIN_DIR ${PROJECT_BINARY_DIR}/palm4MSA CACHE INTERNAL "")
set(FAUST_TEST_BIN_DIR ${PROJECT_BINARY_DIR}/testing/bin CACHE INTERNAL "")
set(FAUST_TESTSRC_BIN_DIR ${PROJECT_BINARY_DIR}/testing/src CACHE INTERNAL "")
set(FAUST_TESTTMP_BIN_DIR ${PROJECT_BINARY_DIR}/testing/tmp CACHE INTERNAL "")
set(FAUST_TESTOUTPUT_BIN_DIR ${PROJECT_BINARY_DIR}/testing/output CACHE INTERNAL "")
set(FAUST_TESTRUNCOMP_BIN_DIR ${PROJECT_BINARY_DIR}/testing/runtime_comparison CACHE INTERNAL "")
set(FAUST_TESTTOOLS_BIN_DIR ${PROJECT_BINARY_DIR}/testing/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACE_BIN_DIR ${FAUST_WRAPPER_BIN_DIR}/mex CACHE INTERNAL "")
set(FAUST_MEXINTERFACETOOLS_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR}/tools CACHE INTERNAL "")
set(FAUST_MEXINTERFACEMEX_BIN_DIR ${FAUST_MEXINTERFACE_BIN_DIR} CACHE INTERNAL "") #  a MODIFIER ET METTRE LES MEXFONCTIONS.mexa DEDANS

#DEMO DIRECTORIES
set(FAUST_DEMO_SRC_DIR ${PROJECT_SOURCE_DIR}/demo CACHE INTERNAL "")
set(FAUST_DEMODATA_SRC_DIR ${PROJECT_SOURCE_DIR}/demo/data CACHE INTERNAL "")
set(FAUST_DEMOTOOLS_SRC_DIR ${PROJECT_SOURCE_DIR}/demo/tools CACHE INTERNAL "")
set(FAUST_DEMO_BIN_DIR ${PROJECT_BINARY_DIR}/demo CACHE INTERNAL "")
set(FAUST_DEMOTOOLS_BIN_DIR ${PROJECT_BINARY_DIR}/tools CACHE INTERNAL "")




#CMDLINE DIRECTORIES

set(FAUST_CMDLINE_SRC_DIR ${PROJECT_SOURCE_DIR}/cmd_line CACHE INTERNAL "")
set(FAUST_CMDLINE_TOOLS_SRC_DIR ${FAUST_CMDLINE_SRC_DIR}/xml_tools CACHE INTERNAL "")
set(FAUST_CMDLINE_FUNC_SRC_DIR ${FAUST_CMDLINE_SRC_DIR}/cmdline_function CACHE INTERNAL "")
set(FAUST_CMDLINE_BIN_SRC_DIR ${FAUST_WRAPPER_BIN_DIR}/cmd_line/src CACHE INTERNAL "")
set(FAUST_CMDLINE_BIN_DIR ${FAUST_WRAPPER_BIN_DIR}/cmd_line/bin CACHE INTERNAL "")


#DOCUMENTATION DIRECTORIES
set(FAUST_DOC_SRC_DIR ${PROJECT_SOURCE_DIR}/gen_doc CACHE INTERNAL "")
set(FAUST_DOC_BIN_DIR ${PROJECT_BINARY_DIR}/documentation CACHE INTERNAL "")



file(MAKE_DIRECTORY ${FAUST_TESTOUTPUT_BIN_DIR})
file(MAKE_DIRECTORY ${FAUST_TESTTMP_BIN_DIR})

# INSTALL DIRECTORIES (RELATIVES PATHS FROM CMAKE_INSTALL_PREFIX)
set(FAUST_INSTALL_MEX ${CMAKE_INSTALL_PREFIX}/mex)
set(FAUST_INSTALL_LIB ${CMAKE_INSTALL_PREFIX}/lib)
set(FAUST_INSTALL_CMDLINE ${CMAKE_INSTALL_PREFIX}/CmdLine)
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
set(FAUST_XML_TARGET faust_xml CACHE INTERNAL "")
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



# Default values for FAUST_USE_SINGLEPRECISION and FAUST_USE_OPENBLAS and USE_MATIO
set(FAUST_USE_SINGLEPRECISION ON CACHE BOOL "Using single precision instead of double precision for matrix and vector computations")
set(FAUST_USE_OPENBLAS ON CACHE BOOL "Using openBLAS for matrix and vector computations")
set(FAUST_USE_MEX OFF CACHE BOOL "Generate Mexfiles")
set(FAUST_USE_MATIO ON CACHE BOOL "Using matio library to read matfile")
set(FAUST_USE_XML ON CACHE BOOL "Using xml configuration file")
set(FAUST_ISVERBOSE OFF CACHE BOOL "Display useful message for debugging")
set(FAUST_GEN_DOC ON CACHE BOOL "Generate html documentation with doxygen")

# find external libraries

if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif()

include(CMake/add_library_path.cmake)
add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/opt" "/opt/local")
add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{MATIODIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/opt" "/opt/local" "/usr/include/libxml2")
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")
message(STATUS "LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")

include(CMake/check_external_libraries.cmake)
check_external_libraries(hdf5 HDF5_LIB_FILE 0)
check_external_libraries(matio MATIO_LIB_FILE 0)
check_external_libraries(xml2 XML2_LIB_FILE 0)
check_external_includes("libxml/parser.h" XML_INC_DIR 0)

if(FAUST_USE_OPENBLAS)
	check_external_libraries(openblas OPENBLAS_LIB_FILE 1)
	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
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

# if(FAUST_USE_SINGLEPRECISION)
	# set(CXX_MEX_FLAGS "${CXX_MEX_FLAGS} -DFAUST_SINGLE")
	# if(UNIX)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# elseif(WIN32)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# endif()
	# message(STATUS "**** SINGLE PRECISION USE *****")
# endif()

# test if executable matlab is in the path
if(UNIX)
message(STATUS "MATLAB_DIR_TMP 1 = ${MATLAB_DIR_TMP}")
	exec_program("which matlab | xargs echo" OUTPUT_VARIABLE MATLAB_DIR_TMP)
message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")
   exec_program("readlink ${MATLAB_DIR_TMP}" OUTPUT_VARIABLE READLINK_TMP)
message(STATUS "READLINK_TMP = ${READLINK_TMP}")
   if(${READLINK_TMP} MATCHES matlab)
      set(MATLAB_DIR_TMP ${READLINK_TMP})
      message(STATUS "MATLAB_DIR_TMP 3 = ${MATLAB_DIR_TMP}")
   endif()
   message(STATUS "MATLAB_DIR_TMP 4 = ${MATLAB_DIR_TMP}")

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
add_subdirectory(${FAUST_DOC_SRC_DIR})
add_subdirectory(${FAUST_CMDLINE_SRC_DIR} ${FAUST_CMDLINE_BIN_DIR})

if (FAUST_USE_MEX)
	add_subdirectory(${FAUST_MEXINTERFACE_SRC_DIR} ${FAUST_MEXINTERFACE_BIN_DIR})
	add_subdirectory(${FAUST_DEMO_SRC_DIR} ${FAUST_DEMO_BIN_DIR})
endif()

file(GLOB MATRIX_HPP_FILES "matrix/*.hpp")
if (FAUST_USE_MEX)
add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}> $<TARGET_OBJECTS:${FAUST_MEXTOOLS_TARGET}> )
else ()
add_library(${FAUST_TARGET} STATIC $<TARGET_OBJECTS:${FAUST_EXCEPTION_TARGET}> $<TARGET_OBJECTS:${FAUST_MATRIX_TARGET}> $<TARGET_OBJECTS:${FAUST_FAUSTCORE_TARGET}> $<TARGET_OBJECTS:${FAUST_PALM4MSA_TARGET}>)
endif()
#target_link_libraries(${FAUST_TARGET} ${MATLAB_ROOT}/bin/glnxa64/libmex.so )


if(FAUST_USE_OPENBLAS)
		target_link_libraries(${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
else()
		target_link_libraries(${FAUST_TARGET} )	
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









