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
CHECK_CXX_COMPILER_FLAG(-std=c++11 COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG(-std=c++0x COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(FATAL_ERROR "The compiler ${CMAKE_CXX_COMPILER} has neither C++11 nor c++0x support.")
endif()


INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")


add_subdirectory(${FAUST_MATRIX_SRC_DIR} ${FAUST_MATRIX_BIN_DIR})
add_subdirectory(${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_FAUSTCORE_BIN_DIR})
add_subdirectory(${FAUST_MATIO_SRC_DIR} ${FAUST_MATIO_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSA_SRC_DIR} ${FAUST_PALM4MSA_BIN_DIR})
add_subdirectory(${FAUST_PALM4MSATEST_SRC_DIR} ${FAUST_PALM4MSATEST_BIN_DIR})
add_subdirectory(${FAUST_MEXINTERFACE_SRC_DIR} ${FAUST_MEXINTERFACE_BIN_DIR})
add_subdirectory(${FAUST_MEXINTERFACETOOLS_SRC_DIR} ${FAUST_MEXINTERFACETOOLS_BIN_DIR})


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








