#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")

file(GLOB CPP_FILES "*.cpp" "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp")
file(GLOB HPP_FILES "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.hpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.hpp")


#if(UNIX AND FAUST_USE_MKL)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#elseif(WIN32 AND FAUST_USE_MKL)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
#endif()
#message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
#add_definitions(-DEIGEN_MPL2_ONLY)

#if(BUILD_OPENBLAS)
#	include_directories( ${OPENBLAS_INC_DIR})
#endif()
#target_link_libraries(${FAUST_FAUSTCORE_TARGET}  ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} )

#include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EIGEN_SRC_DIR})
include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${EIGEN_INC_DIR})

if(BUILD_OPENBLAS)
	#message(STATUS "**********OPENBLAS ACTIF************")
	#message(STATUS "OPENBLAS_INC_DIR ${OPENBLAS_INC_DIR}")
	include_directories( ${OPENBLAS_INC_DIR})
	set(__GEMM_WITH_OPENBLAS__ "1")
	#message(STATUS "**** OPENBLAS USE *****")
endif()

if (FAUST_USE_SINGLEPRECISION)
	set(FAUST_SINGLE "1")
endif()

if (FAUST_USE_PROFILING)
	set(__COMPILE_TIMERS__ "1")
endif()

if (BUILD_VERBOSE)
	set(FAUST_VERBOSE "1")
endif()
	
# permet de definir ou non les variable de compilation __GEMM_WITH_OPENBLAS__ et FAUST_SINGLE
configure_file(${FAUST_SRC_LINEAR_OPERATOR_DIR}/faust_constant.h.in ${FAUST_SRC_LINEAR_OPERATOR_DIR}/faust_constant.h @ONLY)

#message(STATUS "AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
#add_library(${FAUST_FAUSTCORE_TARGET} OBJECT ${HPP_FILES})
#add_library(${FAUST_MATRIX_TARGET} OBJECT ${CPP_FILES})

add_library(${FAUST_LINEAR_OPERATOR_TARGET} OBJECT ${CPP_FILES} ${HPP_FILES})


#set_target_properties(${FAUST_LINEAR_OPERATOR_TARGET} PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(${FAUST_FAUSTCORE_TARGET} PROPERTIES LINKER_LANGUAGE CXX) #comme il n y a que des hpp, force la compilation en  langage C++ pour eviter des erreurs et des warnings





