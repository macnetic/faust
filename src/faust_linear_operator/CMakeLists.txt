#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")

file(GLOB CPP_FILES "*.cpp" "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp")
file(GLOB HPP_FILES "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.hpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.hpp")


if(USE_GPU_MOD)
	include_directories(${GPU_MOD_INCLUDE_DIR} ${GPU_MOD_INCLUDE_DIR}/../build)
	include_directories(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR})
endif()
include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${FAUST_UTILS_SRC_DIR} ${EIGEN_INC_DIR})

if (BUILD_USE_SINGLEPRECISION)
	set(FAUST_SINGLE "1")
endif()

if (BUILD_USE_PROFILING)
	set(__COMPILE_TIMERS__ "1")
endif()

if (BUILD_VERBOSE)
	set(FAUST_VERBOSE "1")
endif()
	
# permet de definir ou non la variable de compilation FAUST_SINGLE
configure_file(${FAUST_SRC_LINEAR_OPERATOR_DIR}/faust_constant.h.in ${FAUST_SRC_LINEAR_OPERATOR_DIR}/faust_constant.h @ONLY)

#message(STATUS "AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
#add_library(${FAUST_FAUSTCORE_TARGET} OBJECT ${HPP_FILES})
#add_library(${FAUST_MATRIX_TARGET} OBJECT ${CPP_FILES})

add_library(${FAUST_LINEAR_OPERATOR_TARGET} OBJECT ${CPP_FILES} ${HPP_FILES})



#set_target_properties(${FAUST_LINEAR_OPERATOR_TARGET} PROPERTIES LINKER_LANGUAGE CXX)
#set_target_properties(${FAUST_FAUSTCORE_TARGET} PROPERTIES LINKER_LANGUAGE CXX) #comme il n y a que des hpp, force la compilation en  langage C++ pour eviter des erreurs et des warnings





