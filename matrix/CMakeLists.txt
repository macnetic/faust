#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB CPP_FILES "*.cpp")


if(UNIX AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
elseif(WIN32 AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()


add_definitions(-DEIGEN_MPL2_ONLY)

include_directories( ${FAUST_EXCEPTION_SRC_DIR} ${EIGEN_INC_DIR})
if(FAUST_USE_OPENBLAS)
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

if (FAUST_ISVERBOSE)
	set(FAUST_VERBOSE "1")
endif()
	
# permet de definir ou non les variable de compilation __GEMM_WITH_OPENBLAS__ et FAUST_SINGLE
configure_file(${FAUST_MATRIX_SRC_DIR}/faust_constant.h.in ${FAUST_MATRIX_SRC_DIR}/faust_constant.h @ONLY)




#message(STATUS "AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
add_library(${FAUST_MATRIX_TARGET} OBJECT ${CPP_FILES})

if(FAUST_USE_OPENBLAS)
	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} ${OPENBLAS_LIB_FILE})
else()
	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET})
endif()
