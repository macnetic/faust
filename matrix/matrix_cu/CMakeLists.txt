#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB CPP_FILES "*.cpp")
file(GLOB H_FILES "*.h" "*.hpp*")
file(GLOB CU_FILES "*.cu")
file(GLOB HU_FILES "*.hu")


find_package(CUDA QUIET REQUIRED)

if(UNIX) # AND FAUST_USE_MKL
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp ")
	#-D__COMPILE_GPU__ -D__COMPILE_SPMAT__
	
	# Pass options to NVCC (-O3)    
	set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3 -Xcompiler -fPIC -use_fast_math  -gencode arch=compute_20,code=compute_20 -gencode arch=compute_30,code=compute_30 -gencode arch=compute_35,code=compute_35 ")
    #-D__COMPILE_GPU__ -D__COMPILE_SPMAT__
	
elseif(WIN32 AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()

#message(STATUS "CUDA_NVCC_FLAGS=${CUDA_NVCC_FLAGS}")
#message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

add_definitions(-DEIGEN_MPL2_ONLY)

include_directories( ${FAUST_MATRIX_CU_SRC_DIR} ${FAUST_PALM4MSA_CU_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR})

#if(FAUST_USE_OPENBLAS)
#	include_directories( ${OPENBLAS_INC_DIR} )
#	set(__GEMM_WITH_OPENBLAS__ "1")
#	#message(STATUS "**** OPENBLAS USE *****")
#endif()


if (FAUST_USE_SINGLEPRECISION)
	set(FAUST_SINGLE "1")
endif()

if (FAUST_USE_PROFILING)
	set(__COMPILE_TIMERS__ "1")
endif()

if (FAUST_ISVERBOSE)
	set(FAUST_VERBOSE "1")
endif()
if (FAUST_USE_GPU)
	set(__COMPILE_GPU__ "1")
	set(__COMPILE_SPMAT__ "1")
endif()

# permet de definir ou non les variable de compilation __GEMM_WITH_OPENBLAS__ et FAUST_SINGLE
configure_file(${FAUST_MATRIX_CU_SRC_DIR}/faust_constant.h.in ${FAUST_MATRIX_CU_SRC_DIR}/faust_constant.h @ONLY)



#cuda_add_library(${FAUST_MATRIX_CU_TARGET} ${CPP_FILES} ${CU_FILES} )

cuda_add_library(${FAUST_MATRIX_CU_TARGET} ${CPP_FILES} ${CU_FILES} )

    
#if(FAUST_USE_OPENBLAS)
##	target_link_libraries(${FAUST_MATRIX_CU_TARGET} ${FAUST_EXCEPTION_TARGET} ${OPENBLAS_LIB_FILE})
#else()
##	target_link_libraries(${FAUST_MATRIX_CU_TARGET} ${FAUST_EXCEPTION_TARGET})
#endif()


#message(STATUS "AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
#add_library(${FAUST_MATRIX_CU_TARGET} OBJECT ${CPP_FILES} ${CU_FILES})

#if(FAUST_USE_OPENBLAS)
#	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} ${OPENBLAS_LIB_FILE})
#else()
#	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET})
#endif()

# Pass options to NVCC
#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O3 -m64 -Xcompiler -fPIC -use_fast_math  -gencode arch=compute_20,code=sm_20 \
#-gencode arch=compute_20,code=sm_21 \
#-gencode arch=compute_30,code=sm_30 \
#-gencode arch=compute_35,code=sm_35 \
#-gencode arch=compute_50,code=sm_50 \
#-D__COMPILE_GPU__ -D__COMPILE_SPMAT__ -D__COMPILE_TIMERS__
#    )

#MATRIX_DIR = ../../matrix/matrix_cu
#PALM_DIR = ../../palm4MSA/palm4MSA_cu
#FAUSTCORE_DIR = ../../faust_core/faust_core_cu
#MATIO_DIR = ../../init_from_matlab
#EXCEPTION_DIR = ../../exception
#TEST_DIR = .
#LIB_MATRIX = ${PROJECT_SOURCE_DIR}/matrix/matrix_cu/libfaustmat.a
#LIB_PALM = ${PROJECT_SOURCE_DIR}/palm4MSA/palm4MSA_cu/libfaustpalm.a
#LIB_MATIO = ${PROJECT_SOURCE_DIR}/init_from_matlab/libfaustmatio.a
#LIB_EXCEPTION = ${PROJECT_SOURCE_DIR}/exception/libfaustexeption.a

# Specify library paths
#link_directories(
#	${PROJECT_SOURCE_DIR}/matrix/matrix_cu
#	${PROJECT_SOURCE_DIR}/palm4MSA/palm4MSA_cu
#	${PROJECT_SOURCE_DIR}/init_from_matlab
#	${PROJECT_SOURCE_DIR}/exception
#)



# For compilation ...
# Specify target & source files to compile it from
#cuda_add_executable(
#    hellocuda
#    hellocuda.cu
#    hellocuda.h
#    kernels/hellokernels.cu
#    kernels/hellokernels.h
#    utility/wrapper.cpp
#    utility/wrapper.h
#    )

# For linking ...
# Specify target & libraries to link it with
#target_link_libraries(
#    hellocuda
#    -libfaustmat
#    -libfaustpalm
#    -libfaustmatio
#    -libfaustexeption
#    )
    
    
    
    
    
