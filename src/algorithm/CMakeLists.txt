##INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
#set(SRC_FILES "")
##file(GLOB CPP_FILES "*.cpp")

#file(GLOB CPP_FILES "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*.cpp")
##message(STATUS "algorithm cpp files ${CPP_FILES}")

#### include directory 
#include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${EIGEN_INC_DIR})
#if (BUILD_USE_GPU)
#	include_directories(${FAUST_SRC_LINEAR_OPERATOR_GPU_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR})
#endif()
#if(BUILD_OPENBLAS)
#	include_directories(${OPENBLAS_INC_DIR})
#endif()


#if(APPLE AND ( CMAKE_GENERATOR STREQUAL "Xcode") )
#	# on compile les sources directement dans le cmaklist principal pour Xcode 
#	#add_library(${FAUST_ALGORITHM_TARGET} ${CPP_FILES})
#else ()
#	add_library(${FAUST_ALGORITHM_TARGET} OBJECT ${CPP_FILES})
#endif()

##target_link_libraries(${FAUST_PALM4MSA_TARGET} ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET})
##target_link_libraries(${FAUST_PALM4MSA_CU_TARGET} ${FAUST_MATRIX_CU_TARGET} ${FAUST_EXCEPTION_TARGET})

foreach(REAL_TYPE IN ITEMS float double)
  if(REAL_TYPE MATCHES double)
    set(SVD_COMPLEX 1)
  else()
    set(SVD_COMPLEX 0)
  endif()
  configure_file(${FAUST_ALGORITHM_SRC_DIR}/factorization/faust_SVDTJ.cpp.in ${FAUST_ALGORITHM_SRC_DIR}/factorization/faust_SVDTJ_${REAL_TYPE}.cpp)
endforeach()

