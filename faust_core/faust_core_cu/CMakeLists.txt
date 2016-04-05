#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB HPP_FILES "*.hpp")

#find_package(CUDA QUIET REQUIRED)

#-D__COMPILE_TIMERS__
if(UNIX) # AND FAUST_USE_MKL
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
elseif(WIN32 AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()

add_definitions(-DEIGEN_MPL2_ONLY)
#-D__COMPILE_GPU__ -D__COMPILE_SPMAT__

include_directories(${FAUST_MATRIX_CU_SRC_DIR} ${FAUST_PALM4MSA_CU_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR})
add_library(${FAUST_FAUSTCORE_CU_TARGET} OBJECT ${HPP_FILES})
set_target_properties(${FAUST_FAUSTCORE_CU_TARGET} PROPERTIES LINKER_LANGUAGE CXX) #comme il n y a que des hpp, force la compilation en langage C++ pour eviter des erreurs et des warnings
 
if(FAUST_USE_OPENBLAS)
	include_directories( ${OPENBLAS_INC_DIR})
endif()
#target_link_libraries(${FAUST_FAUSTCORE_TARGET}  ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} )

