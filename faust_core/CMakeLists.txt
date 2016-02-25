#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB HPP_FILES "*.hpp")
if(UNIX AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
elseif(WIN32 AND FAUST_USE_MKL)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
endif()
add_definitions(-DEIGEN_MPL2_ONLY)
include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EIGEN_SRC_DIR})
add_library(${FAUST_FAUSTCORE_TARGET} OBJECT ${HPP_FILES})
set_target_properties(${FAUST_FAUSTCORE_TARGET} PROPERTIES LINKER_LANGUAGE CXX) #comme il n y a que des hpp, force la compilation en langage C++ pour eviter des erreurs et des warnings
 
if(FAUST_USE_OPENBLAS)
	include_directories( ${OPENBLAS_INC_DIR})
endif()
#target_link_libraries(${FAUST_FAUSTCORE_TARGET}  ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} )

