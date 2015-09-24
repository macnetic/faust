#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB CPP_FILES "*.cpp")

include_directories( ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_EIGEN_SRC_DIR})
if(FAUST_USE_OPENBLAS)
	include_directories( ${OPENBLAS_INC_DIR})
endif()

add_library(${FAUST_MATRIX_TARGET} OBJECT ${CPP_FILES})

if(FAUST_USE_OPENBLAS)
	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} ${OPENBLAS_LIB_FILE})
else()
	#target_link_libraries(${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET})
endif()

