#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB CPP_FILES "*.cpp")

include_directories( $ENV{EIGEN_DIR})
if(FAUST_USE_OPENBLAS)
	include_directories( $ENV{OPENBLASDIR})
endif()

add_library(${FAUST_MATRIX_TARGET} SHARED ${CPP_FILES})

if(FAUST_USE_OPENBLAS)
	target_link_libraries(${FAUST_MATRIX_TARGET} ${OPENBLAS_LIB_PATH}/libopenblas.so)
endif()

