#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")
file(GLOB CPP_FILES "*.cpp")

include_directories( $ENV{EIGEN_DIR})
add_library(${FAUST_MATRIX_TARGET} SHARED ${CPP_FILES})

