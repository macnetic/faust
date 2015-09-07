
#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")


file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MEXINTERFACETOOLS_SRC_DIR}/*.cpp")
include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_MEXINTERFACETOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR} $ENV{EIGEN_DIR})
add_library(${FAUST_MEXTOOLS_TARGET} ${CPP_MEXTOOLS_FILES})
target_link_libraries(${FAUST_MEXTOOLS_TARGET} ${FAUST_FAUSTCORE_TARGET} ${FAUST_PALM4MSA_TARGET}  ${FAUST_MATRIX_TARGET} )


file(GLOB CPP_MEX_FILES "${FAUST_MEXINTERFACEMEX_SRC_DIR}/*.cpp")
string(REGEX REPLACE "\\.cpp(;|$)" ".mexa64\\1" MEXFILE_LIST "${CPP_MEX_FILES}")


add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_LIST} ${FAUST_FAUSTCORE_TARGET} ${FAUST_MEXTOOLS_TARGET} ${FAUST_MATRIX_TARGET} ${FAUST_PALM4MSA_TARGET})

string(STRIP ${CMAKE_CXX_FLAGS} CMAKE_CXX_FLAGS)
string(REGEX REPLACE "( )+" ";" CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS}")

foreach(mexfile ${MEXFILE_LIST})

	string(REGEX REPLACE "\\.mexa64" ".cpp" src_tmp "${mexfile}")
	string(REGEX REPLACE "\\.mexa64" ".o" obj_tmp "${mexfile}")

	add_custom_command(OUTPUT ${mexfile}
		COMMAND ${CMAKE_CXX_COMPILER} 
		ARGS "-c" ${CMAKE_CXX_FLAGS_LIST} "-DNDEBUG -DMATLAB_MEX_FILE" 
		"-I${MATLAB_INCLUDE_DIR}" 
		"-I${FAUST_MATRIX_SRC_DIR}" 
		"-I${FAUST_PALM4MSA_SRC_DIR}" 
		"-I$ENV{EIGEN_DIR}" 
		"-I${FAUST_FAUSTCORE_SRC_DIR}" 
		"-I${FAUST_MEXINTERFACEMEX_SRC_DIR}" 
		"-I${FAUST_MEXINTERFACE_SRC_DIR}" 
		"-I${FAUST_MEXINTERFACETOOLS_SRC_DIR}"
		"${src_tmp}" "-o" "${obj_tmp}"

        	COMMAND ${CMAKE_CXX_COMPILER} 
		ARGS "-pthread" "-Wl,--no-undefined" "-fPIC" "-shared" "-O3" "-Wl,--version-script,${MATLAB_ROOT}/extern/lib/glnxa64/mexFunction.map" 
		"-lstdc++" "-lm"
		"-L${FAUST_FAUSTCORE_BIN_DIR}" "-l${FAUST_FAUSTCORE_TARGET}"
		"-L${FAUST_MATRIX_BIN_DIR}" "-l${FAUST_MATRIX_TARGET}"
		"-L${FAUST_PALM4MSA_BIN_DIR}" "-l${FAUST_PALM4MSA_TARGET}"
		"-L${FAUST_MEXINTERFACE_BIN_DIR}" "-l${FAUST_MEXTOOLS_TARGET}"
		"-L${OPENBLAS_LIB}" "-lopenblas"
		"-L${MATLAB_ROOT}/bin/glnxa64" "-lmx" "-lmex" "-lmat"	
		"${obj_tmp}" "-o" "${mexfile}"
                   COMMENT "Creating mex function ${mexfile}"

		COMMAND rm
		ARGS -f ${obj_tmp} VERBATIM)

endforeach()

