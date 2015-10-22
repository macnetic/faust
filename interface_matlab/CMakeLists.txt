
#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")


file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MEXINTERFACETOOLS_SRC_DIR}/*.cpp")
include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_MEXINTERFACETOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR} ${FAUST_EIGEN_SRC_DIR} )
add_library(${FAUST_MEXTOOLS_TARGET} OBJECT ${CPP_MEXTOOLS_FILES})
#target_link_libraries(${FAUST_MEXTOOLS_TARGET} ${FAUST_TARGET})


file(GLOB CPP_MEX_FILES "${FAUST_MEXINTERFACEMEX_SRC_DIR}/*.cpp")
string(REGEX REPLACE "[a-zA-Z0-9_/:]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")

if(CXX_MEX_FLAGS)
string(STRIP ${CXX_MEX_FLAGS} CXX_MEX_FLAGS)
string(REGEX REPLACE "([ \n\r\t]+)" ";" CMAKE_CXX_FLAGS_LIST "${CXX_MEX_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS_LIST} "-I${MATLAB_INCLUDE_DIR}" "-I${FAUST_EIGEN_SRC_DIR}" "-I${FAUST_MATRIX_SRC_DIR}" "-I${FAUST_PALM4MSA_SRC_DIR}" "-I${FAUST_FAUSTCORE_SRC_DIR}" "-I${FAUST_MEXINTERFACEMEX_SRC_DIR}" "-I${FAUST_EXCEPTION_SRC_DIR}" "-I${FAUST_MEXINTERFACE_SRC_DIR}" "-I${FAUST_MEXINTERFACETOOLS_SRC_DIR}")

if(FAUST_USE_OPENBLAS AND UNIX)
	set(CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS_LIST}" "-I${OPENBLAS_INC_DIR}" "-D__GEMM_WITH_OPENBLAS__")
endif()
#target_link_libraries(${FAUST_FAUSTCORE_TARGET}  ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} )


if(UNIX)
	set(LDD_MEX_FLAGS "-L${FAUST_BIN_DIR}" "-l${FAUST_TARGET}")
	if(FAUST_USE_OPENBLAS)
	   set(LDD_MEX_FLAGS "${LDD_MEX_FLAGS}" "${OPENBLAS_LIB_FILE}")
	endif()
	if(FAUST_USE_MKL)
		set(LDD_MEX_FLAGS "${LDD_MEX_FLAGS}"  "${MKL_LIB_CORE}" "${MKL_LIB_LP64}" "${MKL_LIB_M}" "${MKL_LIB_IOMP5}" "${MKL_LIB_INTEL_THREAD}")
	endif()
elseif(WIN32)

endif()

if(UNIX)
	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_TARGET_LIST} ${FAUST_TARGET} ${FAUST_MEXTOOLS_TARGET} )
elseif(WIN32)
	file(GLOB CPP_SRC_FILES "${FAUST_EXCEPTION_SRC_DIR}/*.cpp" "${FAUST_MATRIX_SRC_DIR}/*.cpp" "${FAUST_FAUSTCORE_SRC_DIR}/*.cpp" "${FAUST_PALM4MSA_SRC_DIR}/*.cpp"  "${FAUST_MEXINTERFACETOOLS_SRC_DIR}/*.cpp")
	string(REGEX REPLACE "([a-zA-Z0-9_/:]+)\\.cpp(;|$)" "\\1\\2" MEXSRC_FULL "${CPP_SRC_FILES}")
	string(REGEX REPLACE "[a-zA-Z0-9_/:]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "${FAUST_TMPMEX_DIR}/\\1.${MEXOBJ_EXT}\\2" MEXOBJ_FILES "${CPP_SRC_FILES}")

	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS MEXFILE_TARGET_LIST_TARGET ) 
	add_custom_target(MEXFILE_TARGET_LIST_TARGET DEPENDS  ${MEXFILE_TARGET_LIST})
	add_custom_target(MEXOBJ_FILES_TARGET DEPENDS ${MEXOBJ_FILES})
	add_dependencies(MEXFILE_TARGET_LIST_TARGET MEXOBJ_FILES_TARGET)

endif()

foreach(mex_target ${MEXFILE_TARGET_LIST})
	if(WIN32)
		add_custom_target(${mex_target} DEPENDS ${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} MEXOBJ_FILES_TARGET )

	elseif(UNIX)
		add_custom_target(${mex_target} DEPENDS ${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} ${FAUST_MEXTOOLS_TARGET} ${FAUST_TARGET})
	endif()


	if(UNIX)
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		COMMAND   ${MATLAB_ROOT}/bin/mex 
		ARGS "-v" "-largeArrayDims" "${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp" 
		"LDFLAGS=$$LDFLAGS -Wl,-rpath=${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}:${FAUST_INSTALL_LIB}:${PROJECT_BINARY_DIR}" 
		"-O" 
		${CMAKE_CXX_FLAGS_LIST} 
		${LDD_MEX_FLAGS}
		COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
	elseif(WIN32)
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		COMMAND   ${MATLAB_ROOT}/bin/mex 
		ARGS "-v" "-largeArrayDims" "${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp" "${FAUST_TMPMEX_DIR}/*.${MEXOBJ_EXT}" 
		${CMAKE_CXX_FLAGS_LIST}
		DEPENDS ${MEXOBJ_FILES}
		COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")

	else()
		message(FATAL_ERROR "Unable to compile mex functions")
	endif()
	install(FILES ${FAUST_MEXINTERFACE_BIN_DIR}/${mex_target}.${MEX_EXT} DESTINATION ${FAUST_INSTALL_MEX})


endforeach()

if(WIN32)
	foreach( mex_src ${MEXSRC_FULL})
		string(REGEX REPLACE "[a-zA-Z0-9_/:]+/([a-zA-Z0-9_]+)" "\\1" mexsrc_short "${mex_src}")
		add_custom_command(OUTPUT ${FAUST_TMPMEX_DIR}/${mexsrc_short}.${MEXOBJ_EXT}
			COMMAND  ${MATLAB_ROOT}/bin/mex 
			ARGS "-v" "-c" "-largeArrayDims" ${mex_src}.cpp  ${CMAKE_CXX_FLAGS_LIST} "-outdir" "${FAUST_TMPMEX_DIR}" "OBJEXT=.${MEXOBJ_EXT}"
			COMMENT "Compiling ${mex_src}.cpp for mex usage")
	endforeach()
endif()


