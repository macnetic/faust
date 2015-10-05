
#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")


file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MEXINTERFACETOOLS_SRC_DIR}/*.cpp")
include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_MEXINTERFACETOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR} ${FAUST_EIGEN_SRC_DIR} )
add_library(${FAUST_MEXTOOLS_TARGET} OBJECT ${CPP_MEXTOOLS_FILES})
#target_link_libraries(${FAUST_MEXTOOLS_TARGET} ${FAUST_TARGET})


file(GLOB CPP_MEX_FILES "${FAUST_MEXINTERFACEMEX_SRC_DIR}/*.cpp")
string(REGEX REPLACE "[a-zA-Z0-9_/:]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")

add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_TARGET_LIST} ${FAUST_TARGET} ${FAUST_MEXTOOLS_TARGET} )

string(STRIP ${CMAKE_CXX_FLAGS} CMAKE_CXX_FLAGS)
string(REGEX REPLACE "([ \n\r\t]+)" ";" CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS}")


set(CXX_MEX_FLAGS "-I${MATLAB_INCLUDE_DIR}" "-I${FAUST_EIGEN_SRC_DIR}" "-I${FAUST_MATRIX_SRC_DIR}" "-I${FAUST_PALM4MSA_SRC_DIR}" "-I${FAUST_FAUSTCORE_SRC_DIR}" "-I${FAUST_MEXINTERFACEMEX_SRC_DIR}" "-I${FAUST_MEXINTERFACE_SRC_DIR}" "-I${FAUST_MEXINTERFACETOOLS_SRC_DIR}")

if(UNIX)
	set(LDD_MEX_FLAGS "-L${FAUST_LIB_DIR}"  "-L${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}" "-l${FAUST_TARGET}"  "-lmx" "-lmex" "-lmat")
	if(FAUST_USE_OPENBLAS)
	   set(LDD_MEX_FLAGS "${LDD_MEX_FLAGS}" "${OPENBLAS_LIB_FILE}")
	endif()
	if(FAUST_USE_MKL)
		set(LDD_MEX_FLAGS "${LDD_MEX_FLAGS}"  "${MKL_LIB_CORE}" "${MKL_LIB_LP64}" "${MKL_LIB_M}" "${MKL_LIB_IOMP5}" "${MKL_LIB_INTEL_THREAD}")
	endif()
elseif(WIN32)
	set(LDD_MEX_FLAGS "/LIBPATH:${FAUST_LIB_DIR}" "/LIBPATH:${MATLAB_ROOT}/extern/lib/${MEX_SUBDIR_LIB}/microsoft" "${FAUST_TARGET}.lib" "libmx.lib" "libmex.lib" "libmat.lib" "kernel32.lib" "user32.lib" "gdi32.lib" "winspool.lib" "comdlg32.lib" "advapi32.lib" "shell32.lib" "ole32.lib" "oleaut32.lib" "uuid.lib" "odbc32.lib" "odbccp32.lib")
endif()


foreach(mex_target ${MEXFILE_TARGET_LIST})

	add_custom_target(${mex_target} ALL DEPENDS ${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} ${FAUST_MEXTOOLS_TARGET} ${FAUST_TARGET})

	if(UNIX)
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		COMMAND   ${MATLAB_ROOT}/bin/mex 
		ARGS "-v" "-largeArrayDims" "${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp" 
		"LDFLAGS=$$LDFLAGS -Wl,-rpath=${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}:${FAUST_MEXINTERFACE_LIB_DIR}" 
		"-O" "${FAUST_COMPILE_SINGLE_FLAG}"
		${CXX_MEX_FLAGS} 
		${LDD_MEX_FLAGS}
		COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
	elseif(WIN32)
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		COMMAND cl 
		ARGS "/c"  "/GR" "/W3" "/EHs" "/nologo" "/MD" "/O2" "/Oy-"
		"/DFAUST_SINGLE" "/D_CRT_SECURE_NO_DEPRECATE" "/D_SCL_SECURE_NO_DEPRECATE" "/D_SECURE_SCL=0" "/DMATLAB_MEX_FILE" "/DNDEBUG"
		${CXX_MEX_FLAGS}	
		"${FAUST_MEXINTERFACEMEX_SRC_DIR}/${mex_target}.cpp"
		"/Fo${FAUST_MEXINTERFACEMEX_BIN_DIR}/${mex_target}.obj"
	
		COMMAND link  
		ARGS "/nologo" "/manifest" "/incremental:NO" "/DLL" "/EXPORT:mexFunction" 
		"${FAUST_MEXINTERFACEMEX_BIN_DIR}/${mex_target}.obj"
		${LDD_MEX_FLAGS}
		"/out:${mex_target}.${MEX_EXT}"
	
		COMMAND mt 
		ARGS "-outputresource:${mex_target}.${MEX_EXT};2" "-manifest" "${mex_target}.${MEX_EXT}.manifest"
		
		COMMAND del 
		ARGS "${mex_target}.exp" "${mex_target}.lib" "${mex_target}.${MEX_EXT}.manifest"
	
		COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
	else()
		message(FATAL_ERROR "Unable to compile mex functions")
	endif()


endforeach()

