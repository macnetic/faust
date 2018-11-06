# CMakeLists used to generated ".mexa64" files. 
# 

#INCLUDE("${PROJECT_SOURCE_DIR}/CMake/define_variables.cmake")

# copy and install the matlab wrapping the underlying C++ class

#file(MAKE_DIRECTORY ${FAUST_MATLAB_TOOLS_INSTALL_DIR})

# matlab tools directory
foreach(MATLAB_TOOLS FaustCore faust_decompose generate_params)
	configure_file(${FAUST_MATLAB_TOOLS_SRC_DIR}/${MATLAB_TOOLS}.m ${FAUST_MATLAB_TOOLS_BIN_DIR}/${MATLAB_TOOLS}.m COPYONLY)
	#install(FILES ${FAUST_MATLAB_TOOLS_BIN_DIR}/${MATLAB_TOOLS}.m DESTINATION ${FAUST_MATLAB_TOOLS_INSTALL_DIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
endforeach()

file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.cpp")
# include directory 
#{${FAUST_EXCEPTION_SRC_DIR}
include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}  ${FAUST_MATLAB_TOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR} ${EIGEN_INC_DIR} )

if(MATIO_INC_DIR)
	include_directories(${MATIO_INC_DIR})
endif()

#modif compatibilite ubuntu avec openblas
if (BUILD_OPENBLAS)
	include_directories(${OPENBLAS_INC_DIR})
endif(BUILD_OPENBLAS)

#message(FATAL_ERROR ${FAUST_MATLAB_MATFAUST_FILE_LIST})
foreach(MATFAUST_FILE IN LISTS FAUST_MATLAB_MATFAUST_FILE_LIST)
	string(REGEX REPLACE "[^+]*\\+matfaust\\/" "" MATFAUST_FILE ${MATFAUST_FILE})
	get_filename_component(PARENT_DIR ${MATFAUST_FILE} DIRECTORY)
	file(MAKE_DIRECTORY ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${PARENT_DIR})
	message(STATUS "Creating folder ${PARENT_DIR}")
	# use script generating inline doc only for real matlab classes
	if(${MATFAUST_FILE} MATCHES "@")
		message(STATUS "matfaust class: " ${MATFAUST_FILE})
		if(WIN32)
			execute_process(COMMAND git-bash ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} ERROR_QUIET)
			#configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} COPYONLY)
		else()
			execute_process(COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} ERROR_QUIET)
			#TODO: copy if error (fallback solution)
		endif()
	else()
		message(STATUS "matfaust function file: " ${MATFAUST_FILE})
		# just copy the .m
		configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} COPYONLY)
	endif()

endforeach()

add_library(${FAUST_MEXTOOLS_TARGET} OBJECT ${CPP_MEXTOOLS_FILES})


foreach(SCALAR_AND_FSUFFIX double:Real std::complex<double>:Cplx) # TODO: float and complex<float> versions
	# configuration de mexFaust.cpp.in
	# specialization of the class mexFaust template scalar with real or complex
	string(REGEX REPLACE "(.*):.*" "\\1" SCALAR ${SCALAR_AND_FSUFFIX})
	string(REGEX REPLACE ".*:(.*)" "\\1" FSUFFIX ${SCALAR_AND_FSUFFIX})
	set(FAUST_SCALAR ${SCALAR})
	set(FAUST_FPP double)
	set(FACT_FPP double) # TODO: we could allow float also
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexFaust.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexFaust${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPalm4MSA.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPalm4MSA${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical_fact.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical_fact${FSUFFIX}.cpp @ONLY)
	# copy the *.m for factorization now, because we have the FSUFFIX in hands
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact${FSUFFIX}.m COPYONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA${FSUFFIX}.m COPYONLY)
endforeach()

file(GLOB CPP_MEX_FILES "${FAUST_MATLAB_MEX_SRC_DIR}/*.cpp")
## modif AL : bug pour les noms de PATH avec un point comme par exemple : ./faust_v1.10/interface_matlab/
## On va chercher les nom des target
#string(REGEX REPLACE "[a-zA-Z0-9_/:.]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")
string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")

if(CXX_MEX_FLAGS)
	string(STRIP ${CXX_MEX_FLAGS} CXX_MEX_FLAGS)
	string(REGEX REPLACE "([ \n\r\t]+)" ";" CMAKE_CXX_FLAGS_LIST "${CXX_MEX_FLAGS}")
endif()

#"-I${FAUST_FAUSTCORE_SRC_DIR}" "-I${FAUST_EXCEPTION_SRC_DIR}"
set(CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS_LIST} "-I${MATIO_INC_DIR}" "-I${MATLAB_INCLUDE_DIR}" "-I${EIGEN_INC_DIR}" "-I${FAUST_SRC_LINEAR_OPERATOR_DIR}" "-I${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}" "-I${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}" "-I${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}" "-I${FAUST_MATLAB_TOOLS_SRC_DIR}")

#if(BUILD_OPENBLAS AND UNIX)
if(BUILD_OPENBLAS)
	set(CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS_LIST}" "-I${OPENBLAS_INC_DIR}")
endif()
#target_link_libraries(${FAUST_FAUSTCORE_TARGET}  ${FAUST_MATRIX_TARGET} ${FAUST_EXCEPTION_TARGET} )

if(UNIX)
    # we always use matio
    if(${USE_MATIO_STATIC_LIBS})
	    set(LDD_MEX_FLAGS ${LDD_MEX_FLAGS} "-L${FAUST_BIN_DIR}" "-l${FAUST_TARGET}" "${MATIO_STATIC_LIB_PATH}" "${HDF5_STATIC_LIB_PATH}" "${Z_STATIC_LIB_PATH}")
    else()
	    string(REGEX REPLACE "^(.+)/.*$" "-L\\1" LDD_MEX_FLAGS "${MATIO_LIB_FILE}")
	    set(LDD_MEX_FLAGS ${LDD_MEX_FLAGS} "-L${FAUST_BIN_DIR}" "-l${FAUST_TARGET}" "-lmatio")
    endif()
elseif(WIN32)
	# MODIF ALALAL
	#set(LDD_MEX_FLAGS "-L${FAUST_BIN_DIR}" "-l${FAUST_TARGET}")
       #windows version uses static lib of matio
       string(REGEX REPLACE "^(.+)/.*$" "\\1/libmatio.a" LDD_MEX_FLAGS "${MATIO_LIB_FILE}")
endif()

if(BUILD_OPENBLAS)
	set(LDD_MEX_FLAGS "${LDD_MEX_FLAGS}" "${OPENBLAS_LIB_FILE}")
endif()

if(UNIX)
	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_TARGET_LIST} ${FAUST_TARGET} ${FAUST_MEXTOOLS_TARGET} )
elseif(WIN32)
	file(GLOB CPP_SRC_FILES "${FAUST_SRC_LINEAR_OPERATOR_DIR}/*.cpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*.cpp"  "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.cpp")
	string(REGEX REPLACE "([a-zA-Z0-9_/:-]+)\\.cpp(;|$)" "\\1\\2" MEXSRC_FULL "${CPP_SRC_FILES}")
	
	#message (STATUS "MEXSRC_FULL ------------- ${MEXSRC_FULL}")
	#message (STATUS "FAUST_TMPMEX_DIR ------------- ${FAUST_TMPMEX_DIR}")
	#message (STATUS "CPP_SRC_FILES ---------${CPP_SRC_FILES}")
	#modif AL
	#string(REGEX REPLACE "[a-zA-Z0-9_/:-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "${FAUST_TMPMEX_DIR}/\\1.${MEXOBJ_EXT}\\2" MEXOBJ_FILES "${CPP_SRC_FILES}")
	string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "${FAUST_TMPMEX_DIR}/\\1.${MEXOBJ_EXT}\\2" MEXOBJ_FILES "${CPP_SRC_FILES}")
	#message (STATUS "MEXOBJ_FILES --------- ${MEXOBJ_FILES}")
	
	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS MEXFILE_TARGET_LIST_TARGET ) 
	add_custom_target(MEXFILE_TARGET_LIST_TARGET DEPENDS ${MEXFILE_TARGET_LIST})
	add_custom_target(MEXOBJ_FILES_TARGET DEPENDS ${MEXOBJ_FILES})
	add_dependencies(MEXFILE_TARGET_LIST_TARGET MEXOBJ_FILES_TARGET)
endif()

if(WIN32)
	foreach( mex_src ${MEXSRC_FULL})
		#modif AL pour gerer les "-" dans le path
		#string(REGEX REPLACE "[a-zA-Z0-9_/:-.]+/([a-zA-Z0-9_]+)" "\\1" mexsrc_short "${mex_src}")
		string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)(;|$)" "\\1\\2" mexsrc_short "${mex_src}")
		#message(STATUS "-------------------------mexsrc_short---${mexsrc_short}")
		add_custom_command(	OUTPUT ${FAUST_TMPMEX_DIR}/${mexsrc_short}.${MEXOBJ_EXT}
							COMMAND  ${MATLAB_ROOT}/bin/mex 
							ARGS "-v" "-c" "-largeArrayDims" ${mex_src}.cpp 
								${CMAKE_CXX_FLAGS_LIST} ${LDD_MEX_FLAGS}
								"-outdir" "${FAUST_TMPMEX_DIR}" "OBJEXT=.${MEXOBJ_EXT}"
							COMMENT "Compiling ${mex_src}.cpp for mex usage")
	#message (STATUS "${FAUST_TMPMEX_DIR}/${mexsrc_short}.${MEXOBJ_EXT}")
	endforeach()
endif()

#file(MAKE_DIRECTORY ${FAUST_MATLAB_MEX_INSTALL_DIR})
#file(MAKE_DIRECTORY ${CMAKE_INSTALL_MATLAB_PREFIX})		
foreach(mex_target ${MEXFILE_TARGET_LIST})
	if(WIN32)
		add_custom_target(${mex_target} DEPENDS ${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} MEXOBJ_FILES_TARGET )
	elseif(UNIX)
		add_custom_target(${mex_target} DEPENDS ${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} ${FAUST_MEXTOOLS_TARGET} ${FAUST_TARGET})
	endif()


	if(UNIX)
		if(APPLE)
			add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		        COMMAND   ${MATLAB_ROOT}/bin/mex 
		        ARGS "-v" "-largeArrayDims" "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" 
		        "-O" 
		        ${CMAKE_CXX_FLAGS_LIST} 
		        ${LDD_MEX_FLAGS}
		        COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
		else(APPLE)
			add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
				COMMAND   ${MATLAB_ROOT}/bin/mex 
				ARGS "-v" "-largeArrayDims" "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" 
					"LDFLAGS=$$LDFLAGS -Wl,-rpath=${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}:${FAUST_INSTALL_LIB}:${PROJECT_BINARY_DIR}" 
					"-O" 
					${CMAKE_CXX_FLAGS_LIST} 
					${LDD_MEX_FLAGS}
				COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
		endif(APPLE)

	elseif(WIN32)		
		#message(STATUS "${MEXOBJ_FILES} ")
		file(GLOB OBJ_MEXTOOLS_FILES "${FAUST_TMPMEX_DIR}/*.${MEXOBJ_EXT}")
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
					COMMAND   ${MATLAB_ROOT}/bin/mex 
					ARGS "-v" "-largeArrayDims" "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" ${OBJ_MEXTOOLS_FILES} ${MEXOBJ_FILES}
						${CMAKE_CXX_FLAGS_LIST} ${LDD_MEX_FLAGS}
					DEPENDS ${MEXOBJ_FILES} 
					COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}"
		)
		#message (STATUS "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp")
	else()
		message(FATAL_ERROR "Unable to compile mex functions for this platform")
	endif()
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/${mex_target}.m ${FAUST_MATLAB_MEX_BIN_DIR}/${mex_target}.m COPYONLY)
	#install(FILES ${FAUST_MATLAB_MEX_BIN_DIR}/${mex_target}.m DESTINATION ${FAUST_MATLAB_MEX_INSTALL_DIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)	
	#install(FILES ${FAUST_MATLAB_MEX_BIN_DIR}/${mex_target}.${MEX_EXT} DESTINATION ${FAUST_MATLAB_MEX_INSTALL_DIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
	
	endforeach()


		# CONFIGURATION DES CHEMINS
	configure_file(${FAUST_MATLAB_SRC_DIR}/setup_FAUST.m ${FAUST_MATLAB_BIN_DIR}/setup_FAUST.m COPYONLY)



	### OLD MATLAB VERSION CODE TO MAKE TIME COMPARISON ########
	#version MATLAB du code (entierement MATLAB, pas de mexfunction)

	set(FAUST_MATLAB_OLD_SRC_DIR ${FAUST_MATLAB_SRC_DIR}/old_matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_BIN_DIR ${FAUST_MATLAB_BIN_DIR}/old_matlab CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_PROX_BIN_DIR ${FAUST_MATLAB_OLD_BIN_DIR}/proxs CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_PROX_SRC_DIR ${FAUST_MATLAB_OLD_SRC_DIR}/proxs CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_UTILS_BIN_DIR ${FAUST_MATLAB_OLD_BIN_DIR}/utils CACHE INTERNAL "")
	set(FAUST_MATLAB_OLD_UTILS_SRC_DIR ${FAUST_MATLAB_OLD_SRC_DIR}/utils CACHE INTERNAL "")



	configure_file(${FAUST_MATLAB_OLD_SRC_DIR}/old_hierarchical_fact.m  ${FAUST_MATLAB_OLD_BIN_DIR}/old_hierarchical_fact.m)
	configure_file(${FAUST_MATLAB_OLD_SRC_DIR}/old_palm4MSA.m  ${FAUST_MATLAB_OLD_BIN_DIR}/old_palm4MSA.m)


	configure_file(${FAUST_MATLAB_OLD_SRC_DIR}/setup_old_FAUST.m ${FAUST_MATLAB_OLD_BIN_DIR}/setup_old_FAUST.m COPYONLY)

	#prox matlab function
	file(GLOB MATLAB_PROX RELATIVE ${FAUST_MATLAB_OLD_PROX_SRC_DIR}  "${FAUST_MATLAB_OLD_PROX_SRC_DIR}/*.m")
	#message(STATUS "MATLAB PROX : ${MATLAB_PROX}")
	foreach(prox ${MATLAB_PROX})
      		configure_file(${FAUST_MATLAB_OLD_PROX_SRC_DIR}/${prox} ${FAUST_MATLAB_OLD_PROX_BIN_DIR}/${prox} COPYONLY)	
   	endforeach()



	file(GLOB MATLAB_PROX RELATIVE ${FAUST_MATLAB_OLD_PROX_SRC_DIR}  "${FAUST_MATLAB_OLD_PROX_SRC_DIR}/*.m")
	#message(STATUS "MATLAB PROX : ${MATLAB_PROX}")
	foreach(prox ${MATLAB_PROX})
      		configure_file(${FAUST_MATLAB_OLD_PROX_SRC_DIR}/${prox} ${FAUST_MATLAB_OLD_PROX_BIN_DIR}/${prox} COPYONLY)	
   	endforeach()

	
	file(GLOB MATLAB_UTILS RELATIVE ${FAUST_MATLAB_OLD_UTILS_SRC_DIR}  "${FAUST_MATLAB_OLD_UTILS_SRC_DIR}/*.m")
	#message(STATUS "MATLAB PROX : ${MATLAB_PROX}")
	foreach(prox ${MATLAB_UTILS})
      		configure_file(${FAUST_MATLAB_OLD_UTILS_SRC_DIR}/${prox} ${FAUST_MATLAB_OLD_UTILS_BIN_DIR}/${prox} COPYONLY)	
   	endforeach()



	#### installation ####
	install(DIRECTORY ${FAUST_MATLAB_BIN_DIR}/ DESTINATION ${CMAKE_INSTALL_MATLAB_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION} DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN "*.${MEX_EXT}" PATTERN "*.m" PATTERN "*.mat" PATTERN "*CMakeFiles*" EXCLUDE PATTERN "*old_matlab*" EXCLUDE)
	

	






