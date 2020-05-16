
# data is retrieved from remote website
if(NOT REMOTE_DATA_URL OR NOT REMOTE_DATA_FILE)
	#file(MAKE_DIRECTORY ${FAUST_MATLAB_TOOLS_INSTALL_DIR})
	# copy BSL data matrices and Laplacians
	file(GLOB BSL_DATA_FILE  RELATIVE ${FAUST_DATA_MAT_DIR}
		"${FAUST_DATA_MAT_DIR}/*MEG*.mat"
		"${FAUST_DATA_MAT_DIR}/Laplacian*.mat"
		"${FAUST_DATA_MAT_DIR}/results_BSL_user.mat")
	foreach(data_BSL ${BSL_DATA_FILE})
		configure_file(${FAUST_DATA_MAT_DIR}/${data_BSL} ${FAUST_MATFAUST_DEMO_DATA_BIN_DIR}/${data_BSL} COPYONLY)
	endforeach()
	# copy data matrix for quick_start.m
	configure_file(${FAUST_DATA_MAT_DIR}/faust_quick_start.mat ${FAUST_MATFAUST_DEMO_DATA_BIN_DIR}/faust_quick_start.mat COPYONLY)
endif()


# matlab tools directory
# copy also license files and readme
file(GLOB MATLAB_TOOLS RELATIVE ${FAUST_MATLAB_TOOLS_SRC_DIR}  "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.m"  "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.txt")
foreach(MATLAB_TOOL ${MATLAB_TOOLS})
	configure_file(${FAUST_MATLAB_TOOLS_SRC_DIR}/${MATLAB_TOOL} ${FAUST_MATLAB_TOOLS_BIN_DIR}/${MATLAB_TOOL} COPYONLY)
	#install(FILES ${FAUST_MATLAB_TOOL_BIN_DIR}/${MATLAB_TOOL}.m DESTINATION ${FAUST_MATLAB_TOOL_INSTALL_DIR} PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)
endforeach()

#FILE(COPY ${FAUST_MATLAB_NAMESPACE_SRC_DIR} DESTINATION ${PROJECT_BINARY_DIR}/tmp_gendoc_matfaust)


#message(FATAL_ERROR ${FAUST_MATLAB_MATFAUST_FILE_LIST})
foreach(MATFAUST_FILE IN LISTS FAUST_MATLAB_MATFAUST_FILE_LIST)
	string(REGEX REPLACE "[^+]*\\+matfaust\\/" "" MATFAUST_FILE ${MATFAUST_FILE})
	get_filename_component(PARENT_DIR ${MATFAUST_FILE} DIRECTORY)
	get_filename_component(MATFAUST_SHORT_FILE ${MATFAUST_FILE} NAME_WE )
	file(MAKE_DIRECTORY ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${PARENT_DIR})
	file(MAKE_DIRECTORY ${FAUST_MATLAB_BIN_DIR}/expfiltered)
	set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${FAUST_MATLAB_BIN_DIR}/expfiltered" "${FAUST_MATLAB_BIN_DIR}/+matfaust")
	if(EXPERIMENTAL_PKG)
		add_custom_command(
			OUTPUT ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMAND ${CMAKE_COMMAND} -E make_directory ${FAUST_MATLAB_BIN_DIR}/expfiltered/${PARENT_DIR}
			COMMAND
			${CMAKE_COMMAND} -E copy
			${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE}
			${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			DEPENDS ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE}
			COMMENT "Copying matfaust ${MATFAUST_FILE} code including experimental code.")
	else()
		add_custom_command(
			OUTPUT ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMAND ${CMAKE_COMMAND} -E echo "mkdir ${FAUST_MATLAB_BIN_DIR}/expfiltered/${PARENT_DIR}"
			COMMAND ${CMAKE_COMMAND} -E make_directory ${FAUST_MATLAB_BIN_DIR}/expfiltered/${PARENT_DIR}
			COMMAND python
			"${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/filter_experimental_code.py"
			${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE}
			${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} "--no-empty"
			DEPENDS ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE}
			COMMENT "Filtering matfaust ${MATFAUST_FILE} experimental code.")
	endif()
	add_custom_target(FILTERING_EXP_CODE_${MATFAUST_SHORT_FILE} ALL DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE})
	# use script generating inline doc only for real matlab classes
	if(${MATFAUST_FILE} MATCHES "@")
		message(STATUS "matfaust class: " ${MATFAUST_FILE})
		if(WIN32)
			set(SH_CMD git-bash)
		endif() # on Unix, count on bash that is in PATH
		add_custom_command(OUTPUT ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMENT "Generating matfaust inline doc for ${MATFAUST_FILE}")
		#execute_process(COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} ERROR_QUIET)
		#TODO: copy if error (fallback solution)
		#configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} COPYONLY)
		add_custom_target(GEN_MATFAUST_INLINE_DOC_${MATFAUST_SHORT_FILE} ALL DEPENDS ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE})
	else()
		message(STATUS "matfaust function file: " ${MATFAUST_FILE})
		# just copy the .m
		#TODO: inline doc should be possible for function files too
		#configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} @ONLY)
		add_custom_command(OUTPUT ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			COMMAND python
			ARGS ${PROJECT_SOURCE_DIR}/gen_doc/matlab_copy_function_file.py ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMENT "Generating matfaust inline doc for ${MATFAUST_FILE} (just copied because it's a function file).")
		add_custom_target(GEN_MATFAUST_INLINE_DOC_${MATFAUST_SHORT_FILE} ALL DEPENDS ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE})
	endif()
endforeach()

file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.cpp")
include_directories(${FAUST_MATLAB_TOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR})
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
	# Givens FGFT is only for double and float scalars
	if(NOT ${FAUST_SCALAR} MATCHES "complex")
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexfgftgivens.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexfgftgivens${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexsvdtj.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexsvdtj${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020${FSUFFIX}.m @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020${FSUFFIX}.m @ONLY)
	endif()
	# copy the *.m for factorization now, because we have the FSUFFIX in hands
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact${FSUFFIX}.m COPYONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA${FSUFFIX}.m COPYONLY)
	if(NOT ${FAUST_SCALAR} MATCHES "complex")
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexfgftgivens.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexfgftgivens${FSUFFIX}.m COPYONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexsvdtj.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexsvdtj${FSUFFIX}.m COPYONLY)
	endif()
endforeach()

file(GLOB CPP_MEX_FILES "${FAUST_MATLAB_MEX_SRC_DIR}/*.cpp")
## modif AL : bug pour les noms de PATH avec un point comme par exemple : ./faust_v1.10/interface_matlab/
## On va chercher les nom des target
#string(REGEX REPLACE "[a-zA-Z0-9_/:.]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")
string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")

if(WIN32)
	set(CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS_LIST} "-I${MATIO_INC_DIR}" "-I${MATLAB_INCLUDE_DIR}" "-I${EIGEN_INC_DIR}" "-I${FAUST_SRC_LINEAR_OPERATOR_DIR}" "-I${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}" "-I${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}" "-I${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}" "-I${FAUST_MATLAB_TOOLS_SRC_DIR}" "-I${FAUST_UTILS_SRC_DIR}")
	
	if(USE_GPU_MOD)
	       set(CMAKE_CXX_FLAGS_LIST ${CMAKE_CXX_FLAGS_LIST} -DUSE_GPU_MOD -I${GPU_MOD_INCLUDE_DIR} -I${GPU_MOD_INCLUDE_DIR}/../build/)
	endif()
	
	if(BUILD_OPENBLAS)
	       set(CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS_LIST}" "-I${OPENBLAS_INC_DIR}")
	endif()

	if(FAUST_TORCH)
		set(CMAKE_CXX_FLAGS_LIST "${CMAKE_CXX_FLAGS_LIST}" "-I${Torch_DIR}/../include" "-I${Torch_DIR}/../include/torch/csrc/api/include")
	endif()

	#windows version uses static lib of matio anyway (so replace any dyn. lib found by the static lib)
	string(REGEX REPLACE "^(.+)/.*$" "\\1/libmatio.lib" LDD_MEX_FLAGS "${MATIO_LIB_FILE}")
	set(LDD_MEX_FLAGS ${LDD_MEX_FLAGS} ${HDF5_STATIC_LIB_PATH} ${Z_STATIC_LIB_PATH})
	# complete static libs by objects of libs if needed (fallback solution from environment)
	set(LDD_MEX_FLAGS ${LDD_MEX_FLAGS} $ENV{WIN_HDF5_OBJS_DIR}/*.obj)
	set(LDD_MEX_FLAGS ${LDD_MEX_FLAGS} $ENV{WIN_SZIP_OBJS_DIR}/*.obj)
	file(GLOB CPP_SRC_FILES "${FAUST_SRC_LINEAR_OPERATOR_DIR}/*.cpp" "${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp" "${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*.cpp"  "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.cpp" "${FAUST_UTILS_SRC_DIR}/*.cpp")
	string(REGEX REPLACE "([a-zA-Z0-9_/:-]+)\\.cpp(;|$)" "\\1\\2" MEXSRC_FULL "${CPP_SRC_FILES}")

	string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "${FAUST_TMPMEX_DIR}/\\1.${MEXOBJ_EXT}\\2" MEXOBJ_FILES "${CPP_SRC_FILES}")

	set(MEX_OBJS_IN_ONE ${FAUST_TMPMEX_DIR}/mexobj_allinone.obj)
	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS MEXFILE_TARGET_LIST_TARGET )
	add_custom_target(MEXFILE_TARGET_LIST_TARGET DEPENDS ${MEXFILE_TARGET_LIST})

	# dirty trick to force the *.obj compilation before mex shared libs which depend on them
	add_custom_target(MEX_OBJS_IN_ONE
		COMMAND  ${MATLAB_ROOT}/bin/mex  -output ${MEX_OBJS_IN_ONE} "-v" "-c" "-largeArrayDims" ${CPP_SRC_FILES}
		${CMAKE_CXX_FLAGS_LIST} ${LDD_MEX_FLAGS}
		"-outdir" "${FAUST_TMPMEX_DIR}" "OBJEXT=.${MEXOBJ_EXT}"
		COMMENT "Compiling ${MEX_OBJS_IN_ONE} for mex usage")

	add_custom_target(MEXOBJ_FILES_TARGET DEPENDS MEX_OBJS_IN_ONE)
	add_dependencies(MEXFILE_TARGET_LIST_TARGET MEXOBJ_FILES_TARGET)
elseif(UNIX)
	find_package(Matlab)
	if(NOT Matlab_FOUND)
		message(FATAL_ERROR "Matlab not found.")
	endif()

	add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_TARGET_LIST})
endif()


foreach(mex_target ${MEXFILE_TARGET_LIST})
	if(UNIX)
		matlab_add_mex(NAME ${mex_target} SHARED SRC "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" R2017b)
		set_target_properties(${mex_target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${FAUST_MATLAB_MEX_BIN_DIR})
		target_link_options(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},LINK_OPTIONS>)
		#target_link_options(${mex_target} PRIVATE ${LDD_MEX_FLAGS})
		target_link_options(${mex_target} PRIVATE ${OpenMP_CXX_LIBRARIES}) # do nothing if BUILD_MULTITHREAD is OFF
		target_compile_options(${mex_target} PRIVATE "-v" "-largeArrayDims")
		target_include_directories(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},INTERFACE_INCLUDE_DIRECTORIES>)
		add_dependencies(${mex_target} ${FAUST_TARGET})
		add_dependencies(${mex_target} ${FAUST_MEXTOOLS_TARGET})
		set_property(TARGET ${mex_target} PROPERTY CXX_STANDARD 11) # TODO: variable in main CMakeLists.txt
		target_compile_definitions(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},COMPILE_DEFINITIONS>)
		target_compile_options(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},COMPILE_OPTIONS>)
		target_link_libraries(${mex_target} ${FAUST_TARGET})
		if(USE_MATIO_STATIC_LIBS)
			target_link_libraries(${mex_target} "${MATIO_STATIC_LIB_PATH}" "${HDF5_STATIC_LIB_PATH}" "${Z_STATIC_LIB_PATH}")
		else()
			target_link_libraries(${mex_target} ${MATIO_LIB_FILE})
		endif()
		#		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
		#			COMMAND   ${MATLAB_ROOT}/bin/mex
		#			ARGS "-v" "-largeArrayDims" "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" "-L${MATLAB_ROOT}/bin/${MEX_SUBDIR_LIB}" "-L${FAUST_INSTALL_LIB}" "-L${PROJECT_BINARY_DIR}"
		#			"LDFLAGS=${OMP_CXX_LFLAGS}"
		#			VERBATIM ${CMAKE_CXX_FLAGS_LIST}
		#			VERBATIM ${LDD_MEX_FLAGS}
		#			"-O"
		#			COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}")
	elseif(WIN32)
		add_custom_target(${mex_target} DEPENDS ${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp ${mex_target}.${MEX_EXT} MEXOBJ_FILES_TARGET )
		if(BUILD_MULTITHREAD)
			if (${CMAKE_GENERATOR} MATCHES "MinGW Makefiles" )
				set(OMP_COMPFLAG -fopenmp)
				set(OMP_LINKFLAG -fopenmp)
			else()
				# assume VS is used
				set(OMP_COMPFLAG /openmp)
				#set(OMP_LINKFLAG /openmp) # unrecognized by msvc cl (useless, compflag suffices)
			endif()
		endif()

		#message(STATUS "${MEXOBJ_FILES} ")
		file(GLOB OBJ_MEXTOOLS_FILES "${FAUST_TMPMEX_DIR}/*.${MEXOBJ_EXT}")
		add_custom_command(OUTPUT ${mex_target}.${MEX_EXT}
			#COMMAND set ARGS CL=/openmp /link /openmp #  doesn't work (https://docs.microsoft.com/en-us/cpp/build/reference/cl-environment-variables?view=vs-2019)
			#COMMAND echo ARGS CL=%CL%
			COMMAND ${MATLAB_ROOT}/bin/mex
			ARGS "-v" "-largeArrayDims" "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" ${OBJ_MEXTOOLS_FILES} ${MEXOBJ_FILES}
			${CMAKE_CXX_FLAGS_LIST} ${LDD_MEX_FLAGS} -DNOMINMAX COMPFLAGS="${OMP_COMPFLAG} $COMPFLAGS" LINKFLAGS="${OMP_LINKFLAG} $LINKFLAGS" # mex will replace -D by /D if VS is used for building
			# (this constant is to contervene the symbol conflicts between std::min/max and min/max macro from Microsoft)
			DEPENDS ${MEXOBJ_FILES}
			COMMENT   "Creating mex function ${mex_target}.${MEX_EXT}"
			)
		#message (STATUS "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp")
	else()
		message(FATAL_ERROR "Unable to compile mex functions for this platform")
	endif()
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/${mex_target}.m ${FAUST_MATLAB_MEX_BIN_DIR}/${mex_target}.m COPYONLY)
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
install(DIRECTORY ${FAUST_MATLAB_BIN_DIR}/ DESTINATION ${CMAKE_INSTALL_MATLAB_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION} DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN "*.${MEX_EXT}" PATTERN "*.m" PATTERN "*.mat" PATTERN "*CMakeFiles*" EXCLUDE PATTERN "*old_matlab*" EXCLUDE PATTERN BitReversalPermutation.m EXCLUDE
	PATTERN CooleyTukeyFact.m EXCLUDE
	PATTERN diagonalization_givens.m EXCLUDE
	PATTERN diagonalization_givens_parall.m EXCLUDE
	PATTERN dvp.m EXCLUDE
	PATTERN normalizeCol.m EXCLUDE
	PATTERN hadamard_mat.m EXCLUDE
	PATTERN generate_params.m EXCLUDE
	PATTERN gen_artificial_faust.m EXCLUDE
	PATTERN TransStockhamFact.m EXCLUDE
	PATTERN test_Givens.m EXCLUDE
	PATTERN test_GivensDiag.m EXCLUDE
	)#PATTERN faust_decompose.m EXCLUDE #TODO: uncomment later when
#  matfaust.demo.hadamard will be updated to use hierarchical)
