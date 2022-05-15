
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
configure_file(${FAUST_MATLAB_SRC_DIR}/+matfaust/version.m.conf ${FAUST_MATLAB_SRC_DIR}/+matfaust/version.m @ONLY)
add_custom_target(PREPARE_MATFAUST_FOR_DOXYDOC ALL)
foreach(MATFAUST_FILE IN LISTS FAUST_MATLAB_MATFAUST_FILE_LIST)
	string(REGEX REPLACE "[^+]*\\+matfaust\\/" "" MATFAUST_FILE ${MATFAUST_FILE})
	get_filename_component(PARENT_DIR ${MATFAUST_FILE} DIRECTORY)
	get_filename_component(MATFAUST_SHORT_FILE ${MATFAUST_FILE} NAME_WE )
	set(MATFAUST_SHORT_FILE ${PARENT_DIR}_${MATFAUST_SHORT_FILE})
	string(REGEX REPLACE "[@/+]" "" MATFAUST_SHORT_FILE ${MATFAUST_SHORT_FILE})
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
	# use script generating inline doc only for real matlab classes and functions
	if(WIN32)
		set(SH_CMD git-bash)
	endif() # on Unix, count on bash that is in PATH
	if(${MATFAUST_FILE} MATCHES "@")
		message(STATUS "matfaust class: " ${MATFAUST_FILE})

		add_custom_command(OUTPUT ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMENT "Generating matfaust inline doc for ${MATFAUST_FILE}")
		#execute_process(COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} ERROR_QUIET)
		#TODO: copy if error (fallback solution)
		#configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} COPYONLY)
	else()
		message(STATUS "matfaust function file: " ${MATFAUST_FILE})
		# just copy the .m
		# generate inline doc for function files too
		#configure_file(${FAUST_MATLAB_NAMESPACE_SRC_DIR}/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} @ONLY)
#		add_custom_command(OUTPUT ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
#			COMMAND python
#			ARGS ${PROJECT_SOURCE_DIR}/gen_doc/matlab_copy_function_file.py ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
#			DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
#			COMMENT "Generating matfaust inline doc for ${MATFAUST_FILE} (just copied because it's a function file).")
		add_custom_command(OUTPUT ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE}
			COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/gen_matlab_inline_doc_from_doxy_blocks.sh ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE} ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} only_func
			DEPENDS ${FAUST_MATLAB_BIN_DIR}/expfiltered/${MATFAUST_FILE}
			COMMENT "Generating matfaust inline doc for ${MATFAUST_FILE}")

		add_custom_target(FILTER_OUT_MATLAB_DOXYBLOCKS_${MATFAUST_SHORT_FILE} ALL
			COMMAND ${PROJECT_SOURCE_DIR}/gen_doc/filterout_matlab_doxy_blocks.sh ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE} DEPENDS doc COMMENT "Removing doxygen blocks from function file.")
	endif()
	add_custom_target(GEN_MATFAUST_INLINE_DOC_${MATFAUST_SHORT_FILE} ALL DEPENDS ${FAUST_MATLAB_NAMESPACE_BIN_DIR}/${MATFAUST_FILE})
	add_dependencies(PREPARE_MATFAUST_FOR_DOXYDOC GEN_MATFAUST_INLINE_DOC_${MATFAUST_SHORT_FILE})
endforeach()

file(GLOB CPP_MEXTOOLS_FILES "${FAUST_MATLAB_TOOLS_SRC_DIR}/*.cpp")
include_directories(${FAUST_MATLAB_TOOLS_SRC_DIR} ${MATLAB_INCLUDE_DIR})
if(USE_GPU_MOD)
	include_directories(${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR} ${GPU_MOD_INCLUDE_DIR} ${GPU_MOD_INCLUDE_DIR}/../build-cu9.2)
endif()
add_library(${FAUST_MEXTOOLS_TARGET} OBJECT ${CPP_MEXTOOLS_FILES})
target_compile_definitions(${FAUST_MEXTOOLS_TARGET} PRIVATE MATLAB_DEFAULT_RELEASE=R2018a) #R2017b) # for complex (deprecated) separate complex API

set(MEX_TYPES_AND_SUFFIXES double:Real std::complex<double>:Cplx)
if(BUILD_FLOAT_MEX)
	list(APPEND MEX_TYPES_AND_SUFFIXES float:RealFloat)
endif()
foreach(SCALAR_AND_FSUFFIX IN LISTS MEX_TYPES_AND_SUFFIXES) # TODO: complex<float> versions
	# configuration de mexFaust.cpp.in
	# specialization of the class mexFaust template scalar with real or complex
	string(REGEX REPLACE "(.*):.*" "\\1" SCALAR ${SCALAR_AND_FSUFFIX})
	string(REGEX REPLACE ".*:(.*)" "\\1" FSUFFIX ${SCALAR_AND_FSUFFIX})
	message(STATUS "${SCALAR} ${FSUFFIX}")
	set(FAUST_SCALAR ${SCALAR})
	if (FSUFFIX MATCHES RealFloat)
		set(FAUST_FPP float)
		set(FACT_FPP float)
	else()
		set(FAUST_FPP double)
		set(FACT_FPP double)
	endif()
	set(DEV Cpu)
	set(DEV_IS_CPU 1)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexFaust.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexFaust${FSUFFIX}.cpp @ONLY)
	if(USE_GPU_MOD)
		set(DEV GPU2)
		unset(DEV_IS_CPU)
		#configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexFaustGPU.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexFaustGPU${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexFaust.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexFaustGPU${FSUFFIX}.cpp @ONLY)
	endif()
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexFaust.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexFaust${FSUFFIX}.m @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPalm4MSA.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPalm4MSA${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical_fact.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical_fact${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexButterfly.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexButterfly${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020${FSUFFIX}.m @ONLY)
	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020${FSUFFIX}.m @ONLY)
	if(USE_GPU_MOD)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020_gpu2.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020_gpu2${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020_gpu2.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPALM4MSA2020_gpu2${FSUFFIX}.m @ONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020_gpu2.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical2020_gpu2${FSUFFIX}.m @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020_gpu2.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020_gpu2${FSUFFIX}.cpp @ONLY)
	else()
		file(REMOVE ${FAUST_MATLAB_MEX_SRC_DIR}/mexPALM4MSA2020_gpu2${FSUFFIX}.cpp)  # delete a possible previous file generated with USE_GPU_MOD = ON (because otherwise it fails the building)
		file(REMOVE ${FAUST_MATLAB_MEX_SRC_DIR}/mexHierarchical2020_gpu2${FSUFFIX}.cpp)  # delete a possible previous file generated with USE_GPU_MOD = ON (because otherwise it fails the building)
	endif()

	# Givens FGFT is only for double and float scalars
	if(NOT ${FAUST_SCALAR} MATCHES "complex")
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexfgftgivens.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexfgftgivens${FSUFFIX}.cpp @ONLY)
		configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexsvdtj.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexsvdtj${FSUFFIX}.cpp @ONLY)
	endif()
	# copy the *.m for factorization now, because we have the FSUFFIX in hands
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexHierarchical_fact${FSUFFIX}.m COPYONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPalm4MSA${FSUFFIX}.m COPYONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexButterfly.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexButterfly${FSUFFIX}.m @ONLY)
	if(NOT ${FAUST_SCALAR} MATCHES "complex")
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexfgftgivens.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexfgftgivens${FSUFFIX}.m COPYONLY)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexsvdtj.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexsvdtj${FSUFFIX}.m COPYONLY)
	endif()

	configure_file(${FAUST_MATLAB_MEX_SRC_DIR}/mexPoly.cpp.in ${FAUST_MATLAB_MEX_SRC_DIR}/mexPoly${FSUFFIX}.cpp @ONLY)
	configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/mexPoly.m.in ${FAUST_MATLAB_DOC_SRC_DIR}/mexPoly${FSUFFIX}.m @ONLY)
endforeach()

file(GLOB CPP_MEX_FILES "${FAUST_MATLAB_MEX_SRC_DIR}/*.cpp")
## modif AL : bug pour les noms de PATH avec un point comme par exemple : ./faust_v1.10/interface_matlab/
## On va chercher les nom des target
#string(REGEX REPLACE "[a-zA-Z0-9_/:.]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")
string(REGEX REPLACE "[a-zA-Z0-9_/:.-]+/([a-zA-Z0-9_]+)\\.cpp(;|$)" "\\1\\2" MEXFILE_TARGET_LIST "${CPP_MEX_FILES}")

find_package(Matlab)
if(NOT Matlab_FOUND)
	message(FATAL_ERROR "Matlab not found.")
endif()

add_custom_target(${FAUST_MEX_TARGET} ALL DEPENDS ${MEXFILE_TARGET_LIST})

if(MATLAB_2017_API)
	set(MATLAB_API_VER "R2017b")
else()
	set(MATLAB_API_VER "R2018a")
endif()

#message(FATAL_ERROR "MEXFILE_TARGET_LIST=${MEXFILE_TARGET_LIST}")
foreach(mex_target ${MEXFILE_TARGET_LIST})
	matlab_add_mex(NAME ${mex_target} SHARED SRC "${FAUST_MATLAB_MEX_SRC_DIR}/${mex_target}.cpp" ${MATLAB_API_VER})
	set_target_properties(${mex_target} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${FAUST_MATLAB_MEX_BIN_DIR})
	target_link_options(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},LINK_OPTIONS>)
	target_compile_definitions(${mex_target} PRIVATE MATLAB_DEFAULT_RELEASE=${MATLAB_API_VER}) # R2017b for complex (deprecated) separate complex API
	if(MATLAB_2017_API)
		target_compile_definitions(${mex_target} PRIVATE __RELEASE_VERSION_DETECTOR__) # avoid mex.h to override the constants below
		target_compile_definitions(${mex_target} PRIVATE TARGET_API_VERSION=700) # set target matlab version for compiling to R2017b
		# to continue using mxGetImagData which is deprecated since R2018a # https://fr.mathworks.com/help/matlab/apiref/mxgetimagdata.html
		target_compile_definitions(${mex_target} PRIVATE MATLAB_TARGET_API_VERSION=700) # consistency expected with the constant above
	endif()
	if(WIN32)

		#target_compile_options(${mex_target} PRIVATE "${OMP_CXX_FLAGS}")
		# previous call is not necessary, to verify, you can check the build file generated for msbuild
		# for example wrapper/matlab/mex/mexFaustReal.vcxproj (it should contain a tag OpenMPSupport to True)
		target_compile_options(${mex_target} PRIVATE "/VERBOSE") # doesn't work at least with MSVC
		# complete static libs (linked later) by objects of libs if needed (fallback solution from environment)
		target_link_options(${mex_target} PRIVATE $ENV{WIN_HDF5_OBJS_DIR}/*.obj  $ENV{WIN_SZIP_OBJS_DIR}/*.obj)
		target_compile_definitions(${mex_target} PRIVATE NOMINMAX)
	else()
		target_link_options(${mex_target} PRIVATE ${OpenMP_CXX_LIBRARIES}) # do nothing if BUILD_MULTITHREAD is OFF
		target_compile_definitions(${mex_target} PRIVATE MATLAB_DEFAULT_RELEASE=${MATLAB_API_VER})
	endif()
	if(BUILD_MULTITHREAD)
		target_compile_definitions(${mex_target} PRIVATE _MUL_OMP_) # cf. ./src/faust_linear_operator/CPU/faust_prod_opt.hpp
		target_compile_options(${mex_target} PRIVATE "${OMP_CXX_FLAGS}")
	endif()
	target_include_directories(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},INTERFACE_INCLUDE_DIRECTORIES>)
	add_dependencies(${mex_target} ${FAUST_TARGET})
	add_dependencies(${mex_target} ${FAUST_MEXTOOLS_TARGET})
	set_property(TARGET ${mex_target} PROPERTY CXX_STANDARD 11) # TODO: variable in main CMakeLists.txt
	target_compile_definitions(${mex_target} PRIVATE $<TARGET_PROPERTY:${FAUST_TARGET},COMPILE_DEFINITIONS>)
	target_compile_options(${mex_target} PUBLIC $<TARGET_PROPERTY:${FAUST_TARGET},COMPILE_OPTIONS>)
	target_link_libraries(${mex_target} ${FAUST_TARGET})
	if(USE_MATIO_STATIC_LIBS OR WIN32)
		# Windows version must use static lib of matio anyway # TODO: FATAL_ERROR if MATIO_STATIC_LIB_PATH is undefined
		target_link_libraries(${mex_target} "${MATIO_STATIC_LIB_PATH}" "${HDF5_STATIC_LIB_PATH}" "${Z_STATIC_LIB_PATH}")
	else()
		target_link_libraries(${mex_target} ${MATIO_LIB_FILE})
	endif()
	if(WIN32)
		# not perfect patch to fix the issue of *.mexw64 located in Release subdir. with associated *.m only in the parent directory
		# that's because VS must be configured in Release mode for another reason (debug symbols missing in other libraries)
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/${mex_target}.m ${FAUST_MATLAB_MEX_BIN_DIR}/Release/${mex_target}.m COPYONLY)
	else()
		configure_file(${FAUST_MATLAB_DOC_SRC_DIR}/${mex_target}.m ${FAUST_MATLAB_MEX_BIN_DIR}/${mex_target}.m COPYONLY)
	endif()
endforeach()


# CONFIGURATION DES CHEMINS
configure_file(${FAUST_MATLAB_SRC_DIR}/setup_FAUST.m ${FAUST_MATLAB_BIN_DIR}/setup_FAUST.m COPYONLY)



#### installation ####
install(DIRECTORY ${FAUST_MATLAB_BIN_DIR}/ DESTINATION ${CMAKE_INSTALL_MATLAB_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION} DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN "*.${MEX_EXT}" PATTERN "*.m" PATTERN "*.mat" PATTERN "*CMakeFiles*" EXCLUDE PATTERN "*old_matlab*" EXCLUDE PATTERN "*expfiltered*" EXCLUDE PATTERN BitReversalPermutation.m EXCLUDE
	PATTERN CooleyTukeyFact.m EXCLUDE
	PATTERN diagonalization_givens.m EXCLUDE
	PATTERN diagonalization_givens_parall.m EXCLUDE
    #	PATTERN dvp.m EXCLUDE # needed for the BSL experiment
    #	PATTERN normalizeCol.m EXCLUDE # needed for the BSL experiment
	PATTERN hadamard_mat.m EXCLUDE
	PATTERN generate_params.m EXCLUDE
	PATTERN gen_artificial_faust.m EXCLUDE
	PATTERN TransStockhamFact.m EXCLUDE
	PATTERN test_Givens.m EXCLUDE
	PATTERN test_GivensDiag.m
	EXCLUDE PATTERN *m.bak
	)#PATTERN faust_decompose.m EXCLUDE #TODO: uncomment later when
#  matfaust.demo.hadamard will be updated to use hierarchical)
