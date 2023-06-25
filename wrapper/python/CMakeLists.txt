if(NOT PYTHON_ENCODING)
	set(PYTHON_ENCODING utf-8)
endif()

set(FAUST_PYTHON_LIB_DIR "'${FAUST_BIN_DIR}'")

if(${CMAKE_GENERATOR} MATCHES "Visual")
	# needed because depending on the configuration the faust.lib is put in subfolder Release, Debug, etc.
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR}, '${FAUST_BIN_DIR}/${CMAKE_BUILD_TYPE}'")
endif()

# set the variable for cython extension sources
set(PY_EXT_SOURCES "['${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx']")

set(FAUST_PYTHON_LIBS "'faust'")

# PYFAUST_PYTHON_REQUIREMENTS is defined in the project root CMakeLists.txt (because it is needed for win32 faust.nsi installer script)
string(REGEX REPLACE " " "', '" PYFAUST_PYTHON_REQUIREMENTS_AS_PYLIST ${PYFAUST_PYTHON_REQUIREMENTS})
set(PYFAUST_PYTHON_REQUIREMENTS_AS_PYLIST "'${PYFAUST_PYTHON_REQUIREMENTS_AS_PYLIST}'")

if(BUILD_READ_MAT_FILE AND NOT NO_MATIO)
	if(${USE_MATIO_STATIC_LIBS})
		# Build using matio static library and dependencies
		if(NOT WIN32)
			set(PYTHON_EXT_EXTRA_OBJECTS "[ '${MATIO_STATIC_LIB_PATH}', '${HDF5_STATIC_LIB_PATH}', '${Z_STATIC_LIB_PATH}',  '${SZIP_STATIC_LIB_PATH}']") # the order matters for compilation
		else() # WIN32 (tested on Visual Studio)
			# compile the pyx and all faust .cpp needed (instead to link to faust.lib)
			set(PY_EXT_SOURCES "${PY_EXT_SOURCES}+glob('${FAUST_PYTHON_SRC_SRC_DIR}/*.cpp')")
			# extract the name (without .lib suffix) and the path of each lib
			# 1/ matio
			get_filename_component(MATIO_LIB_DIR ${MATIO_STATIC_LIB_PATH} DIRECTORY)
			get_filename_component(MATIO_LIB_NAME ${MATIO_STATIC_LIB_PATH} NAME_WE)
			#set(FAUST_PYTHON_LIBS "'${MATIO_LIB_NAME}'") # override variable (ignoring faust lib., recompiling from code)
			set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${MATIO_LIB_NAME}'") # override variable (ignoring faust lib., recompiling from code)
			set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${MATIO_LIB_DIR}'")
			# 2/ zlib
			get_filename_component(Z_LIB_DIR ${Z_STATIC_LIB_PATH} DIRECTORY)
			get_filename_component(Z_LIB_NAME ${Z_STATIC_LIB_PATH} NAME_WE)
			set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${Z_LIB_NAME}'")
			set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${Z_LIB_DIR}'")
			# 3/ hdf5
			get_filename_component(HDF5_LIB_DIR ${HDF5_STATIC_LIB_PATH} DIRECTORY)
			get_filename_component(HDF5_LIB_NAME ${HDF5_STATIC_LIB_PATH} NAME_WE)
			set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${HDF5_LIB_NAME}'")
			set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${HDF5_LIB_DIR}'")
			set(PYTHON_EXT_EXTRA_OBJECTS "None")
		endif()
	else()
		# Build using matio dynamic library or static standalone library (including defs for zlib, hdf5 and szip -- which is normally in hdf5 static lib)
		set(PYTHON_EXT_EXTRA_OBJECTS "None")
		get_filename_component(MATIO_LIB_DIR ${MATIO_LIB_FILE} DIRECTORY)
		get_filename_component(MATIO_LIB_NAME ${MATIO_LIB_FILE} NAME_WE)
		if((NOT WIN32 OR ${CMAKE_GENERATOR} MATCHES "MinGW Makefiles") AND NOT ${MATIO_LIB_FILE} MATCHES ".*\\.a")
			string(REGEX REPLACE "^lib" "" MATIO_LIB_NAME ${MATIO_LIB_NAME})
		endif()
		set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${MATIO_LIB_NAME}'")
		set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${MATIO_LIB_DIR}'")
	endif()
else()

	set(PYTHON_EXT_EXTRA_OBJECTS "None")
endif()

set(PROC_LIST CPU)
if(USE_GPU_MOD)
	list(APPEND PROC_LIST GPU)
endif()


# remove old generated files (in case cmake config has changed -- BUILD_COMPLEX_PYX or USE_GPU_MOD)
foreach(PROC IN ITEMS CPU GPU)
	file(GLOB FC_CPP_FILES ${FAUST_PYTHON_BIN_DIR}/FaustCoreCpp${PROC}*.cpp)
	list(LENGTH FC_CPP_FILES LEN)
	if(${LEN} GREATER 0)
		file(REMOVE ${FC_CPP_FILES})
	endif()
endforeach()

# Generate cpp files from generic cpp.in
foreach(PROC IN LISTS PROC_LIST)
	set(TYPES double)
	if(BUILD_COMPLEX_PYX)
		list(APPEND TYPES complex<double>)
	endif()
	foreach(TYPE IN LISTS TYPES)
		string(REGEX REPLACE "::" "" TYPE_NAME ${TYPE})
		string(REGEX REPLACE "[<>]" "" TYPE_NAME ${TYPE_NAME})
		configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCpp${PROC}.cpp.in ${FAUST_PYTHON_BIN_DIR}/FaustCoreCpp${PROC}${TYPE_NAME}.cpp)
	endforeach()
endforeach()

# TODO: support complex float (when ok put this in the foreach loop above)
if(BUILD_FLOAT_PYX)
	set(TYPE_NAME float)
	set(PROC CPU)
	set(TYPE float)
	configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCpp${PROC}.cpp.in ${FAUST_PYTHON_BIN_DIR}/FaustCoreCpp${PROC}${TYPE_NAME}.cpp)
	if(USE_GPU_MOD)
		set(TYPE_NAME float)
		set(PROC GPU)
		set(TYPE float)
		configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCpp${PROC}.cpp.in ${FAUST_PYTHON_BIN_DIR}/FaustCoreCpp${PROC}${TYPE_NAME}.cpp)
	endif()
endif()


set(PY_EXT_SOURCES "${PY_EXT_SOURCES}+glob('${FAUST_PYTHON_BIN_DIR}/FaustCoreCppCPU*.cpp')+['${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCpp.cpp']")

if(USE_GPU_MOD)
	set(PY_EXT_SOURCES "${PY_EXT_SOURCES}+glob('${FAUST_PYTHON_BIN_DIR}/FaustCoreCppGPU*.cpp')")
endif()

set(FAUST_PYTHON_INCLUDE_DIR "'${EIGEN_INC_DIR}',
'${FAUST_SRC_LINEAR_OPERATOR_DIR}',
'${FAUST_PYTHON_SRC_DIR}',
'${FAUST_PYTHON_SRC_SRC_DIR}',
'${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}',
'${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}',
'${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}',
'${FAUST_UTILS_SRC_DIR}',
'${FAUST_UTILS_SRC_DIR}/matio',
'${MATIO_INC_DIR}'")


if(BUILD_OPENBLAS)
	# in the setup.py (Python Makefile) , we need the directory where the library openblas) is stored
	# not the full path to the lib file
	get_filename_component(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE} PATH)

	set(FAUST_PYTHON_INCLUDE_DIR "${FAUST_PYTHON_INCLUDE_DIR},
	'${OPENBLAS_INC_DIR}'")
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},
	'${OPENBLAS_LIB_DIR}'")
	set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'openblas'")
endif(BUILD_OPENBLAS)

# read the license header to insert in top of each pyfaust python file
file(READ ${FAUST_PYTHON_SRC_DIR}/PYFAUST_LICENSE_HEADER.txt PYFAUST_LICENSE_HEADER)
if(PYTORCH_EXP)
	file(READ ${FAUST_PYTHON_SRC_DIR}/pytorch_faust_exp.py PYTORCH_EXP_CODE)
endif()


if(WIN32)
	set(FAUST_SETUP_PY_CFLAGS "'/w', '/EHsc', '/DNOMINMAX'")
	set(FAUST_SETUP_PY_LFLAGS "")
else() #UNIX
	set(FAUST_SETUP_PY_CFLAGS "'-std=c++11', '-O3', '-fno-wrapv'") # without -fno-wrapv the butterfly factorization is twice slower
	set(FAUST_SETUP_PY_LFLAGS "'-s'")
endif()

if(BUILD_MULTITHREAD)
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '${OMP_CXX_FLAGS}'")
	# about msvc/cl.exe: "You may use either a forward slash (/) or a dash (-) to specify a compiler option."
	# (https://docs.microsoft.com/en-us/cpp/build/reference/compiler-options?view=vs-2015)
	# cf ./src/faust_linear_operator/CPU/faust_prod_opt.hpp
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-D_MUL_OMP_'")
	string(REGEX REPLACE "([^'])[ \t\r\n]([^'])" "\\1','\\2" FAUST_SETUP_PY_CFLAGS ${FAUST_SETUP_PY_CFLAGS})
	if(DEFINED OMP_CXX_LFLAGS)
		set(FAUST_SETUP_PY_LFLAGS "${FAUST_SETUP_PY_LFLAGS}, '${OMP_CXX_LFLAGS}'")
		string(REGEX REPLACE "([^'])[ \t\r\n]([^'])" "\\1','\\2" FAUST_SETUP_PY_LFLAGS ${FAUST_SETUP_PY_LFLAGS})
	endif()
endif()

if(FAUST_TORCH)
	get_filename_component(Torch_DIR ${TORCH_LIBRARY} DIRECTORY)
	set(FAUST_PYTHON_INCLUDE_DIR "${FAUST_PYTHON_INCLUDE_DIR}, '${Torch_DIR}/../include', '${Torch_DIR}/../include/torch/csrc/api/include'")
	set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS}, 'torch', 'c10'")
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR}, '${Torch_DIR}'")
	#variable Torch_DIR is initialized by Torch CMake scripts (or rather from main CMakeLists.txt because torch script needs two executions to set the Torch_DIR...)
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DFAUST_TORCH', '-D_GLIBCXX_USE_CXX11_ABI=1'")
	# read the code to load torch libraries in python
	file(READ ${FAUST_PYTHON_SRC_DIR}/torch_loading.py TORCH_LIBS_LOADING)
	set(PYFAUST_PKG_SUFFIX "_torch")
endif()


if(BUILD_MULTITHREAD AND (APPLE  OR ${CMAKE_CXX_COMPILER} MATCHES .*clang.* OR ${CMAKE_GENERATOR} MATCHES "Visual"))
	# insert the libomp py code loader
	file(READ ${FAUST_PYTHON_SRC_DIR}/libomp_loading.py OMP_LIB_LOADING)
endif()

if(USE_OPENBLASO AND UNIX)
	file(READ ${FAUST_PYTHON_SRC_DIR}/blasomp_loading.py BLASOMP_LOADING)
	set(FAUST_SETUP_PY_LFLAGS "${FAUST_SETUP_PY_LFLAGS}, '-lopenblaso'")
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DEIGEN_USE_BLAS'")
	set(PYFAUST_PKG_SUFFIX "-openblaso")
endif()

if(NO_MATIO)
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DNO_MATIO'")
endif()

if(USE_GPU_MOD)
	set(FAUST_PYTHON_INCLUDE_DIR "${FAUST_PYTHON_INCLUDE_DIR}, '${GPU_MOD_INCLUDE_DIR}', '${GPU_MOD_INCLUDE_DIR}/../build-cu11.4/', '${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR}'")
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DUSE_GPU_MOD'")

	set(GPU_MOD_LIB_PPATH ${FAUST_PYTHON_BIN_DIR}/pyfaust/lib)
	get_filename_component(GPU_MOD_EXT ${GPU_MOD_LIB} LAST_EXT)
	get_filename_component(GPU_MOD_NAME ${GPU_MOD_LIB} NAME_WLE)
	# TODO: used versions of CUDA and include dir should be CMake variables
	# do not support CUDA 9.2 on Windows (require too old version VS <= 2017)
	set(CUDA_VERSIONS 12.1;11.4)
	foreach(CUDA_CONF IN LISTS CUDA_VERSIONS)
		set(GPU_MOD_BUILD_DIR ${GPU_MOD_INCLUDE_DIR}/../build-cu${CUDA_CONF})
		if(NOT EXISTS ${GPU_MOD_BUILD_DIR})
			message(WARNING "the gpu_mod build directory ${GPU_MOD_BUILD_DIR} must exist and the library compiled in order to support CUDA ${CUDA_CONF} in the pyfaust wrapper.")
		endif()
		set(CUDA_CONF_GPU_MOD_NAME ${GPU_MOD_NAME}-cu${CUDA_CONF}${GPU_MOD_EXT})
		if(NOT BUILD_WRAPPER_MATLAB)
			# we are packaging a pip package or at least using python setup_tools (and not cmake)
			file(COPY ${GPU_MOD_BUILD_DIR}/${GPU_MOD_LIB} DESTINATION ${GPU_MOD_LIB_PPATH})
			message(STATUS COPY "${GPU_MOD_BUILD_DIR}/${GPU_MOD_LIB} DESTINATION ${GPU_MOD_LIB_PPATH}")
			file(RENAME ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} ${GPU_MOD_LIB_PPATH}/${CUDA_CONF_GPU_MOD_NAME})
			message(STATUS "RENAME ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} ${GPU_MOD_LIB_PPATH}/${CUDA_CONF_GPU_MOD_NAME}")
			# if the link already exists it will be overwritten
#			file(CREATE_LINK ${CUDA_CONF_GPU_MOD_NAME} ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} RESULT SYMLINK_OK SYMBOLIC)
#			message(STATUS "CREATE_LINK ${CUDA_CONF_GPU_MOD_NAME} ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} RESULT SYMLINK_OK SYMBOLIC")
			# in fact the symlink is useless because defaulty wheel packages don't support symlink (a wheel is a zip), it is converted to a copy resulting to a unnecessary larger whl file
		else()
			# MATLAB wrapper is compiled too, it means we share the libgm between wrappers and have to symlink it in pyfaust
			file(MAKE_DIRECTORY ${GPU_MOD_LIB_PPATH})
			file(CREATE_LINK ../../../lib/${CUDA_CONF_GPU_MOD_NAME} ${GPU_MOD_LIB_PPATH}/${CUDA_CONF_GPU_MOD_NAME} RESULT SYMLINK_OK SYMBOLIC)
			# the relative path we'll be correct in a faust typical system package tree (relatively to CMAKE_INSTALL_PREFIX/python/pyfaust/lib/)
			message(STATUS SIMLINK ../../../lib/${CUDA_CONF_GPU_MOD_NAME} ${GPU_MOD_LIB_PPATH}/${CUDA_CONF_GPU_MOD_NAME})
			install(FILES ${GPU_MOD_LIB_PPATH}/${CUDA_CONF_GPU_MOD_NAME} DESTINATION ${CMAKE_INSTALL_PYTHON_PREFIX}/pyfaust/lib)
			file(CREATE_LINK ${CUDA_CONF_GPU_MOD_NAME} ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} RESULT SYMLINK_OK SYMBOLIC)
			install(FILES ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} DESTINATION ${CMAKE_INSTALL_PYTHON_PREFIX}/pyfaust/lib)
		endif()
		# verify here that classifiers are available for the used CUDA versions: https://pypi.org/classifiers/
		string(REGEX REPLACE "\\..*" "" CUDA_MAJOR_VERSION ${CUDA_CONF})
		set(GPU_CUDA_VERSION_PYPI_CLASSIFIER "${GPU_CUDA_VERSION_PYPI_CLASSIFIER} 'Environment :: GPU :: NVIDIA CUDA :: ${CUDA_MAJOR_VERSION}', ")
	endforeach()
endif()

#if(BUILD_COMPLEX_PYX)
#	set(FAUST_CORE_CPLX_TYPE FaustCoreCplx) #TODO: rename variable to CORE_CLASS and class to *ComplexDouble
#else()
#	set(FAUST_CORE_CPLX_TYPE FaustCoreCplxDummy)
#endif()

# configure the setup.py.in into setup.py (equivalent of Makefile for Python)
configure_file(${FAUST_PYTHON_SRC_DIR}/setup.py.in ${FAUST_PYTHON_BIN_DIR}/setup.py @ONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/_FaustCorePy.pyx ${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx @ONLY)


############### generate pxd file
set(REAL_FPP double)
set(CPP_CORE_CLASS FaustCoreCppCPU)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCy.pxd ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd COPYONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreGenCy.pxd.in ${FAUST_PYTHON_BIN_DIR}/FaustCoreGenCyCPU.pxd @ONLY)
file(READ ${FAUST_PYTHON_BIN_DIR}/FaustCoreGenCyCPU.pxd CORE_CPU_PXD)
file(APPEND ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd ${CORE_CPU_PXD})
file(READ ${FAUST_PYTHON_SRC_SRC_DIR}/FaustAlgo.pxd ALGO_PXD)
file(APPEND ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd ${ALGO_PXD})

if(USE_GPU_MOD)
	set(CPP_CORE_CLASS FaustCoreCppGPU)
	configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreGenCy.pxd.in ${FAUST_PYTHON_BIN_DIR}/FaustCoreGenCyGPU.pxd @ONLY)
	file(READ ${FAUST_PYTHON_BIN_DIR}/FaustCoreGenCyGPU.pxd GPU_CORE_PXD)
	file(APPEND ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd ${GPU_CORE_PXD})
	file(READ ${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCyGPU.pxd GPU_CORE_PXD)
	file(APPEND ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd ${GPU_CORE_PXD})
endif()
########## end of pxd file generation

# This function generates pyx files from generic ones (pyx.in)
function (Generate_pyx TYPE_NAME TYPE REAL_TYPE_NAME REAL_TYPE PROC CORE_OBJ_SUFFIX FUNC_TYPE_SUFFIX CORE_OBJ_REAL_SUFFIX)
	# warning: TYPE is used in configure_file
	set(CORE_CLASS_CPU FaustCoreGen${TYPE_NAME}CPU)
	set(CORE_CLASS_GPU FaustCoreGen${TYPE_NAME}GPU)
	if(${PROC} MATCHES CPU)
		set(OTHER_PROC GPU)
		set(CORE_CLASS_OTHER_PROC ${CORE_CLASS_GPU})
	else()
		set(OTHER_PROC CPU)
		set(CORE_CLASS_OTHER_PROC ${CORE_CLASS_CPU})
	endif()
	string(TOLOWER ${OTHER_PROC} OTHER_PROC_LOWER)
	string(TOLOWER ${PROC} PROC_LOWER)

	set(CPP_CORE_CLASS FaustCoreCpp${PROC})
	set(CORE_CLASS FaustCoreGen${TYPE_NAME}${PROC})
	set(CORE_CLASS_REAL FaustCoreGen${REAL_TYPE_NAME}${PROC})
	set(CORE_OBJ core_faust_${CORE_OBJ_SUFFIX})
	set(CORE_OBJ_REAL core_faust_${CORE_OBJ_REAL_SUFFIX})
	set(PYX_LIST _FaustCoreGen)
	if(NOT DISABLE_ALGO_PYX)
		list(APPEND PYX_LIST _FaustAlgoGenProc)
		if(${PROC} MATCHES CPU)
			# these pyx are only on CPU
			list(APPEND PYX_LIST _FaustAlgoGen)
			if(BUILD_EIGTJ_SVDTJ)
				list(APPEND PYX_LIST _FaustAlgoGenGivens)
			endif()
		endif()
		if(USE_GPU_MOD)
			list(APPEND PYX_LIST _FaustCoreGenNonMemberFuncs)
		endif()
	endif()
	foreach(PYX_NAME IN LISTS PYX_LIST)
		configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/${PYX_NAME}.pyx.in ${FAUST_PYTHON_BIN_DIR}/${PYX_NAME}${TYPE_NAME}${PROC}.pyx @ONLY)
		file(READ ${FAUST_PYTHON_BIN_DIR}/${PYX_NAME}${TYPE_NAME}${PROC}.pyx PYX_VAR)
		file(APPEND ${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx ${PYX_VAR})
	endforeach()
endfunction()

foreach(PROC IN LISTS PROC_LIST)
	Generate_pyx(Dbl double Dbl double ${PROC} dbl "" dbl)
	if(BUILD_COMPLEX_PYX)
		Generate_pyx(CplxDbl complex Dbl double ${PROC} cplx _cplx dbl)
		# Generate and include the complex-only class
		configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/_FaustAlgoCplxDblGenProc.pyx.in ${FAUST_PYTHON_BIN_DIR}/_FaustAlgoCplxDblGen${PROC}.pyx @ONLY)
		file(READ ${FAUST_PYTHON_BIN_DIR}/_FaustAlgoCplxDblGen${PROC}.pyx PYX_VAR)
		file(APPEND ${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx ${PYX_VAR})
	endif()
endforeach()

# TODO: support complex float and GPU (when ok put this in the foreach loop above)
if(BUILD_FLOAT_PYX)
	Generate_pyx(Flt float Flt float CPU flt "" flt)
	if(USE_GPU_MOD)
		Generate_pyx(Flt float Flt float GPU flt "" flt)
	endif()
endif()

if(NOT REMOTE_DATA_URL OR NOT REMOTE_DATA_FILE)
	# install data files used by pyfaust.demo
	# they will be included only in pip/easyinstall packages (see install cmake function/rule in the end of file)
	FILE(GLOB DEMO_MAT_FILES RELATIVE ${FAUST_DATA_MAT_DIR} ${FAUST_DATA_MAT_DIR}/*MEG*.mat ${FAUST_DATA_MAT_DIR}/*quick_start*.mat ${FAUST_DATA_MAT_DIR}/Laplacian*.mat)
	foreach(MAT_FILE ${DEMO_MAT_FILES})
		message(STATUS "Copying demo mat files into ${FAUST_PYTHON_BIN_DIR}/data")
		configure_file(${FAUST_DATA_MAT_DIR}/${MAT_FILE} ${FAUST_PYTHON_BIN_DIR}/pyfaust/data/${MAT_FILE} COPYONLY)
	endforeach()
endif()

# embed data downloader script in pyfaust package
configure_file(${FAUST_MISC_DIR}/faust_data_downloader/faust_data_downloader.py ${FAUST_PYTHON_BIN_DIR}/pyfaust/datadl.py @ONLY)

# copy license files which setup.py will need to package
configure_file(${PROJECT_SOURCE_DIR}/license.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/license_matio.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license_matio.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/license_hdf5.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license_hdf5.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/license_zlib.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license_zlib.txt COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/license_scipy.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license_scipy.txt COPYONLY)
if(USE_OPENBLASO)
	configure_file(${PROJECT_SOURCE_DIR}/license_openblas.txt ${FAUST_PYTHON_BIN_DIR}/pyfaust/license_openblas.txt COPYONLY)
endif()


#TODO: externalize this function to use it for matlab filetering too (e.g. in ${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/filter_experimental_code.cmake)
#TODO: take care that PYTHON_EXES is undefined when BUILD_WRAPPER_PYTHON=OFF
function(filter_experimental_code PY_INPUT PY_OUTPUT)
	file(TO_NATIVE_PATH "${PROJECT_SOURCE_DIR}/misc/binpkg_postinst_scripts/filter_experimental_code.py" FILTER_SCRIPT)
	file(TO_NATIVE_PATH ${PY_INPUT} PY_TO_FILTER)
	file(TO_NATIVE_PATH ${PY_OUTPUT} FILTERED_OUT_PY)
	message(STATUS "Filtering out exp. code from ${PY_TO_FILTER} into ${FILTERED_OUT_PY}")
	if(NOT EXISTS ${FILTER_SCRIPT})
		message(FATAL_ERROR "filter script not found (${FILTER_SCRIPT}).")
	endif()
	configure_file(${PY_TO_FILTER} ${FILTERED_OUT_PY} @ONLY)
	list(GET PYTHON_EXES -1 PY_EXE)
	execute_process(COMMAND
		${PY_EXE}
		"${FILTER_SCRIPT}"
		"${FILTERED_OUT_PY}"
		"${FILTERED_OUT_PY}"
		"--no-empty" OUTPUT_VARIABLE FILTER_OUTPUT)
	if (NOT FILTER_OUTPUT MATCHES "deleting" AND NOT EXISTS ${FILTERED_OUT_PY}) # if the file is empty it is deleted by the python script
		message(FATAL_ERROR "Filtering failed.")
	endif()
endfunction(filter_experimental_code)

message(STATUS "Installing pyfaust.py as pyfaust/__init__.py")
if(EXPERIMENTAL_PKG)
	configure_file(${FAUST_PYTHON_SRC_DIR}/pyfaust.py ${FAUST_PYTHON_BIN_DIR}/pyfaust/__init__.py @ONLY)
else()
	# filter experimental code out of pyfaust modules
	filter_experimental_code("${FAUST_PYTHON_SRC_DIR}/pyfaust.py" "${FAUST_PYTHON_BIN_DIR}/pyfaust/__init__.py")
endif()

FILE(GLOB_RECURSE PYFAUST_FILES RELATIVE ${FAUST_PYTHON_SRC_DIR} ${FAUST_PYTHON_SRC_DIR}/pyfaust/*py)
foreach(PYFAUST_FILE ${PYFAUST_FILES})
	message(STATUS "Installing ${PYFAUST_FILE}")
	if(EXPERIMENTAL_PKG)
		configure_file(${FAUST_PYTHON_SRC_DIR}/${PYFAUST_FILE} ${FAUST_PYTHON_BIN_DIR}/${PYFAUST_FILE} @ONLY)
	else()
		message(STATUS "Filtering exp. code from ${FAUST_PYTHON_SRC_DIR}/${PYFAUST_FILE} into ${FAUST_PYTHON_BIN_DIR}/${PYFAUST_FILE}")
		filter_experimental_code("${FAUST_PYTHON_SRC_DIR}/${PYFAUST_FILE}" "${FAUST_PYTHON_BIN_DIR}/${PYFAUST_FILE}")
	endif()
endforeach()

set(PY_CPP_TARGET python_cpp)
file(GLOB PY_CPP_FILES "${FAUST_PYTHON_BIN_DIR}/F*cpp")
add_library(${PY_CPP_TARGET} STATIC ${PY_CPP_FILES})
target_include_directories(${PY_CPP_TARGET} PRIVATE ${FAUST_PYTHON_SRC_SRC_DIR}) # py include dir is set later
if(USE_GPU_MOD)
	target_include_directories(${PY_CPP_TARGET} PRIVATE ${GPU_MOD_INCLUDE_DIR} PRIVATE ${GPU_MOD_INCLUDE_DIR}/../build-cu11.4/ PRIVATE ${FAUST_SRC_LINEAR_OPERATOR_GPU2_DIR})
	target_compile_definitions(${PY_CPP_TARGET} PRIVATE -DUSE_GPU_MOD)
endif()
string(REGEX REPLACE "\\.cpp" ".o" PY_CPP_OBJS ${PY_EXT_SOURCES})
set(PYTHON_EXT_EXTRA_OBJECTS ${PYTHON_EXT_EXTRA_OBJECTS}"+${PY_CPP_OBJS}") # TODO: it should be the static library not the objects but who cares?
set(PY_EXT_SOURCES "['${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx']") # reset src without cpp objs # TODO: add an option to do one way or another

add_custom_target(${FAUST_PYTHON_TARGET} ALL DEPENDS ${FAUST_TARGET} ${PY_CPP_TARGET})

list(LENGTH PYTHON_EXES PYTHON_EXES_LEN)
list(LENGTH CYTHON_EXES CYTHON_EXES_LEN)
#TODO: check that the two LENGHTS are equal

#message(STATUS "PYTHON_EXES_LEN=${PYTHON_EXES_LEN}, PYTHON_EXES=${PYTHON_EXES}")
#message(STATUS "CYTHON_EXES_LEN=${CYTHON_EXES_LEN}, CYTHON_EXES=${CYTHON_EXES}")
if(UNIX)
	foreach(PYTHON_EXE IN LISTS PYTHON_EXES)
		exec_program("${PYTHON_EXE}" ARGS "--version" OUTPUT_VARIABLE PY_VER_ RETURN_VALUE PY_RES)
		#message(STATUS "PY_VER=${PY_VER}")
		string(REGEX REPLACE "Python ([0-9]\\.[0-9]+).*$" "\\1" PY_VER_WITH_POINTS ${PY_VER_})
		string(REGEX REPLACE "Python ([0-9])\\.([0-9]+).*$" "\\1\\2" PY_VER ${PY_VER_})
		if(PY_VER MATCHES 3)
			set(PY3_VER ${PY_VER_WITH_POINTS} PARENT_SCOPE)
			string(REGEX REPLACE "Python ([0-9])\\.([0-9]+).*$" "\\2" PY3_MINOR_VER ${PY_VER_})
			set(PY3_MINOR_VER ${PY3_MINOR_VER} PARENT_SCOPE)
		endif()
		add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
			DEPENDS faust ${FAUST_PYTHON_SRC_DIR}/_FaustCorePy.pyx.in
			COMMAND   ${PYTHON_EXE}
			ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
			COMMENT   "Creating pyfaust${PY_VER}")
		if(APPLE)
			set(LIBEXT dylib)
		else()
			set(LIBEXT so)
		endif()
		if(APPLE AND BUILD_MULTITHREAD OR ${CMAKE_CXX_COMPILER} MATCHES .*clang.*)
			add_custom_target(EMBED_LIBOMP_IN_PYFAUST ALL
				COMMAND cmake -P ${PROJECT_SOURCE_DIR}/CMake/embed_sys_libomp_into_pyfaust.txt DEPENDS ${FAUST_PYTHON_TARGET} COMMENT "Embedding the libomp library in pyfaust/lib directory")
		endif()
		if(USE_OPENBLASO)
			add_custom_target(EMBED_OPENBLASO_IN_PYFAUST ALL
				COMMAND cmake -P ${PROJECT_SOURCE_DIR}/CMake/embed_sys_openblaso_into_pyfaust.txt DEPENDS ${FAUST_PYTHON_TARGET} COMMENT "Embedding the libopenblaso library in pyfaust/lib directory")
		endif()

		exec_program("${PYTHON_EXE}" ARGS "-c '\r\nfrom sysconfig import get_paths;print\(get_paths\(\)[\"include\"]\)'" OUTPUT_VARIABLE PY_C_INC_DIR RETURN_VALUE PY_RES)
		include_directories(${PY_CPP_TARGET} ${PY_C_INC_DIR})
	endforeach()
elseif(WIN32)
	if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
		foreach(PYTHON_EXE IN LISTS PYTHON_EXES)
			exec_program("${PYTHON_EXE}" ARGS "--version" OUTPUT_VARIABLE PY_VER_ RETURN_VALUE PY_RES)
			#message(STATUS "PY_VER=${PY_VER}")
			string(REGEX REPLACE "Python ([0-9]\\.[0-9]+).*$" "\\1" PY_VER_WITH_POINTS ${PY_VER_})
			string(REGEX REPLACE "Python ([0-9])\\.([0-9]+).*$" "\\1\\2" PY_VER ${PY_VER_})
			# ONLY Python 3 for Faust builds with Visual Studio
			if(PY_VER MATCHES 3)
				set(PY3_VER ${PY_VER_WITH_POINTS} PARENT_SCOPE)
				add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
					COMMAND   ${PYTHON_EXE}
					ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
					COMMENT   "Creating pyfaust${PY_VER}")
				string(REGEX REPLACE "Python ([0-9])\\.([0-9]+).*$" "\\1\\2" PY3_MINOR_VER ${PY_VER_})
				set(PY3_MINOR_VER ${PY3_MINOR_VER} PARENT_SCOPE)
			endif()
			if(BUILD_MULTITHREAD)
				# embed openmp dll in pyfaust
				if(NOT VCOMPLIB_PATH)
					message(FATAL_ERROR "The CMake option VCOMPLIB_PATH must have been set to the vcomp*.dll lib to embed in pyfaust")
				endif()
				# we do not use embed_sys_libomp_into_pyfaust because it doesn't indicate full paths for all libraries on Windows with VS
				configure_file(${VCOMPLIB_PATH} ${FAUST_PYTHON_BIN_DIR}/pyfaust/lib COPYONLY)
			endif()

			exec_program("${PYTHON_EXE}" ARGS "-c '\r\nfrom sysconfig import get_paths;print\(get_paths\(\)[\"include\"]\)'" OUTPUT_VARIABLE PY_C_INC_DIR RETURN_VALUE PY_RES)
		include_directories(${PY_CPP_TARGET} ${PY_C_INC_DIR})
		endforeach()
	endif()
else()
	message(FATAL_ERROR "Not handled system.")
endif()

#install the python wrapper
install(DIRECTORY ${FAUST_PYTHON_BIN_DIR}/ DESTINATION
	${CMAKE_INSTALL_PYTHON_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION}
	DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN
	"pyfaust/*.py" PATTERN "pyfaust/tests/*.py" PATTERN "*.so" PATTERN "*.pyd"
	PATTERN "libomp.*"
	PATTERN "build" EXCLUDE PATTERN "dist" EXCLUDE PATTERN "setup*" EXCLUDE PATTERN "CMakeFiles*" EXCLUDE)



