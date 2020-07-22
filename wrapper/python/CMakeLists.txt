if(NOT PYTHON_ENCODING)
	set(PYTHON_ENCODING utf-8)
endif()

set(FAUST_PYTHON_LIB_DIR "'${FAUST_BIN_DIR}'")

# set the variable for cython extension sources
set(PY_EXT_SOURCES "['${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx']")

set(FAUST_PYTHON_LIBS "'faust'")

if(${USE_MATIO_STATIC_LIBS})
	if(NOT WIN32)
		set(PYTHON_EXT_EXTRA_OBJECTS "[ '${MATIO_STATIC_LIB_PATH}', '${HDF5_STATIC_LIB_PATH}', '${Z_STATIC_LIB_PATH}']") # the order matters for compilation
	else() # WIN32 # intended for Visual Studio!
		# compile the pyx and all faust .cpp needed (instead to link to faust.lib)
		set(PY_EXT_SOURCES "${PY_EXT_SOURCES}+glob('${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}/*.cpp')+glob('${FAUST_SRC_LINEAR_OPERATOR_DIR}/*.cpp')+glob('${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}/*.cpp')+glob('${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}/*.cpp')+glob('${FAUST_UTILS_SRC_DIR}/*.cpp')")
		# extract the name (without .lib suffix) and the path of each libs:
		# 1/ matio
		get_filename_component(MATIO_LIB_DIR ${MATIO_STATIC_LIB_PATH} DIRECTORY)
		get_filename_component(MATIO_LIB_NAME ${MATIO_STATIC_LIB_PATH} NAME_WE)
		set(FAUST_PYTHON_LIBS "'${MATIO_LIB_NAME}'") # override variable (ignoring faust lib., recompiling from code)
		set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${MATIO_LIB_DIR}'")
		# 2/ zlib
		get_filename_component(Z_LIB_DIR ${Z_STATIC_LIB_PATH} DIRECTORY)
		get_filename_component(Z_LIB_NAME ${Z_STATIC_LIB_PATH} NAME_WE)
		set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${Z_LIB_NAME}'")
		set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${Z_LIB_DIR}'")
		# 3/ for HDF5 we get directly objects files (not the static lib because of linking errors) # temporary solution
		# the paths must be set in environment
		# SZIP is a dependency of HDF5
		set(PYTHON_EXT_EXTRA_OBJECTS "glob('$ENV{WIN_HDF5_OBJS_DIR}/*.obj')+glob('$ENV{WIN_SZIP_OBJS_DIR}/*.obj')")
	endif()
else() # in fact this is only for WIN32 MINGW
	set(PYTHON_EXT_EXTRA_OBJECTS "None")
	get_filename_component(MATIO_LIB_DIR ${MATIO_LIB_FILE} DIRECTORY)
	get_filename_component(MATIO_LIB_NAME ${MATIO_LIB_FILE} NAME_WE)
	string(REGEX REPLACE "^lib" "" MATIO_LIB_NAME ${MATIO_LIB_NAME})
	set(FAUST_PYTHON_LIBS "${FAUST_PYTHON_LIBS},'${MATIO_LIB_NAME}'")
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${MATIO_LIB_DIR}'")

endif()



set(FAUST_PYTHON_INCLUDE_DIR "'${EIGEN_INC_DIR}',
'${FAUST_SRC_LINEAR_OPERATOR_DIR}',
'${FAUST_PYTHON_SRC_DIR}',
'${FAUST_PYTHON_SRC_SRC_DIR}',
'${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}',
'${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}',
'${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}',
'${FAUST_UTILS_SRC_DIR}',
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


if(WIN32)
	set(FAUST_SETUP_PY_CFLAGS "'/w', '/EHsc', '/DNOMINMAX'")
	set(FAUST_SETUP_PY_LFLAGS "")
else() #UNIX
	set(FAUST_SETUP_PY_CFLAGS "'-std=c++11', '-O3'")
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
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DFAUST_TORCH'")
	# read the code to load torch libraries in python
	file(READ ${FAUST_PYTHON_SRC_DIR}/torch_loading.py TORCH_LIBS_LOADING)
	set(PYFAUST_PKG_SUFFIX "_torch")
endif()

if(USE_GPU_MOD)
	set(FAUST_PYTHON_INCLUDE_DIR "${FAUST_PYTHON_INCLUDE_DIR}, '${GPU_MOD_INCLUDE_DIR}', '${GPU_MOD_INCLUDE_DIR}/../build/'")
	set(FAUST_SETUP_PY_CFLAGS "${FAUST_SETUP_PY_CFLAGS}, '-DUSE_GPU_MOD'")

	set(GPU_MOD_LIB_PPATH ${FAUST_PYTHON_BIN_DIR}/pyfaust/lib)
	if(NOT BUILD_WRAPPER_MATLAB)
		# we are packaging a pip package or at least using python setup_tools (and not cmake)
		file(COPY ${GPU_MOD_INCLUDE_DIR}/../build/${GPU_MOD_LIB} DESTINATION ${GPU_MOD_LIB_PPATH})
	else()
		# MATLAB wrapper is compiled too, it means we share the libgm between wrappers and have to symlink it pyfaust
		file(MAKE_DIRECTORY ${GPU_MOD_LIB_PPATH})
		file(CREATE_LINK ../../../lib/${GPU_MOD_LIB} ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} RESULT SYMLINK_OK SYMBOLIC)
		# the relative path we'll be correct in a faust typical system package tree (relatively to CMAKE_INSTALL_PREFIX/python/pyfaust/lib/)
		message(STATUS SIMLINK ../../../lib/${GPU_MOD_LIB} ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB})
		install(FILES ${GPU_MOD_LIB_PPATH}/${GPU_MOD_LIB} DESTINATION ${CMAKE_INSTALL_PYTHON_PREFIX}/pyfaust/lib)
	endif()
endif()

# configure the setup.py.in into setup.py (equivalent of Makefile for Python)
configure_file(${FAUST_PYTHON_SRC_DIR}/setup.py.in ${FAUST_PYTHON_BIN_DIR}/setup.py @ONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/_FaustCorePy.pyx ${FAUST_PYTHON_BIN_DIR}/_FaustCorePy.pyx COPYONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCy.pxd ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd COPYONLY)


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
		"--no-empty")
	if(NOT EXISTS ${FILTERED_OUT_PY})
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

add_custom_target(${FAUST_PYTHON_TARGET} ALL DEPENDS ${FAUST_TARGET})

list(LENGTH PYTHON_EXES PYTHON_EXES_LEN)
list(LENGTH CYTHON_EXES CYTHON_EXES_LEN)
#TODO: check that the two LENGHTS are equal

#message(STATUS "PYTHON_EXES_LEN=${PYTHON_EXES_LEN}, PYTHON_EXES=${PYTHON_EXES}")
#message(STATUS "CYTHON_EXES_LEN=${CYTHON_EXES_LEN}, CYTHON_EXES=${CYTHON_EXES}")
if(UNIX)
	foreach(PYTHON_EXE IN LISTS PYTHON_EXES)
		exec_program("${PYTHON_EXE}" ARGS "--version" OUTPUT_VARIABLE PY_VER RETURN_VALUE PY_RES)
		#message(STATUS "PY_VER=${PY_VER}")
		string(REGEX REPLACE "Python ([0-9]\\.[0-9]).*$" "\\1" PY_VER_WITH_POINTS ${PY_VER})
		string(REGEX REPLACE "Python ([0-9])\\.([0-9]).*$" "\\1\\2" PY_VER ${PY_VER})
		if(PY_VER MATCHES 3)
			set(PY3_VER ${PY_VER_WITH_POINTS} PARENT_SCOPE)
			string(REGEX REPLACE ".*([0-9])$" "\\1" PY3_MINOR_VER ${PY_VER})
			set(PY3_MINOR_VER ${PY3_MINOR_VER} PARENT_SCOPE)
		endif()
		add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
			COMMAND   ${PYTHON_EXE}
			ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
			COMMENT   "Creating pyfaust${PY_VER}")
	endforeach()
elseif(WIN32)
	if(${CMAKE_GENERATOR} MATCHES "Visual Studio")
		foreach(PYTHON_EXE IN LISTS PYTHON_EXES)
			exec_program("${PYTHON_EXE}" ARGS "--version" OUTPUT_VARIABLE PY_VER RETURN_VALUE PY_RES)
			#message(STATUS "PY_VER=${PY_VER}")
			string(REGEX REPLACE "Python ([0-9]\\.[0-9]).*$" "\\1" PY_VER_WITH_POINTS ${PY_VER})
			string(REGEX REPLACE "Python ([0-9])\\.([0-9]).*$" "\\1\\2" PY_VER ${PY_VER})
			# ONLY Python 3 for Faust built with Visual Studio
			if(PY_VER MATCHES 3)
				set(PY3_VER ${PY_VER_WITH_POINTS} PARENT_SCOPE)
				add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
					COMMAND   ${PYTHON_EXE}
					ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
					COMMENT   "Creating pyfaust${PY_VER}")
				string(REGEX REPLACE ".*([0-9])$" "\\1" PY3_MINOR_VER ${PY_VER})
				set(PY3_MINOR_VER ${PY3_MINOR_VER} PARENT_SCOPE)
			endif()
		endforeach()
	endif()
else()
	message(FATAL_ERROR "Not handled system.")
endif()

#install the python wrapper
install(DIRECTORY ${FAUST_PYTHON_BIN_DIR}/ DESTINATION ${CMAKE_INSTALL_PYTHON_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION} DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN "pyfaust/*.py" PATTERN "*.so" PATTERN "*.pyd" PATTERN "build" EXCLUDE PATTERN "setup*" EXCLUDE PATTERN "CMakeFiles*" EXCLUDE)



