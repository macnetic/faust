
set(FAUST_PYTHON_LIB_DIR "'${FAUST_BIN_DIR}'")

if(${PYTHON_USES_MATIO_STATIC_LIBS})
	if(NOT MATIO_STATIC_LIB_PATH OR NOT Z_STATIC_LIB_PATH OR NOT HDF5_STATIC_LIB_PATH)
		message(FATAL_ERROR "When you set PYTHON_USES_MATIO_STATIC_LIBS to ON the variables MATIO_STATIC_LIB_PATH, Z_STATIC_LIB_PATH and HDF5_STATIC_LIB_PATH must be set to the full paths of static libs (.a)")
	endif()
	if(NOT EXISTS ${MATIO_STATIC_LIB_PATH})
		message(FATAL_ERROR "The filepath for matio doesn't exist: ${MATIO_STATIC_LIB_PATH}")
	elseif(NOT EXISTS ${Z_STATIC_LIB_PATH})
		message(FATAL_ERROR "The filepath for zlib doesn't exist: ${Z_STATIC_LIB_PATH}")
	elseif(NOT EXISTS ${HDF5_STATIC_LIB_PATH})
		message(FATAL_ERROR "The filepath for hdf5 lib doesn't exist: ${HDF5_STATIC_LIB_PATH}")
	endif()
	set(PYTHON_EXT_EXTRA_OBJECTS "[ '${MATIO_STATIC_LIB_PATH}', '${HDF5_STATIC_LIB_PATH}', '${Z_STATIC_LIB_PATH}']") # the order matters for compilation
else()
	set(PYTHON_EXT_EXTRA_OBJECTS "None")
	get_filename_component(MATIO_LIB_DIR ${MATIO_LIB_FILE} DIRECTORY)
	get_filename_component(MATIO_LIB_NAME ${MATIO_LIB_FILE} NAME_WE)
	string(REGEX REPLACE "^lib" "" MATIO_LIB_NAME ${MATIO_LIB_NAME})
	set(FAUST_PYTHON_LIBS ",'${MATIO_LIB_NAME}'")
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},'${MATIO_LIB_DIR}'")
endif()



set(FAUST_PYTHON_INCLUDE_DIR "'${EIGEN_INC_DIR}',
			'${FAUST_SRC_LINEAR_OPERATOR_DIR}',
			'${FAUST_PYTHON_SRC_DIR}',
			'${FAUST_PYTHON_SRC_SRC_DIR}',
			'${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}',
			'${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}',
			'${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}',
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

#move quickstart demo file
configure_file(${FAUST_PYTHON_SRC_DIR}/quickstart.py ${FAUST_PYTHON_BIN_DIR}/quickstart.py COPYONLY)

# configure the setup.py.in into setup.py (equivalent of Makefile for Python)
configure_file(${FAUST_PYTHON_SRC_DIR}/setup.py.in ${FAUST_PYTHON_BIN_DIR}/setup.py @ONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCorePy.pyx ${FAUST_PYTHON_BIN_DIR}/FaustCorePy.pyx COPYONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustCoreCy.pxd ${FAUST_PYTHON_BIN_DIR}/FaustCoreCy.pxd COPYONLY)
configure_file(${FAUST_PYTHON_SRC_DIR}/FaustPy.py ${FAUST_PYTHON_BIN_DIR}/FaustPy.py COPYONLY)

add_custom_target(${FAUST_PYTHON_TARGET} ALL DEPENDS ${FAUST_TARGET})

foreach(PYTHON_EXE IN LISTS PYTHON_EXES)
	add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
		COMMAND   ${PYTHON_EXE}
		ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
		COMMENT   "Creating PyFaust")
endforeach()

#install the python wrapper
install(DIRECTORY ${FAUST_PYTHON_BIN_DIR}/ DESTINATION ${CMAKE_INSTALL_PYTHON_PREFIX} FILE_PERMISSIONS  ${INSTALL_FILE_PERMISSION} DIRECTORY_PERMISSIONS ${INSTALL_DIR_PERMISSION} FILES_MATCHING PATTERN "*.${MEX_EXT}" PATTERN "*.py" PATTERN "*.so" PATTERN "build" EXCLUDE PATTERN "setup*" EXCLUDE PATTERN "CMakeFiles*" EXCLUDE)
