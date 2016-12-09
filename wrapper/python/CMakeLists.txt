## CMAKELIST PRINCIPAL
#option(BUILD_WRAPPER_PYTHON "Enable building PyFaust (librairie) Python." ON)
#if(BUILD_WRAPPER_PYTHON)
#	set(FAUST_PYTHON_SRC_DIR ${FAUST_SRC_WRAPPER_DIR}/python CACHE INTERNAL "")
#	set(FAUST_PYTHON_SRC_SRC_DIR ${FAUST_PYTHON_SRC_DIR}/src/ CACHE INTERNAL "")
#	set(FAUST_PYTHON_BIN_DIR ${FAUST_BIN_WRAPPER_DIR}/python CACHE INTERNAL "")
#	set(FAUST_PYTHON_TARGET faust_python CACHE INTERNAL "")
#	add_subdirectory(${FAUST_PYTHON_SRC_DIR} ${FAUST_PYTHON_BIN_DIR})
	
#endif(BUILD_WRAPPER_PYTHON)




#message(STATUS "FFFFFFAUST_PYTHONNNNN ${EIGEN_INC_DIR}")
#message(STATUS "EIGEN DIR ${EIGEN_INC_DIR}")
#message(STATUS "WRAPPER PYTHON BIN DIR : ${FAUST_PYTHON_BIN_DIR}")
#message(STATUS "FAUST SRC DIR ${FAUST_SRC_LINEAR_OPERATOR_DIR}")
set(FAUST_PYTHON_INCLUDE_DIR"'${EIGEN_INC_DIR}',
			'${FAUST_SRC_LINEAR_OPERATOR_DIR}',
			'${FAUST_PYTHON_SRC_DIR}',
			'${FAUST_PYTHON_SRC_SRC_DIR}',
			'${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR}',
			'${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR}', 
			'${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR}'")


set(FAUST_PYTHON_LIB_DIR "'${FAUST_BIN_DIR}'")

if(BUILD_OPENBLAS)
	
	# in the setup.py (Python Makefile) , we need the directory where the library openblas) is stored
	# not the full path to the lib file
	get_filename_component(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE} PATH)

	set(FAUST_PYTHON_INCLUDE_DIR "${FAUST_PYTHON_INCLUDE_DIR},
			'${OPENBLAS_INC_DIR}'")
	set(FAUST_PYTHON_LIB_DIR "${FAUST_PYTHON_LIB_DIR},
			'${OPENBLAS_LIB_DIR}'")
	set(FAUST_PYTHON_LIB ",'openblas'")
endif(BUILD_OPENBLAS)




#configure_file(${FAUST_PYTHON_SRC_DIR}/setup.py.in ${FAUST_PYTHON_BIN_DIR}/setup.py @ONLY)
configure_file(${FAUST_PYTHON_SRC_DIR}/setup.py.in ${FAUST_PYTHON_BIN_DIR}/setup.py @ONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/PyFaust.pyx ${FAUST_PYTHON_BIN_DIR}/PyFaust.pyx COPYONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/CyFaust.pxd ${FAUST_PYTHON_BIN_DIR}/CyFaust.pxd COPYONLY)
configure_file(${FAUST_PYTHON_SRC_SRC_DIR}/FaustPy.py ${FAUST_PYTHON_BIN_DIR}/FaustPy.py COPYONLY)
#configure_file(${FAUST_PYTHON_SRC_DIR}/main.py ${FAUST_PYTHON_BIN_DIR}/main.py COPYONLY)

add_custom_target(${FAUST_PYTHON_TARGET} ALL DEPENDS ${FAUST_TARGET} )

add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
		        COMMAND   ${PYTHON_EXE} 
		        ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py"  "build_ext" "--inplace" 
		        COMMENT   "Creating PyFaust")


#add_custom_command(TARGET ${FAUST_PYTHON_TARGET}
#		        COMMAND   python 
#		        ARGS "${FAUST_PYTHON_BIN_DIR}/setup.py" "build_ext" "--build-lib" "${FAUST_PYTHON_BIN_DIR}" "--build-temp"  "${FAUST_PYTHON_BIN_DIR}"  
#		        COMMENT   "Creating PyFaust")
