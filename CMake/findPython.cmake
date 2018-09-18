##############################################################################
##                              Description:                                ##
##  cmake script to find Python executable.            			            ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##############################################################################





# This module defines: 
#  PYTHON_EXE: Python executable
#  CYTHON_EXE: Cython executable

function(check_py27 PYPATH)
	message(STATUS PYPATH=${PYPATH})
	if(EXISTS ${PYPATH})
		execute_process(COMMAND ${PYPATH} "--version" OUTPUT_VARIABLE PY2_VER ERROR_VARIABLE PY2_VER2 RESULT_VARIABLE RES)
#		message(STATUS PY2_VER=${PY2_VER} PY2_VER2=${PY2_VER2})
		if(${RES} EQUAL 0 AND ("${PY2_VER}" MATCHES "Python 2\\.7\\..*" OR "${PY2_VER2}" MATCHES "Python 2\\.7\\..*") )
			SET(PYTHON2_EXE ${PYPATH} PARENT_SCOPE)
		else()
			UNSET(PYTHON2_EXE)
		endif()
	endif()
endfunction(check_py27)

message(STATUS "------------------------------------------------")
message(STATUS "------------ Looking for Python PATH -----------")
message(STATUS "------------------------------------------------")


if(UNIX OR WIN32)


	message(STATUS "INFO- If you want to choose an other version of Python,")
	message(STATUS "INFO- please add an environment variable PYTHON_PATH or PYTHON2_PATH and PYTHON3_PATH for using the both major versions. ")
	message(STATUS "INFO- Example : PYTHON_PATH=/usr/bin/python")
	if("$ENV{PYTHON2_PATH}" MATCHES python OR "$ENV{PYTHON3_PATH}" MATCHES python)
		#message(STATUS $ENV{PYTHON2_PATH})
		set(PYTHON_EXES $ENV{PYTHON2_PATH};$ENV{PYTHON3_PATH})
		#message(STATUS "PYTHON_EXES=${PYTHON_EXES} set from environment.")
		set(PYTHON_EXECUTABLE $ENV{PYTHON2_PATH}) # necessary for doxypypy
	elseif ("$ENV{PYTHON_PATH}" MATCHES python)
		set(PYTHON_EXE $ENV{PYTHON_PATH})
		message(STATUS "PYTHON_EXE=$ENV{PYTHON_PATH} defined from environment variable")
		set(PYTHON_EXES ${PYTHON_EXE})
	elseif (${PYTHON_PATH} MATCHES python)
		set(PYTHON_EXE ${PYTHON_PATH})
		message(STATUS "PYTHON_EXE=${PYTHON_EXE} defined from local input variable")
		set(PYTHON_EXES ${PYTHON_EXE})
	else()  # auto-find pythons 2 and 3
		#set(Python_ADDITIONAL_VERSIONS 3)
		set(Python_ADDITIONAL_VERSIONS 3)
		find_package(PythonInterp 3 EXACT)# REQUIRED)
		find_package(PythonLibs 3)# REQUIRED)
		if(PYTHONINTERP_FOUND)
			set(PYTHON3_EXE ${PYTHON_EXECUTABLE})
		endif(PYTHONINTERP_FOUND)
		# search python2 anyway (even if we found python3)
		# cmake fails to find py2 if py3 was found
		# try to find python2 in same path as python3
		set(Python_ADDITIONAL_VERSIONS 2)
		find_package(PythonInterp 2 EXACT)# REQUIRED)
		find_package(PythonLibs 2 REQUIRED)
		if(PYTHONINTERP_FOUND AND NOT PYTHON3_EXE)
			set(PYTHON2_EXE ${PYTHON_EXECUTABLE})
		else()
			string(REGEX REPLACE "(.*)/python.*" "\\1/python2" TMP_PYTHON2_EXE ${PYTHON3_EXE})
			# message(STATUS "PYTHON3_EXE=${PYTHON3_EXE} TMP_PYTHON2_EXE=${TMP_PYTHON2_EXE}")
			check_py27(${TMP_PYTHON2_EXE})
			# if py27 still not found use find_program
			if(NOT PYTHON2_EXE)
				find_program(TMP_PYTHON2_EXE_ python2)
				check_py27(${TMP_PYTHON2_EXE_})
				message(STATUS "TMP_PYTHON2=${PYTHON2_EXE}")
			endif()
		endif()
		if(NOT PYTHON2_EXE)
			message(FATAL_ERROR "python is not installed or at least not present in your PATH ; Please insert in the PATH environment.")
		endif()
		set(PYTHON_EXES ${PYTHON2_EXE};${PYTHON3_EXE})
	endif()
else()
	message(FATAL_ERROR "Python wrapper Unsupported OS (only compatible with Unix System (Linux or Mac Os X) or Windows.")
endif()
#message(STATUS PYTHON_EXES=${PYTHON_EXES})
foreach(PYTHON_EXE IN LISTS PYTHON_EXES)

	string(REGEX REPLACE "[a-zA-Z0-9_/:.]+/p(ython.?)" "c\\1" CYTHON_BIN_NAME "${PYTHON_EXE}")
	string(REGEX REPLACE "(.*)(/|\\\\).*" "\\1" PYTHON_PARENT_PATH ${PYTHON_EXE}) # / separator for Unices, \ for Win
	if(NOT EXISTS ${PYTHON_PARENT_PATH}/${CYTHON_BIN_NAME})
		set(CYTHON_BIN_NAME cython)
	endif()

	message(STATUS "------------------------------------------------")
	message(STATUS "PYTHON_EXE has been found : ${PYTHON_EXE}")
	message(STATUS "Probed cython program name: ${CYTHON_BIN_NAME}")
	message(STATUS "------------------------------------------------")


	message(STATUS " ")
	message(STATUS "------------------------------------------------------")
	message(STATUS "--- Looking for Python module (cython,numpy,scipy) ---")
	message(STATUS "------------------------------------------------------")
	exec_program("${PYTHON_EXE} ${PROJECT_SOURCE_DIR}/CMake/check_python.py" OUTPUT_VARIABLE LIST_PYTHON_MODULE  RETURN_VALUE PYTHON_MODULE_MISSING)
	#message(STATUS "${LIST_PYTHON_MODULE}")

	set(PYTHON_MODULE_SCIPY ON)
	if(${PYTHON_MODULE_MISSING} EQUAL -1)
		message(FATAL_ERROR "necessary python module (numpy or cython) are missing !!!")
	elseif(${PYTHON_MODULE_MISSING} EQUAL 1)
		message(STATUS "optional python module scipy is missing, no time comparison with scipy will be made")
		set(PYTHON_MODULE_SCIPY OFF)
	else(${PYTHON_MODULE_MISSING})
		message(STATUS "all the python modules are installed (numpy,cython and scipy)")
	endif()

	message(STATUS "------------------------------------------------")
	message(STATUS " ")
	##################################################################
	message(STATUS "------------------------------------------------")
	message(STATUS "------------ Looking for Cython PATH -----------")
	message(STATUS "------------------------------------------------")


	message(STATUS "INFO- If you want to choose an other version of Cython,")
	message(STATUS "INFO- please add an environment variable CYTHON_PATH. ")
	message(STATUS "INFO- Example : CYTHON_PATH=/usr/bin/cython")
	if ($ENV{CYTHON_PATH}} MATCHES ${CYTHON_BIN_NAME})
		set(CYTHON_TMP $ENV{CYTHON_PATH})
		message(STATUS "CYTHON_PATH=$ENV{CYTHON_PATH} defined from environment variable")
	elseif (${CYTHON_PATH} MATCHES ${CYTHON_BIN_NAME})
		set(CYTHON_TMP ${CYTHON_PATH})
		message(STATUS "CYTHON_PATH=${CYTHON_PATH} defined from local input variable")
	else() # PYTHON_EXE_DIR
		set(CYTHON_TMP CYTHON_TMP-NOTFOUND) # forcing find_program to search again (don't rely on older value)
		if(UNIX)
			find_program(CYTHON_TMP ${CYTHON_BIN_NAME})
		elseif(WIN32)
			#message(STATUS "Python parent path ${PYTHON_PARENT_PATH}")
			find_program(CYTHON_TMP ${CYTHON_BIN_NAME} PATHS ${PYTHON_PARENT_PATH} PATH_SUFFIXES Scripts NO_DEFAULT_PATH)
			#message(STATUS "CYTHON_TMP=${CYTHON_TMP}")
		endif()
		if(${CYTHON_TMP} MATCHES ".*-NOTFOUND")
			message(FATAL_ERROR "Cython is not present in your PATH ; Please insert in the PATH environment.")
		else()
			message(STATUS "Cython path: ${CYTHON_TMP}")
		endif()
	endif()
	string(REGEX REPLACE "([a-zA-Z0-9_/:.]+)(\\|/)${CYTHON_BIN_NAME}.?" "\\1" CYTHON_BIN_DIR "${CYTHON_TMP}")

	message(STATUS "CYTHON_BIN_DIR has been found : ${CYTHON_BIN_DIR}")
	message(STATUS "------------------------------------------------")
	list(APPEND CYTHON_EXES ${CYTHON_BIN_DIR})
endforeach()
