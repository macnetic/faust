##############################################################################
##                              Description:                                ##
##  cmake script to find Python executable.            			            ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.gforge.inria.fr>                         ##
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


message(STATUS "------------------------------------------------")
message(STATUS "------------ Looking for Python PATH -----------")
message(STATUS "------------------------------------------------")


if(UNIX)


	message(STATUS "INFO- If you want to choose an other version of Python,")
	message(STATUS "INFO- please add an environment variable PYTHON_PATH. ")
	message(STATUS "INFO- Example : PYTHON_PATH=/usr/bin/python")
	if ($ENV{PYTHON_PATH}} MATCHES python)
		set(PYTHON_EXE $ENV{PYTHON_PATH})
		message(STATUS "PYTHON_EXE=$ENV{PYTHON_PATH} defined from environment variable")
	elseif (${PYTHON_PATH} MATCHES python)
		set(PYTHON_EXE ${PYTHON_PATH})
		message(STATUS "PYTHON_EXE=${PYTHON_EXE} defined from local input variable")
	else() # PYTHON_EXE_DIR
        #set(Python_ADDITIONAL_VERSIONS 3)
        find_package(PythonInterp 3)# REQUIRED)
        find_package(PythonLibs 3)# REQUIRED)

        if(NOT PYTHONINTERP_FOUND)
            find_package(PythonInterp 2)
            find_package(PythonInterp 2)
            if(NOT PYTHONINTERP_FOUND)
                message(FATAL_ERROR "python is not installed or at least not present in your PATH ; Please insert in the PATH environment.")
            endif()
        endif(PYTHONINTERP_FOUND)
        set(PYTHON_EXE ${PYTHON_EXECUTABLE})
	endif()
else(UNIX)
	message(FATAL_ERROR "Python wrapper Unsupported OS (only compatible with Unix System (Linux or Mac Os X)")
endif()

string(REGEX REPLACE "[a-zA-Z0-9_/:.]+/p(ython.?)" "c\\1" CYTHON_BIN_NAME "${PYTHON_EXE}")

message(STATUS "------------------------------------------------")
message(STATUS "PYTHON_EXE has been found : ${PYTHON_EXE}")
message(STATUS "Probed cython program name: ${CYTHON_BIN_NAME}")
message(STATUS "------------------------------------------------")


message(STATUS " ")
message(STATUS "------------------------------------------------------")
message(STATUS "--- Looking for Python module (cython,numpy,scipy) ---")
message(STATUS "------------------------------------------------------")
exec_program("${PYTHON_EXE} ${PROJECT_SOURCE_DIR}/CMake/check_python.py" OUTPUT_VARIABLE LIST_PYTHON_MODULE  RETURN_VALUE PYTHON_MODULE_MISSING)
message("${LIST_PYTHON_MODULE}")

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


if(UNIX)
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
        find_program(CYTHON_TMP ${CYTHON_BIN_NAME})
        if(${CYTHON_TMP} MATCHES ".*-NOTFOUND")
            message(FATAL_ERROR "Cython is not present in your PATH ; Please insert in the PATH environment.")
        else()
            message(STATUS "Cython path: ${CYTHON_TMP}")
        endif()
	endif()
    string(REGEX REPLACE "([a-zA-Z0-9_/:.]+)/${CYTHON_BIN_NAME}.?" "\\1" CYTHON_BIN_DIR "${CYTHON_TMP}")
else(UNIX)
	message(FATAL_ERROR "Cython wrapper Unsupported OS (only compatible with Unix System (Linux or Mac Os X)")
endif()

message(STATUS "CYTHON_BIN_DIR has been found : ${CYTHON_BIN_DIR}")
message(STATUS "------------------------------------------------")


