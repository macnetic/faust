##############################################################################
##                              Description:                                ##
##  cmake script to find Python executable.            			    ##
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

		#message(STATUS "PYTHON_DIR_TMP 1 = ${PYTHON_DIR_TMP}")
		exec_program("which python | xargs echo" OUTPUT_VARIABLE PYTHON_EXE)
		#message(STATUS "PYTHON_DIR_TMP 2 = ${PYTHON_DIR_TMP}")
		
		if (${PYTHON_EXE} MATCHES "which: no python in")			
			message(FATAL_ERROR "python is not present in your PATH ; Please insert in the PATH environment.")
		endif()

	   	
		
	endif() 
else(UNIX)
	message(FATAL_ERROR "Python wrapper Unsupported OS (only compatible with Unix System (Linux or Mac Os X)")		
endif()


message(STATUS "PYTHON EXE ${PYTHON_EXE}")	
message(STATUS "------------------------------------------------")
##################################################################

