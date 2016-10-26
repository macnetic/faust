##############################################################################
##                              Description:                                ##
##  	    							    	    ## 
##  Define 2 macro finding libraries :					    ## 
##   check_external_libraries find the library file of the sepcified lib    ##
##   check_external_includes find the include directory of the specified lib##		 
##									    ##      
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

macro(check_external_libraries LIBRARY_NAME LIBRARY_FILE IS_NECESSARY)
	if(NOT ${ARGC} EQUAL 3)
		message(FATAL_ERROR "Wrong number of arguments in check_external_libraries macro")
	endif()
	find_library(${LIBRARY_FILE} ${LIBRARY_NAME} PATHS ${LIBRARY_PATH_LIST} NO_DEFAULT_PATH)
	if(${LIBRARY_FILE})
		message(STATUS "${LIBRARY_NAME} library found")
	else()
		#message(STATUS "${LIBRARY_FILE} has not been set")
		if(${IS_NECESSARY})
			message(FATAL_ERROR "${LIBRARY_NAME} library not found !!!")
		else()
			message(STATUS "${LIBRARY_NAME} library not found in the selected directory.")
		endif()
	endif()
endmacro()


macro(check_external_includes INCLUDE_NAME INCLUDE_FILE IS_NECESSARY)
	if(NOT ${ARGC} EQUAL 3)
		message(FATAL_ERROR "Wrong number of arguments in check_external_includes macro")
	endif()
	find_path(${INCLUDE_FILE} NAMES ${INCLUDE_NAME} PATHS ${INCLUDE_PATH_LIST} NO_DEFAULT_PATH)
	if(${INCLUDE_FILE})
		message(STATUS "${INCLUDE_NAME} header found")
	else()
		#message(STATUS "${INCLUDE_FILE} has not been set")
		if(${IS_NECESSARY})
			message(FATAL_ERROR "${INCLUDE_NAME} header not found !!!")
		else()
			message(STATUS "${INCLUDE_NAME} header not found in the selected directory.")
		endif()
	endif()
endmacro()


