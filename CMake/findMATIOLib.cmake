##############################################################################
##                              Description:                                ##
##  cmake script to check and find MATIO C++ library.                       ##
##  Library which allow to read and write 				    ##
##  files with extension ".mat" (used in Matlab)   			    ## 
##  2 output variable are assigned :					    ##
##    - MATIO_INC_DIR which is the path 				    ##
##		to the include directory of MATIO  			    ##
##    -MATIO_LIB_FILE which is the path to library file of MATIO	    ##
##         for instance path/libmatio.so for linux, 			    ##
##	                path/libmatio.dylib for Mac,			    ##
##			path/libmatio.dll for Windows		    	    ##
##      								    ##	                  
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


set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

check_external_includes("matio.h" MATIO_INC_DIR 0)
check_external_libraries(matio MATIO_LIB_FILE 0)

if ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
		message(STATUS "matio lib is here : ${MATIO_LIB_FILE}")
		message(STATUS "matio include is here : ${MATIO_INC_DIR}")
else ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
	if(UNIX)
		find_library(MATIO matio REQUIRED)
	elseif(WIN32)
		# TODO
	else(UNIX)
		message(WARNING "Unknown type of plateform for library OpenBlas")	
	endif(UNIX)

	set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_MATIO}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_MATIO})# CACHE PATH "List of include paths used as PATH parameter in find_path")
	check_external_includes("matio.h" MATIO_INC_DIR 0)
	check_external_libraries(matio MATIO_LIB_FILE 0)

	if ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
		message(STATUS "matio lib is here : ${MATIO_LIB_FILE}")
		message(STATUS "matio include is here : ${MATIO_INC_DIR}")
	else()
		message(STATUS "ERROR !!! matio is not installed !!!!!")
		message(FATAL_ERROR "matio lib is not installed on your system. Please check matio install.")	
	endif()

message(STATUS "------------------------------------------------")
message(STATUS "------------------------------------------------")	
################################################################
endif()



