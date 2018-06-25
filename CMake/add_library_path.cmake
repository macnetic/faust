##############################################################################
##                              Description:                                ##
##  	    							    	    ## 
##  Define 2 macro setting path :					    ## 
##   add_library_path set the path where the library file will be searched  ##
##   add_include_path set the path where the  include directory will be     ## 
##		  searched						    ##		 
##									    ##      
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

function(add_library_path LIB_TMP_PATH)
	if(${ARGC} LESS 2)
		message(FATAL_ERROR "ARGC=${ARGC} : Wrong number of arguments in add_library_path macro")
	endif()
	string ( REGEX REPLACE "\\\\" "/" ARGN2 "${ARGN}")
	#set(list_var "${ARGN}")
	#foreach(loop_var IN LISTS list_var)
	foreach(loop_var IN LISTS ARGN2)
		if(EXISTS "${loop_var}")
			if(EXISTS "${loop_var}/lib64")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib64 )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			elseif(EXISTS "${loop_var}/lib")
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}/lib )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			else()
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} ${loop_var}  )
				set(${LIB_TMP_PATH} ${${LIB_TMP_PATH}} PARENT_SCOPE )
			endif()
		else()
			#message(STATUS "${loop_var} directory has been ignored, as it does not exist")
		endif()
	endforeach()
endfunction()


function(add_include_path INC_TMP_PATH)
	if(${ARGC} LESS 2)
		message(FATAL_ERROR "ARGC=${ARGC} : Wrong number of arguments in add_include_path macro")
	endif()
	string ( REGEX REPLACE "\\\\" "/" ARGN2 "${ARGN}")
	#set(list_var "${ARGN}")
	#foreach(loop_var IN LISTS list_var)
	foreach(loop_var IN LISTS ARGN2)
		if(EXISTS "${loop_var}")
			if(EXISTS "${loop_var}/include")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} "${loop_var}/include")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			elseif(EXISTS "${loop_var}/inc")
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} "${loop_var}/inc" )
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			else()
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} ${loop_var} )
				set(${INC_TMP_PATH} ${${INC_TMP_PATH}} PARENT_SCOPE )
			endif()
		else()
			#message(STATUS "${loop_var} directory has been ignored, as it does not exist")
		endif()
	endforeach()
endfunction()





