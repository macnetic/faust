##############################################################################
##                              Description:                                ##
##  cmake script to check and find/install Eigen C++ library.		    ##
##   (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)             ##
##  1 output variable are assigned :					    ##
##    - EIGEN_INC_DIR which is the path 				    ##
##		to the include directory of EIGEN  			    ##		    
##  No library file is required to use Eigen, only header file.     	    ##							    ##									    ##	                  
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
if(UNIX)
	find_package(Eigen3 3.4 REQUIRED)
	if(EIGEN3_INCLUDE_DIR)
		set(EIGEN_INC_DIR ${EIGEN3_INCLUDE_DIR})
		message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
		check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)
	else()
		message(FATAL_ERROR "Eigen not found, won't build.")
	endif()
elseif(WIN32)
	### find and/or install eigen lib (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)
	# TODO: set a cmake/env variable to vary the windows install path
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT} "C:/Program Files (x86)/Eigen3/include/eigen3") # CACHE PATH "List of include paths used as PATH parameter in find_path")
	check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0) # on test avec INCLUDE_PATH_LIST_TMP_DEFAULT
	if (EIGEN_INC_DIR)
		message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
	else (EIGEN_INC_DIR)
		message(FATAL_ERROR "Didn't find eigen installed on your system, it should be in C:\Program Files (x86)\Eigen3")
	endif (EIGEN_INC_DIR)
else()
	message(WARNING "Unknown type of plateform for library Eigen")
endif()
