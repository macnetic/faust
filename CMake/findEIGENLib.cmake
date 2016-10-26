##############################################################################
##                              Description:                                ##
##  cmake script to check and find/install Eigen C++ library.		    ##
##   (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)             ##
##  1 output variable are assigned :					    ##
##    - EIGEN_INC_DIR which is the path 				    ##
##		to the include directory of EIGEN  			    ##		    
##  No library file is required to use Eigen, only header file.     	    ##							    ##									    ##	                  
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

### find and/or install eigen lib (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of include paths used as PATH parameter in find_path")
check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0) # on test avec INCLUDE_PATH_LIST_TMP_DEFAULT
if (EIGEN_INC_DIR)
	message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
else (EIGEN_INC_DIR)
	message(STATUS "------------------------------------------------")		
	message(STATUS "------------------------------------------------")
	message(STATUS "------------- EIGEN LIB INSTALLATION -----------")
	message(STATUS "------------------------------------------------")
	message(STATUS "------------------------------------------------")
		
	if(UNIX)
		#exec_program("which hg | xargs echo" OUTPUT_VARIABLE MERCURIAL_HG_PATH_DIR)
	   	#exec_program("readlink ${MERCURIAL_HG_PATH_DIR}" OUTPUT_VARIABLE READLINK_TMP)
		#sudo dnf install eigen3-devel
		#exec_program("hg clone https://bitbucket.org/eigen/eigen/ ${CMAKE_SOURCE_DIR}/externals/eigen")
		
		# Eigen library come from following link. You can modifying svn version of eigen for latest stable release version.
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/tmp https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2")

		set(EIGEN_LIB_NAME "3.2.8.tar.bz2") 
		exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${EIGEN_LIB_NAME} -C ${CMAKE_SOURCE_DIR}/externals/unix")
		exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/eigen")
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/eigen-* ${CMAKE_SOURCE_DIR}/externals/unix/eigen")

		add_include_path(INCLUDE_PATH_LIST_TMP_EIGEN "${PROJECT_SOURCE_DIR}/externals/unix/eigen")
	elseif(WIN32)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win/eigen http://bitbucket.org/eigen/eigen/get/3.2.8.zip")
		set(EIGEN_LIB_NAME "eigen-eigen-07105f7124f9")
		#exec_program("rd ${CMAKE_SOURCE_DIR}/externals/win/eigen")
		set(7zBinDIR ${CMAKE_SOURCE_DIR}/externals/win/7z/x64/7za.exe)
		if ( 7zBinDIR )
			#message(STATUS "7za.exe is locate in your external directory")
			exec_program("${CMAKE_SOURCE_DIR}/externals/win/7z/x64/7za x ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/${EIGEN_LIB_NAME}.zip -o${CMAKE_SOURCE_DIR}/externals/win -y")
		else()
			message(STATUS "7z.exe must be present in your PATH.")
			exec_program("7z x ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/${EIGEN_LIB_NAME}.zip -o${CMAKE_SOURCE_DIR}/externals/win -y")
		endif()
		#exec_program("move ${CMAKE_SOURCE_DIR}/externals/win/${EIGEN_LIB_NAME} ${CMAKE_SOURCE_DIR}/externals/win/eigen")
		add_include_path(INCLUDE_PATH_LIST_TMP_EIGEN "${PROJECT_SOURCE_DIR}/externals/win/${EIGEN_LIB_NAME}")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library Eigen")	
	endif(UNIX)
	

	#add_include_path(INCLUDE_PATH_LIST_TMP2 "${PROJECT_SOURCE_DIR}/externals/eigen")
	#set(INCLUDE_PATH_LIST2 ${INCLUDE_PATH_LIST_TMP2} CACHE PATH "List of include paths used as PATH parameter in find_path")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_EIGEN})# CACHE PATH "List of include paths used as PATH parameter in find_path")
	#message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")
	check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)

if ( EIGEN_INC_DIR )
	message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
else(EIGEN_INC_DIR)	
	message(STATUS "Eigen library is not available here: ${EIGEN_INC_DIR}")
	message(FATAL_ERROR "Eigen library is not installed on your system. Please check Eigen install.")	
endif(EIGEN_INC_DIR)



	
message(STATUS "------------------------------------------------")
message(STATUS "------------------------------------------------")	
endif (EIGEN_INC_DIR)
#############################################################

