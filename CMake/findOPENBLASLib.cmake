##############################################################################
##                              Description:                                ##
##  cmake script to check and find/install OPENBLAS C++ library.            ##
##   (cf http://www.openblas.net/)					    ##	
##  2 output variable are assigned :					    ##
##    - OPENBLAS_INC_DIR which is the path 				    ##
##		to the include directory of OPENBLAS  			    ##
##    -OPENBLAS_LIB_FILE which is the path to library file of OPENBLAS	    ##
##         for instance path/libopenblas.so for linux, 			    ##
##	                path/libopenblas.dylib for Mac,			    ##
##			path/libopenblas.dll for Windows		    ##
##      								    ##	                  
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


set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)

if ( (OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR) )
	message(STATUS "OpenBlas library is installed here : ${OPENBLAS_LIB_FILE}")
	message(STATUS "OpenBlas include is installed here : ${OPENBLAS_INC_DIR}")
else ( (OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR) )
	if(UNIX)
		message(STATUS "------------------------------------------------")
		message(STATUS "------------ OPENBLAS LIB INSTALLATION ---------")
		message(STATUS "------------------------------------------------")
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")

		set(OPENBLAS_LIB_NAME "v0.2.18.tar.gz")
		exec_program("tar xzf ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${OPENBLAS_LIB_NAME} -C ${CMAKE_SOURCE_DIR}/externals/unix")
		exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
		exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS")		
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS* ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
		exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS && make --quiet TARGET=NEHALEM && make install PREFIX='${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS'")
		
		add_include_path(INCLUDE_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS")
		add_library_path(LIBRARY_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS")
	
#exec_program(" ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
	elseif(WIN32)
		#set(OPENBLAS_LIB_NAME "OpenBLAS-0.2.19")
		#message(STATUS "----------------------${CMAKE_SIZEOF_VOID_P}")
		if(CMAKE_SIZEOF_VOID_P MATCHES "4") #windows 32 bit 
			message(FATAL_ERROR "OpenBlas Library is available in externals only for windows 64 bit system. Please adjust OpenBlas version in directory externals/win/zipLibs, and set the corresponding name of library in the file findOPENBLASLib.cmake")
		else()
			set(OPENBLAS_LIB_NAME "OpenBLAS-v0.2.14-Win64-int64") 
		endif()

		# Download openblas http://www.openblas.net/
		# dezip in ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/${OPENBLAS_LIB_NAME}
		# https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio 
		# section 2. CMake and Visual Studio
		# do this from powershell so cmake can find visual studio
		# install perl on system (with exe): padre-on-strawberry https://code.google.com/archive/p/padre-perl-ide/downloads
		# make sure perl is in your path "perl -v".
		# cmake -G "Visual Studio 12 Win64" .
		# in visual studio : generate ALL_BUILD 

		message(STATUS "------------------------------------------------")		
		message(STATUS "------------ Looking for OPENBLAS LIB ---------")
		message(STATUS "------------------------------------------------")
		
		exec_program("${CMAKE_SOURCE_DIR}/externals/win/7z/x64/7za x ${PROJECT_SOURCE_DIR}/externals/win/zipLibs/${OPENBLAS_LIB_NAME}.zip -o${PROJECT_SOURCE_DIR}/externals/win -y")

		add_include_path(INCLUDE_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/win/${OPENBLAS_LIB_NAME}")
		add_library_path(LIBRARY_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/win/${OPENBLAS_LIB_NAME}")
		
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		#exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/v0.2.18.tar.bz -C ${CMAKE_SOURCE_DIR}/externals/win")
		#exec_program("mv ${CMAKE_SOURCE_DIR}/externals/win/sdk_openBlas* ${CMAKE_SOURCE_DIR}/externals/win/OpenBLAS")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library OpenBlas")	
	endif(UNIX)
	
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_OPENBLAS}) # CACHE PATH "List of include paths used as PATH parameter in find_path")
	set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_OPENBLAS}) # CACHE PATH "List of include paths used as PATH parameter in find_path")

	message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")
		
	check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)

if ( (OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR) )
	message(STATUS "OpenBlas library is installed here : ${OPENBLAS_LIB_FILE}")
	message(STATUS "OpenBlas include is installed here : ${OPENBLAS_INC_DIR}")
else((OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR))	
	message(STATUS "OpenBlas library is not available here : ${OPENBLAS_LIB_FILE}")
	message(STATUS "OpenBlas include is not available here : ${OPENBLAS_INC_DIR}")
	message(FATAL_ERROR "openBLAS lib is not installed on your system. Please check openBLAS install.")	
endif((OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR))

message(STATUS "------------------------------------------------")
message(STATUS "------------------------------------------------")	
################################################################
endif()


##################################################################
#if(BUILD_OPENBLAS)
#	check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
#	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
#endif(BUILD_OPENBLAS)


#find_path(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE})
#find_path(EIGEN_LIB_DIR ${EIGEN_LIB_FILE})
# if(BUILD_USE_SINGLEPRECISION)
	# set(CXX_MEX_FLAGS "${CXX_MEX_FLAGS} -DFAUST_SINGLE")
	# if(UNIX)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# elseif(WIN32)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# endif()
	# message(STATUS "**** SINGLE PRECISION USE *****")
# endif()
