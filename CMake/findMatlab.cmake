# - Try to find a version of Matlab and headers/library required by the 
#   used compiler. It determines the right MEX-File extension and add 
#   a macro to help the build of MEX-functions.
#
# This module detects a Matlab's version between Matlab 7.0
# and Matlab 7.9
#
# This module defines: 
#  MATLAB_ROOT: Matlab installation path
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib
#  MATLAB_LIBRARIES:   required libraries: libmex, libmx, libeng
#  MEX_EXTENSION: MEX extension required for the current plateform
#  MATLAB_CREATE_MEX: macro to build a MEX-file
#
# The macro MATLAB_CREATE_MEX requires in this order:
#  - function's name which will be called in Matlab;
#  - C/C++ source files;
#  - third libraries required.

message(STATUS "------------------------------------------------")
message(STATUS "------------ Looking for Matlab PATH -----------")
message(STATUS "------------------------------------------------")
	
###### test if executable matlab is in the path ######
#if (BUILD_MATLAB_MEX_FILES)
if(UNIX)

	#SET(MATLAB_EXE_DIR " " CACHE STRING "force the directory of your expected matlab binary" FORCE )
	message(STATUS "INFO- If you want to choose an other version of Matlab,") 
	message(STATUS "INFO- please add an environment variable MATLAB_EXE_DIR. ")
	message(STATUS "INFO- Example : MATLAB_EXE_DIR=F:\\ProgramFile\\MATLAB\\R2002b\\bin\\matlab.exe ")
	if ($ENV{MATLAB_EXE_DIR}} MATCHES matlab)
		set(MATLAB_DIR_TMP $ENV{MATLAB_EXE_DIR})
		message(STATUS "MATLAB_DIR_TMP=$ENV{MATLAB_EXE_DIR} defined from environment variable")
	elseif (${MATLAB_EXE_DIR} MATCHES matlab)
		set(MATLAB_DIR_TMP ${MATLAB_EXE_DIR})
		message(STATUS "MATLAB_DIR_TMP=${MATLAB_EXE_DIR} defined from local input variable")
	else() # MATLAB_EXE_DIR

		#message(STATUS "MATLAB_DIR_TMP 1 = ${MATLAB_DIR_TMP}")
		exec_program("which matlab | xargs echo" OUTPUT_VARIABLE MATLAB_DIR_TMP)
		#message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")

		if (${MATLAB_DIR_TMP} MATCHES "which: no matlab in")
			set(MATLAB_DIR_TMP "")			
			message(FATAL_ERROR "matlab is not present in your PATH ; Please insert in the PATH environment.")
		endif()

	   	exec_program("readlink ${MATLAB_DIR_TMP}" OUTPUT_VARIABLE READLINK_TMP)
		#message(STATUS "READLINK_TMP = ${READLINK_TMP}")
		if(${READLINK_TMP} MATCHES matlab)
			set(MATLAB_DIR_TMP ${READLINK_TMP})
			message(STATUS "MATLAB_DIR_TMP = ${MATLAB_DIR_TMP}")
	   	endif()
		#message(STATUS "MATLAB_DIR_TMP 4 = ${MATLAB_DIR_TMP}")
	endif() # MATLAB_EXE_DIR
elseif(WIN32)
	message(STATUS "INFO- If you want to choose an other version of Matlab,") 
	message(STATUS "INFO- please add an environment variable MATLAB_EXE_DIR. ")
	message(STATUS "INFO- Example : MATLAB_EXE_DIR=F:\\ProgramFile\\MATLAB\\R2002b\\bin\\matlab.exe ")
	#SET(MATLAB_EXE_DIR " " CACHE STRING "force the directory of your expected matlab binary" FORCE )
	
	exec_program("${CMAKE_SOURCE_DIR}/CMake/find_matlab_path.bat" ${PROJECT_BINARY_DIR})
	FILE(READ "${PROJECT_BINARY_DIR}/tmp/tmpPathMatlab.txt" contents)
	#STRING(REGEX REPLACE "\n" "" contents "${contents}")

	# On enregistre la premiére ligne du fichier logPath.txt comme chemin de matlab.
	# Si plusieurs versions de matlab, WARNING
	string(REGEX REPLACE "(\n)[a-zA-Z0-9_/\\:.\n\ ]+" "\\1" contents1 "${contents}")
	string(REGEX REPLACE "\n" "" contents1 "${contents1}")

	# On garde la 2éme ligne
	#string(REGEX REPLACE "(\n)[a-zA-Z0-9_/\\:.]+(\n)" "\\1" contents2tmp "${contents}")
	#string(REGEX REPLACE "(\n)[a-zA-Z0-9_/\\:.\n]+" "\\1" contents2 "${contents2tmp}")
	#string(REGEX REPLACE "\n" "" contents2 "${contents2}")
	
	#message(STATUS "contents=${contents}")
	#message(STATUS "contents1=${contents1}")
	#message(STATUS "contents2tmp=${contents2tmp}")
	#message(STATUS "contents2=${contents2}")
	
	set(MATLAB_EXE_DIR_FILE_PATH "${contents1}")
	#set(MATLAB_EXE_DIR_TMP2 "${contents2}")
	#message(STATUS "MATLAB_EXE_DIR=${MATLAB_EXE_DIR_TMP}")
	#message(STATUS "MATLAB_EXE_DIR_TMP2=${MATLAB_EXE_DIR_TMP2}")
			
	if ($ENV{MATLAB_EXE_DIR}} MATCHES matlab)
		set(MATLAB_DIR_TMP $ENV{MATLAB_EXE_DIR})
		message(STATUS "MATLAB_DIR_TMP=$ENV{MATLAB_EXE_DIR} defined from environment variable")
	elseif (${MATLAB_EXE_DIR} MATCHES matlab)
		set(MATLAB_DIR_TMP ${MATLAB_EXE_DIR})
		message(STATUS "MATLAB_DIR_TMP=${MATLAB_EXE_DIR} defined from local input variable")
	elseif (${MATLAB_EXE_DIR_FILE_PATH} MATCHES matlab)
		set(MATLAB_DIR_TMP ${MATLAB_EXE_DIR_FILE_PATH})
	else()
		#message(STATUS "MATLAB_DIR_TMP = $ENV{MATLAB_EXE_DIR}")
		message(STATUS "MATLAB_EXE_DIR is not available. It corresponds to the path of matlab.exe }")
		message(FATAL_ERROR "Unknown path of matlab.exe. Please add environment variable MATLAB_EXE_DIR ")
	endif()
	
	#set(MATLAB_DIR_TMP "F:\\programFiles\\MATLAB\\R2014b\\bin\\matlab.exe")
	#message(STATUS "MATLAB_DIR_TMP 2222 = ${MATLAB_DIR_TMP}")	
else()
	message(WARNING "Unknown type of plateform for matlab")		
endif()


if( ${MATLAB_DIR_TMP} MATCHES "matlab")
	if(UNIX)
		#message(STATUS "MATLAB_DIR_TMP ${MATLAB_DIR_TMP}")
		#message(STATUS "MATLAB_ROOT TMP: ${MATLAB_ROOT}")
		# string(REGEX REPLACE "([a-zA-Z0-9_/:]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
		string(REGEX REPLACE "([a-zA-Z0-9_/:.]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}") # sous mac on a un point ds le path .app
		#message(STATUS "MATLAB_DIR_TMP2 ${MATLAB_DIR_TMP}")
		#message(STATUS "MATLAB_ROOT TMP2: ${MATLAB_ROOT}")		
		if( ${MATLAB_ROOT} MATCHES "matlab") # Dans ce cas, le remplacment ne s'est pas fait.. on essaie avec la double separation "//"
			string(REGEX REPLACE "([a-zA-Z0-9_/:.]+)/bin//matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}") # ds le cas ou on a /bin// matlab
			message(STATUS "MATLAB_DIR_TMP=${MATLAB_DIR_TMP}")
			message(STATUS "MATLAB_ROOT TMP=${MATLAB_ROOT}")
		endif()

	elseif(WIN32)
		string(REGEX REPLACE "([a-zA-Z0-9_\\:.]+)\\\\bin\\\\matlab.exe" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	else()
		message(WARNING "Unknown type of plateform for matlab")	
	endif()
	set(MATLAB_ROOT ${MATLAB_ROOT} CACHE PATH "Matlab root directory")
	
	
	#message(STATUS "MATLAB_ROOT has been found : ${MATLAB_ROOT}") # example : "/usr/local/MATLAB/R2014b"
	
	message(STATUS "MATLAB_ROOT has been found : ${MATLAB_ROOT}")

	set(MATLAB_INCLUDE_DIR "${MATLAB_ROOT}/extern/include" CACHE INTERNAL "Matlab include directory")
	set(MATLAB_ARCH_FILE "${FAUST_TMP_BUILD_DIR}/matlab_arch.txt")
	# LINUX AND APPLE METHOD ARE VERY SIMILAR, CODE COULD BE FACTORIZED
	if(UNIX) 
		if(APPLE)
			exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i mac" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
                	if("${MEX_SUBDIR_LIB}" STREQUAL "maci64")
                        	set(MEX_EXT  "mexmaci64")
                	elseif("${MEX_SUBDIR_LIB}" STREQUAL "maci")
                        	set(MEX_EXT  "mexmaci")
                	else()
						message(WARNING "No extension for mex function is available. (see ./CMAke/findMatlab.cmake)")
					endif()		
		# METHODE 1 (without using matlab)
		else(APPLE)		
			exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i glnx" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
			if("${MEX_SUBDIR_LIB}" STREQUAL "glnxa64")
				set(MEX_EXT  "mexa64")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "glnx86")
				set(MEX_EXT  "mexa32")
           	else()
				message(WARNING "No extension for mex function is available. (see ./CMAke/findMatlab.cmake)")						
			endif()
		endif(APPLE)	
	

		# METHODE 2 (using matlab)
		#exec_program("matlab -wait -nodesktop -nojvm -nodisplay -r \"fid=fopen('${MATLAB_ARCH_FILE}','w');fprintf(fid,'%s\\n%s\\n',computer('arch'),mexext);fclose(fid);exit\" > ${FAUST_TMP_BUILD_DIR}/matlab_output.log")
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 1" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 2 | tail -n 1" OUTPUT_VARIABLE MEX_EXT)

	elseif(WIN32)
		# METHODE 1 (without using matlab)
		#methode sure
		if ("$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "AMD64")
			set(MEX_EXT  "mexw64")
		elseif ("$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "x86")
			set(MEX_EXT  "mexw32")	
		else() 
			execute_process(COMMAND ${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat "2" "${MATLAB_ROOT}" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
			string(REGEX REPLACE "\n" "" MEX_SUBDIR_LIB ${MEX_SUBDIR_LIB})
			if ("${MEX_SUBDIR_LIB}" STREQUAL "win64")
				set(MEX_EXT  "mexw64")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "win32")
				set(MEX_EXT  "mexw32")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "win32win64")
				set(MEX_EXT  "mexw64")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "win64win32")
				set(MEX_EXT  "mexw64")
			else()
				message(WARNING "No extension for mex function is available. (see ./CMAke/findMatlab.cmake)")
			endif("${MEX_SUBDIR_LIB}" STREQUAL "win64")	
		endif("$ENV{PROCESSOR_ARCHITECTURE}" STREQUAL "AMD64")	
		#message(STATUS "TEST ALALALALALALALAL MATLAB_ROOT has been found : ${MEX_SUBDIR_LIB} ${MEX_EXT}")
		# METHODE 2 (using matlab)
		#exec_program("${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat 0 \"${MATLAB_ARCH_FILE}\"" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		#exec_program("${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat 1 \"${MATLAB_ARCH_FILE}\"" OUTPUT_VARIABLE MEX_EXT)
	else()
		message(WARNING "Unknown type of plateform for matlab")
	endif()
else()
	set(MATLAB_ROOT "" CACHE PATH "Matlab root directory")
	message(WARNING "Matlab executable seems not to be in the path. So \"MATLAB_ROOT\" and \"MATLAB_INCLUDE_DIR\" wont'be defined and mex files won't be compiled. if matlab is installed on your computer, please add matlabroot/bin tu the PATH and try again.")	
endif()
message(STATUS "------------------------------------------------")
##################################################################

