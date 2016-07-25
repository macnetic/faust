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

###### test if executable matlab is in the path ######
#if (BUILD_MATLAB_MEX_FILES)
if(UNIX)
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
		#message(STATUS "MATLAB_DIR_TMP 3 = ${MATLAB_DIR_TMP}")
   	endif()
	#message(STATUS "MATLAB_DIR_TMP 4 = ${MATLAB_DIR_TMP}")

elseif(WIN32)
	#message(STATUS "MATLAB_DIR_TMP = ${MATLAB_DIR_TMP}")
	#message(STATUS "where /R \"C:\\Program Files\\MATLAB\" matlab.exe")
	#exec_program("where /R \"C:\\Program Files\\MATLAB\" matlab.exe" OUTPUT_VARIABLE MATLAB_DIR_TMP)
	#message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")
	set(MATLAB_DIR_TMP "C:\\Program Files\\MATLAB\\R2015b\\bin\\matlab.exe")
	#message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")	
else()
	message(WARNING "Unknown type of plateform for matlab")		
endif()


if( ${MATLAB_DIR_TMP} MATCHES "matlab")
	if(UNIX)
		# string(REGEX REPLACE "([a-zA-Z0-9_/:]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
		string(REGEX REPLACE "([a-zA-Z0-9_/:.]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}") # sous mac on a un point ds le path .app
	elseif(WIN32)
		string(REGEX REPLACE "([a-zA-Z0-9_\\:.]+)\\\\bin\\\\matlab.exe" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	else()
		message(WARNING "Unknown type of plateform for matlab")	
	endif()
	set(MATLAB_ROOT ${MATLAB_ROOT} CACHE PATH "Matlab root directory")
	
	#message(STATUS "MATLAB_DIR_TMP ${MATLAB_DIR_TMP}")
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
                	endif()		
		# METHODE 1 (without using matlab)
		else(APPLE)		
			exec_program("ls ${MATLAB_ROOT}/extern/lib | grep -i glnx" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
			if("${MEX_SUBDIR_LIB}" STREQUAL "glnxa64")
				set(MEX_EXT  "mexa64")
			elseif("${MEX_SUBDIR_LIB}" STREQUAL "glnx86")
				set(MEX_EXT  "mexa32")
			endif()
		endif(APPLE)	
	

		# METHODE 2 (using matlab)
		#exec_program("matlab -wait -nodesktop -nojvm -nodisplay -r \"fid=fopen('${MATLAB_ARCH_FILE}','w');fprintf(fid,'%s\\n%s\\n',computer('arch'),mexext);fclose(fid);exit\" > ${FAUST_TMP_BUILD_DIR}/matlab_output.log")
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 1" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		#exec_program("cat ${MATLAB_ARCH_FILE} | head -n 2 | tail -n 1" OUTPUT_VARIABLE MEX_EXT)

	elseif(WIN32)
		# METHODE 1 (without using matlab)
		execute_process(COMMAND ${PROJECT_SOURCE_DIR}/CMake/matlab_arch.bat "2" "${MATLAB_ROOT}" OUTPUT_VARIABLE MEX_SUBDIR_LIB)
		string(REGEX REPLACE "\n" "" MEX_SUBDIR_LIB ${MEX_SUBDIR_LIB})
		if("${MEX_SUBDIR_LIB}" STREQUAL "win64")
			set(MEX_EXT  "mexw64")
		elseif("${MEX_SUBDIR_LIB}" STREQUAL "win32")
			set(MEX_EXT  "mexw32")
		endif()	
		
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
##################################################################

