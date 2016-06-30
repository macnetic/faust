###### find external libraries ######
if(WIN32)
	set(CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES} "lib" )
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} ".a")
endif(WIN32)


# functions tools is looking for library and include files in the environment variables and/or externals Path using the cmake find command.
include(CMake/add_library_path.cmake)
include(CMake/check_external_libraries.cmake)


#set(OPENBLAS_DEFAULT_DIR=/opt/OpenBLAS CACHE INTERNAL "")
#set(EIGEN_DEFAULT_DIR=/usr/include/eigen3 CACHE INTERNAL "")

#set(OPENBLAS_EXTERNALS_DIR=${PROJECT_SOURCE_DIR}/externals/OpenBLAS CACHE INTERNAL "")
#set(EIGEN_EXTERNALS_DIR=${PROJECT_SOURCE_DIR}/externals/eigen CACHE INTERNAL "") 

#message(STATUS "DEFAULT EIGEN_DEFAULT_DIR=${EIGEN_DEFAULT_DIR}")
#message(STATUS "DEFAULT EIGEN_EXTERNALS_DIR=${EIGEN_EXTERNALS_DIR}")
# /sw /opt/local
#add_library_path(LIBRARY_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{HDF5_ROOT_DIR}" "$ENV{OPENBLASDIR}" "/usr" "/usr/local" "/usr/local/lib" "/opt" "/opt/local" "/usr/lib/x86_64-linux-gnu/" "${PROJECT_SOURCE_DIR}/externals" )
#add_include_path(INCLUDE_PATH_LIST_TMP "$ENV{CUDADIR}" "$ENV{MATIODIR}" "$ENV{MATIODIRINC}" "$ENV{OPENBLASDIR}" "$ENV{EIGENDIR}" "/usr" "/usr/local" "/usr/include/libxml2" "/opt"  "/opt/local" "${PROJECT_SOURCE_DIR}/externals" )

#add_library_path(LIBRARY_PATH_LIST_TMP "${OPENBLAS_DEFAULT_DIR}" "${OPENBLAS_EXTERNALS_DIR}")
#add_include_path(INCLUDE_PATH_LIST_TMP "${OPENBLAS_DEFAULT_DIR}" "${OPENBLAS_EXTERNALS_DIR}" "${EIGEN_DEFAULT_DIR}" "${EIGEN_EXTERNALS_DIR}")


message(STATUS "******* Check externals library ***********")

# Default path library (where the library is automatically install)
add_library_path(LIBRARY_PATH_LIST_TMP 	"/opt/OpenBLAS"
										"${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS"
)

add_include_path(INCLUDE_PATH_LIST_TMP 	"/opt/OpenBLAS" 
										"${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS"
										"/usr/include/eigen3_"
										"${PROJECT_SOURCE_DIR}/externals/unix/eigen"
)


set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP} CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP} CACHE PATH "List of include paths used as PATH parameter in find_path")
#message(STATUS "DEFAULT LIBRARY_PATH_LIST=${LIBRARY_PATH_LIST}")
#message(STATUS "DEFAULT INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")

######## EIGEN ##############
include(CMake/findEIGENLib.cmake)
######## OPENBLAS ##############
if (FAUST_USE_OPENBLAS)
	include(CMake/findOPENBLASLib.cmake)
endif(FAUST_USE_OPENBLAS)

############################


###### test if executable matlab is in the path ######
if (FAUST_USE_MEX)
if(UNIX)
	#message(STATUS "MATLAB_DIR_TMP 1 = ${MATLAB_DIR_TMP}")
	exec_program("which matlab | xargs echo" OUTPUT_VARIABLE MATLAB_DIR_TMP)
	#message(STATUS "MATLAB_DIR_TMP 2 = ${MATLAB_DIR_TMP}")
   	exec_program("readlink ${MATLAB_DIR_TMP}" OUTPUT_VARIABLE READLINK_TMP)
	#message(STATUS "READLINK_TMP = ${READLINK_TMP}")
	if(${READLINK_TMP} MATCHES matlab)
		set(MATLAB_DIR_TMP ${READLINK_TMP})
		#message(STATUS "MATLAB_DIR_TMP 3 = ${MATLAB_DIR_TMP}")
   	endif()
	#message(STATUS "MATLAB_DIR_TMP 4 = ${MATLAB_DIR_TMP}")

elseif(WIN32)
	exec_program("where matlab.exe" OUTPUT_VARIABLE MATLAB_DIR_TMP)
else()
	message(WARNING "Unknown type of plateform for matlab")		
endif()
if( ${MATLAB_DIR_TMP} MATCHES "matlab")
	if(UNIX)
		string(REGEX REPLACE "([a-zA-Z0-9_/:]+)/bin/matlab" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	elseif(WIN32)
		string(REGEX REPLACE "([a-zA-Z0-9_\\:]+)\\\\bin\\\\matlab.exe" "\\1" MATLAB_ROOT "${MATLAB_DIR_TMP}")
	else()
		message(WARNING "Unknown type of plateform for matlab")	
	endif()
	set(MATLAB_ROOT ${MATLAB_ROOT} CACHE PATH "Matlab root directory")
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
endif(FAUST_USE_MEX)
##################################################################



check_external_libraries(matio MATIO_LIB_FILE 0)
check_external_libraries(xml2 XML2_LIB_FILE 0)
check_external_libraries(hdf5 HDF5_LIB_FILE 0)

#LDFLAGS = -L$(CUDADIR)/lib64 -L$(MATIODIR)/lib -lpthread -lm -lcublas -lcudart -lcusparse -lstdc++ -lgfortran -lz -lmatio -lhdf5
check_external_includes("libxml/parser.h" XML_INC_DIR 0)
check_external_includes("matio.h" MATIO_INC_DIR 0)



if (FAUST_USE_GPU)
	check_external_libraries(cublas CUBLAS_LIB_FILE 1)
	check_external_libraries(cudart CUDART_LIB_FILE 1)
	check_external_libraries(cusparse CUSPARSE_LIB_FILE 1)

	check_external_includes("cublas_v2.h" CUBLAS_V2_INC_DIR 1)
	check_external_includes("cusparse.h" CUSPARSE_INC_DIR 1)
	check_external_includes("cuda.h" CUDA_INC_DIR 1)
	check_external_includes("cuda_runtime_api.h" CUDA_RUNTIME_API_INC_DIR 1)
endif(FAUST_USE_GPU)


