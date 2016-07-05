###### chek and find Openblas library
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)

if ( (OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR) )
	message(STATUS "OpenBlas library is installed here : ${OPENBLAS_LIB_FILE}")
	message(STATUS "OpenBlas include is installed here : ${OPENBLAS_INC_DIR}")
else ( (OPENBLAS_LIB_FILE) AND (OPENBLAS_INC_DIR) )
	if(UNIX)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		message(STATUS "------------------------------------------------")		
		message(STATUS "------------------------------------------------")
		message(STATUS "------------ OPENBLAS LIB INSTALLATION ---------")
		message(STATUS "------------------------------------------------")
		message(STATUS "------------------------------------------------")
		set(OPENBLAS_LIB_NAME "v0.2.18.tar.gz")
		exec_program("tar xzf ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${OPENBLAS_LIB_NAME} -C ${CMAKE_SOURCE_DIR}/externals/unix")
		if(EXISTS ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS)		
			exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
		endif(EXISTS)
		if(EXISTS ${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS)
			exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS")		
		endif(EXISTS)		
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS* ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
		exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS && make TARGET=NEHALEM && make install PREFIX='${CMAKE_SOURCE_DIR}/externals/unix/OpenBLAS'")

#exec_program(" ${CMAKE_SOURCE_DIR}/externals/unix/sdk_OpenBLAS")
	elseif(WIN32)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/v0.2.18.tar.bz -C ${CMAKE_SOURCE_DIR}/externals/win")
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/win/sdk_openBlas* ${CMAKE_SOURCE_DIR}/externals/win/OpenBLAS")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library OpenBlas")	
	endif(UNIX)

	add_include_path(INCLUDE_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS")
	add_library_path(LIBRARY_PATH_LIST_TMP_OPENBLAS "${PROJECT_SOURCE_DIR}/externals/unix/OpenBLAS")
	
	#message(STATUS "INCLUDE_PATH_LIST_TMP2=${INCLUDE_PATH_LIST_TMP2}")
	#message(STATUS "LIBRARY_PATH_LIST_TMP2=${LIBRARY_PATH_LIST_TMP2}")
	
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_OPENBLAS}) # CACHE PATH "List of include paths used as PATH parameter in find_path")
	set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_OPENBLAS}) # CACHE PATH "List of include paths used as PATH parameter in find_path")

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
#if(FAUST_USE_OPENBLAS)
#	check_external_libraries(openblas OPENBLAS_LIB_FILE 0)
#	check_external_includes("cblas.h" OPENBLAS_INC_DIR 0)
#endif(FAUST_USE_OPENBLAS)


#find_path(OPENBLAS_LIB_DIR ${OPENBLAS_LIB_FILE})
#find_path(EIGEN_LIB_DIR ${EIGEN_LIB_FILE})
# if(FAUST_USE_SINGLEPRECISION)
	# set(CXX_MEX_FLAGS "${CXX_MEX_FLAGS} -DFAUST_SINGLE")
	# if(UNIX)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# elseif(WIN32)
		# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DFAUST_SINGLE" CACHE STRING "compile flags" FORCE)
	# endif()
	# message(STATUS "**** SINGLE PRECISION USE *****")
# endif()
