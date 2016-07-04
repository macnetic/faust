###### chek and find Openblas library
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

check_external_includes("matio.h" MATIO_INC_DIR 0)
check_external_libraries(matio MATIO_LIB_FILE 0)

if ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
		message(STATUS "matio lib is here : ${MATIO_LIB_FILE}")
		message(STATUS "matio include is here : ${MATIO_INC_DIR}")
else ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
	if(UNIX)

		message(STATUS "------------------------------------------------")
		message(STATUS "------------------------------------------------")
		message(STATUS "------------ MATIO LIB INSTALLATION ------------")
		message(STATUS "------------------------------------------------")
		message(STATUS "------------------------------------------------")	
	
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		set(MATIO_LIB_NAME "matio-1.5.7.7z")
		exec_program("7z x ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${MATIO_LIB_NAME} -o${CMAKE_SOURCE_DIR}/externals/unix")
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/matio-* ${CMAKE_SOURCE_DIR}/externals/unix/matio")
		exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/matio && chmod -R 777 ./ && ./configure && make") # && make check
		# NOTE WARNING : WE DON'T run make install because it is not a root user install.  We used directly the lib and include in sdk_matio package. 
		# NOTE IF you have not 7z -->mac :  brew install p7zip
		add_include_path(INCLUDE_PATH_LIST_TMP_MATIO "${PROJECT_SOURCE_DIR}/externals/unix/matio/src")
		add_library_path(LIBRARY_PATH_LIST_TMP_MATIO "${PROJECT_SOURCE_DIR}/externals/unix/matio/src/.libs")	
	elseif(WIN32)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		#exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/v0.2.18.tar.bz -C ${CMAKE_SOURCE_DIR}/externals/win")
		#exec_program("mv ${CMAKE_SOURCE_DIR}/externals/win/sdk_openBlas* ${CMAKE_SOURCE_DIR}/externals/win/OpenBLAS")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library OpenBlas")	
	endif(UNIX)

	set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_MATIO}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_MATIO})# CACHE PATH "List of include paths used as PATH parameter in find_path")
	check_external_includes("matio.h" MATIO_INC_DIR 0)
	check_external_libraries(matio MATIO_LIB_FILE 0)

	if ( (MATIO_LIB_FILE) AND (MATIO_INC_DIR) )
		message(STATUS "matio lib is here : ${MATIO_LIB_FILE}")
		message(STATUS "--testAL----------------------------------------------")
		message(STATUS "matio include is here : ${MATIO_INC_DIR}")
	else()
		message(STATUS "ERROR !!! matio is not installed !!!!!")
	endif()

message(STATUS "------------------------------------------------")
message(STATUS "------------------------------------------------")	
################################################################
endif()



