###### chek and find Openblas library
set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_DEFAULT})# CACHE PATH "List of include paths used as PATH parameter in find_path")

check_external_libraries(xml2 XML2_LIB_FILE 0)
check_external_includes("libxml/parser.h" XML_INC_DIR 0)

if ( (XML2_LIB_FILE) AND (XML_INC_DIR) )
		message(STATUS "xml2 lib is here : ${XML2_LIB_FILE}")
		message(STATUS "libxml/parser.h is here : ${XML_INC_DIR}")
else ( (XML2_LIB_FILE) AND (XML_INC_DIR) )
	if(UNIX)
		message(STATUS "------------------------------------------------")
		message(STATUS "------------------------------------------------")
		message(STATUS "------------ XML2 LIB INSTALLATION ------------")
		message(STATUS "------------------------------------------------")
		message(STATUS "------------------------------------------------")
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		set(XML2_LIB_NAME "libxml2-2.9.4.tar.gz")
		exec_program("tar xzf ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${XML2_LIB_NAME} -C ${CMAKE_SOURCE_DIR}/externals/unix")
		exec_program("rm -r ${CMAKE_SOURCE_DIR}/externals/unix/libxml2")
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/libxml2-* ${CMAKE_SOURCE_DIR}/externals/unix/libxml2")
		exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/libxml2 && chmod -R 777 ./ && ./configure && make ")
		exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/libxml2 && sudo make install")
		#exec_program("cd ${CMAKE_SOURCE_DIR}/externals/unix/libxml2 && make install PREFIX='${CMAKE_SOURCE_DIR}/externals/unix/libxml2' ")
		# NOTE WARNING : WE DON'T run make install because it is not a root user install.  We used directly the lib and include in sdk_matio package. 
		# NOTE libxml2 install :  brew install libxml2 / apt_get / dnf /...
		add_include_path(INCLUDE_PATH_LIST_TMP_XML2 "/usr/include/libxml2" "${PROJECT_SOURCE_DIR}/externals/unix/libxml2/include")
		add_library_path(LIBRARY_PATH_LIST_TMP_XML2 "/usr/local/lib" "${PROJECT_SOURCE_DIR}/externals/unix/libxml2/.libs")	

		#add_include_path(INCLUDE_PATH_LIST_TMP_XML2 "/usr/include/libxml2")
		#add_library_path(LIBRARY_PATH_LIST_TMP_XML2 "/usr/lib")

	elseif(WIN32)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win http://github.com/xianyi/OpenBLAS/archive/v0.2.18.tar.gz")
		#exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/v0.2.18.tar.bz -C ${CMAKE_SOURCE_DIR}/externals/win")
		#exec_program("mv ${CMAKE_SOURCE_DIR}/externals/win/sdk_openBlas* ${CMAKE_SOURCE_DIR}/externals/win/OpenBLAS")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library OpenBlas")	
	endif(UNIX)

	set(LIBRARY_PATH_LIST ${LIBRARY_PATH_LIST_TMP_XML2}) # CACHE PATH "List of library paths used as PATH parameter in find_library")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_XML2})# CACHE PATH "List of include paths used as PATH parameter in find_path")
	check_external_libraries(xml2 XML2_LIB_FILE 0)
	check_external_includes("libxml/parser.h" XML_INC_DIR 0)

	if (  (XML2_LIB_FILE) AND (XML_INC_DIR) )
		message(STATUS "xml2 lib is here : ${XML2_LIB_FILE}")
		message(STATUS "libxml/parser.h is here : ${XML_INC_DIR}")
	else()
		message(STATUS "ERROR !!! xml2 is not installed !!!!!")
		message(FATAL_ERROR "xml2 lib is not installed on your system. Please check xml2 install.")	
	endif()

message(STATUS "------------------------------------------------")
message(STATUS "------------------------------------------------")	
################################################################
endif()



