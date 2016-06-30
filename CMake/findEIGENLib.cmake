### find and/or install eigen lib (cf. http://eigen.tuxfamily.org/index.php?title=Main_Page)
set(INCLUDE_PATH_LIST ${LIBRARY_PATH_LIST_TMP_DEFAULT}) # CACHE PATH "List of include paths used as PATH parameter in find_path")
check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0) # on test avec INCLUDE_PATH_LIST_TMP_DEFAULT
if (EIGEN_INC_DIR)
	message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
else (EIGEN_INC_DIR)
	if(UNIX)
		#exec_program("which hg | xargs echo" OUTPUT_VARIABLE MERCURIAL_HG_PATH_DIR)
	   	#exec_program("readlink ${MERCURIAL_HG_PATH_DIR}" OUTPUT_VARIABLE READLINK_TMP)
		#sudo dnf install eigen3-devel
		#exec_program("hg clone https://bitbucket.org/eigen/eigen/ ${CMAKE_SOURCE_DIR}/externals/eigen")
		
		# Eigen library come from following link. You can modifying svn version of eigen for latest stable release version.
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/tmp https://bitbucket.org/eigen/eigen/get/3.2.8.tar.bz2")
		set(EIGEN_LIB_NAME "3.2.8.tar.bz2") 
		exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/unix/tarLibs/${EIGEN_LIB_NAME} -C ${CMAKE_SOURCE_DIR}/externals/unix")
		exec_program("mv ${CMAKE_SOURCE_DIR}/externals/unix/eigen-* ${CMAKE_SOURCE_DIR}/externals/unix/eigen")

		add_include_path(INCLUDE_PATH_LIST_TMP_EIGEN "${PROJECT_SOURCE_DIR}/externals/unix/eigen")
	elseif(WIN32)
		#exec_program("wget -P ${CMAKE_SOURCE_DIR}/externals/win/eigen http://bitbucket.org/eigen/eigen/get/3.2.8.zip")
		set(EIGEN_LIB_NAME "eigen-eigen-07105f7124f9.zip")
		exec_program("tar jxf ${CMAKE_SOURCE_DIR}/externals/win/zipLibs/${EIGEN_LIB_NAME} ${CMAKE_SOURCE_DIR}/externals/win/eigen/")
		add_include_path(INCLUDE_PATH_LIST_TMP_EIGEN "${PROJECT_SOURCE_DIR}/externals/win/eigen")
	else(UNIX)
		message(WARNING "Unknown type of plateform for library Eigen")	
	endif(UNIX)
	

	#add_include_path(INCLUDE_PATH_LIST_TMP2 "${PROJECT_SOURCE_DIR}/externals/eigen")
	#set(INCLUDE_PATH_LIST2 ${INCLUDE_PATH_LIST_TMP2} CACHE PATH "List of include paths used as PATH parameter in find_path")
	set(INCLUDE_PATH_LIST ${INCLUDE_PATH_LIST_TMP_EIGEN})# CACHE PATH "List of include paths used as PATH parameter in find_path")
	#message(STATUS "INCLUDE_PATH_LIST=${INCLUDE_PATH_LIST}")
	check_external_includes("Eigen/Dense" EIGEN_INC_DIR 0)
	message(STATUS "Eigen library is available here: ${EIGEN_INC_DIR}")
endif (EIGEN_INC_DIR)
#############################################################

