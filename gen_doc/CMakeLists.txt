##Generation of the documentation using doxygen
# add a target to generate API documentation with Doxygen
if(BUILD_DOCUMENTATION)
	find_package(Doxygen)
	if(DOXYGEN_FOUND)
		configure_file(${FAUST_DOC_SRC_DIR}/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/Doxyfile @ONLY)
		# ./gen_doc/images/* files is duplicated in doc/html/ to call images documentation in the source code with relative path of image's files, from build directory.     
		if(UNIX)
			configure_file(${FAUST_DOC_SRC_DIR}/py_filter.sh.in ${PROJECT_BINARY_DIR}/doc/py_filter @ONLY)
		elseif(WIN32)
			configure_file(${FAUST_DOC_SRC_DIR}/py_filter.bat.in ${PROJECT_BINARY_DIR}/doc/py_filter.bat @ONLY)
		endif()
		file(COPY ${FAUST_DOC_SRC_DIR}/images/ DESTINATION ${PROJECT_BINARY_DIR}/doc/html/ )
		add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc COMMENT "Generating API documentation with Doxygen" VERBATIM)
		#install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION doc/share/doc)
		#file(MAKE_DIRECTORY ${FAUST_INSTALL_DOC})
		#install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION ${FAUST_INSTALL_DOC} FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)


	else()
		message(WARNING "Doxygen not found, cannot generate the documentation")
	endif(DOXYGEN_FOUND)
endif(BUILD_DOCUMENTATION) 





#IF (DOXYGEN_FOUND)
#    IF (NOT DOXY_IN_FILE_UPDATED)
#		IF(DOXYGEN_DOT_FOUND)
#    		IF(BUILD_DOCUMENTATION_GUI)
#				CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/doc/doxygen_configGUIWITHDOT.in ${PROJECT_BINARY_DIR}/doc/doxygen_config.in @ONLY) 
#			ELSE(BUILD_DOCUMENTATION_GUI)
#				CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/doc/doxygen_configSTDWITHDOT.in ${PROJECT_BINARY_DIR}/doc/doxygen_config.in @ONLY)    
#			ENDIF(BUILD_DOCUMENTATION_GUI)
#		ELSE(DOXYGEN_DOT_FOUND)
#			IF(BUILD_DOCUMENTATION_GUI)
#				CONFIGURE_FILE(${PROJECT_BINARY_DIR}/doc/doxygen_configGUI.in ${PROJECT_BINARY_DIR}/doc/doxygen_config.in @ONLY) 
#			ELSE(BUILD_DOCUMENTATION_GUI)
#				CONFIGURE_FILE(${PROJECT_BINARY_DIR}/doc/doxygen_configSTD.in ${PROJECT_BINARY_DIR}/doc/doxygen_config.in @ONLY) 
#			ENDIF(BUILD_DOCUMENTATION_GUI)
#		ENDIF(DOXYGEN_DOT_FOUND)
#		SET(DOXY_IN_FILE_UPDATED 1 INTERNAL) 
#	ENDIF (NOT DOXY_IN_FILE_UPDATED)
#	ADD_CUSTOM_TARGET(refman_doxygen ALL echo DEPENDS faust)
#	ADD_CUSTOM_COMMAND(
#  		TARGET refman_doxygen
#		COMMAND   ${DOXYGEN}
#		ARGS "${PROJECT_BINARY_DIR}/doc/doxygen_config.in"
#		COMMENT   "Doxygen generate the documentation")
#ELSE (DOXYGEN_FOUND)
#	MESSAGE(STATUS "Doxygen not found, cannot generate the documentation")
#ENDIF (DOXYGEN_FOUND)
