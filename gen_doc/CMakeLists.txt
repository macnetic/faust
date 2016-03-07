# add a target to generate API documentation with Doxygen
if(FAUST_GEN_DOC)
	find_package(Doxygen)
	if(DOXYGEN_FOUND)
		configure_file(${FAUST_DOC_SRC_DIR}/Doxyfile.in ${PROJECT_BINARY_DIR}/doc/Doxyfile @ONLY)
# ./gen_doc/images/* files is duplicated in doc/html/ to call images documentation in the source code with relative path of image's files, from build directory.     
		file(COPY ${FAUST_DOC_SRC_DIR}/images/ DESTINATION ${PROJECT_BINARY_DIR}/doc/html/ )
		add_custom_target(doc ALL ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doc/Doxyfile WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/doc COMMENT "Generating API documentation with Doxygen" VERBATIM)
		install(DIRECTORY ${PROJECT_BINARY_DIR}/doc/html DESTINATION doc/share/doc)
	endif(DOXYGEN_FOUND)
endif(FAUST_GEN_DOC) 



