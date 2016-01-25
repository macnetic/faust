# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
configure_file(${FAUST_DOC_SRC_DIR}/Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile @ONLY)
add_custom_target(doc ALL
${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
COMMENT "Generating API documentation with Doxygen" VERBATIM)
install(DIRECTORY ${PROJECT_BINARY_DIR}/html DESTINATION share/doc)


endif(DOXYGEN_FOUND)

 