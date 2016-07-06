if (XML_INC_DIR AND XML2_LIB_FILE AND BUILD_READ_XML_FILE)	

include_directories( ${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${EIGEN_INC_DIR} ${FAUST_XML_SRC_DIR} ${XML_INC_DIR} ${FAUST_CMDLINE_TYPE_FORMAT_XML_SRC_DIR})


#compatibilite ubuntu avec openblas
if (BUILD_OPENBLAS)
	include_directories(${OPENBLAS_INC_DIR} )
endif(BUILD_OPENBLAS)

	foreach(algofile launch_hierarchical_fact launch_palm4MSA)
		configure_file(${FAUST_CMDLINE_FUNC_SRC_DIR}/${algofile}.cpp ${FAUST_BIN_CMDLINE_SRC_DIR}/${algofile}.cpp COPYONLY)
		add_executable(${algofile} ${FAUST_BIN_CMDLINE_SRC_DIR}/${algofile}.cpp ${FAUST_CMDLINE_TYPE_FORMAT_XML_SRC_DIR}/xml_utils.cpp)
		if(BUILD_OPENBLAS)
			target_link_libraries(${algofile} ${FAUST_TARGET} ${OPENBLAS_LIB_FILE} ${XML2_LIB_FILE})
		else()
			target_link_libraries(${algofile} ${FAUST_TARGET} ${XML2_LIB_FILE})
		endif()
	
	endforeach()
	#install(TARGETS launch_hierarchical_fact launch_palm4MSA DESTINATION ${FAUST_BIN_CMDLINE_BIN_DIR})
   
endif()  
