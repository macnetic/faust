if (XML_INC_DIR AND XML2_LIB_FILE AND FAUST_USE_XML)	

include_directories( ${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${EIGEN_INC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_XML_SRC_DIR} ${XML_INC_DIR} ${FAUST_CMDLINE_TOOLS_SRC_DIR})


#compatibilite ubuntu avec openblas
if (FAUST_USE_OPENBLAS)
	include_directories(${OPENBLAS_INC_DIR} )
endif(FAUST_USE_OPENBLAS)

	foreach(algofile launch_hierarchical_fact launch_palm4MSA)
		configure_file(${FAUST_CMDLINE_FUNC_SRC_DIR}/${algofile}.cpp ${FAUST_CMDLINE_BIN_SRC_DIR}/${algofile}.cpp COPYONLY)
		add_executable(${algofile} ${FAUST_CMDLINE_BIN_SRC_DIR}/${algofile}.cpp ${FAUST_CMDLINE_TOOLS_SRC_DIR}/xml_utils.cpp)
		if(FAUST_USE_OPENBLAS)
			target_link_libraries(${algofile}  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE} ${XML2_LIB_FILE})
		else()
			target_link_libraries(${algofile}  ${FAUST_TARGET} ${XML2_LIB_FILE})
		endif()
	
	endforeach()
	#install(TARGETS launch_hierarchical_fact launch_palm4MSA DESTINATION ${FAUST_CMDLINE_BIN_DIR})
   
endif()  
