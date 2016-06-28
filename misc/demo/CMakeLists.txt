file(GLOB MATLAB_SCRIPTS "*.m")
file(GLOB MATLAB_SCRIPTS "*.m")







if(FAUST_USE_MEX)
	
	############# Brain source localization demo ####################
	foreach(matlab_script BSL Fig_BSL)
      		configure_file(${FAUST_DEMO_BSL_SRC_DIR}/${matlab_script}.m ${FAUST_DEMO_BSL_BIN_DIR}/${matlab_script}.m COPYONLY)
      		install(FILES ${FAUST_DEMO_BSL_BIN_DIR}/${matlab_script}.m DESTINATION ${FAUST_DEMO_BSL_INSTALL_DIR})	
   	endforeach()
	
	
	foreach(data_BSL curv M_6 M_8 M_16 M_25 X_meg)
      		configure_file(${FAUST_DEMO_BSL_DATA_SRC_DIR}/${data_BSL}.mat ${FAUST_DEMO_BSL_DATA_BIN_DIR}/${data_BSL}.mat COPYONLY)
      		install(FILES ${FAUST_DEMO_BSL_DATA_BIN_DIR}/${data_BSL}.mat DESTINATION ${FAUST_DEMO_BSL_DATA_INSTALL_DIR})	
   	endforeach()
	

	


endif()
	
	#################### Tools directory ############################
	foreach(tools bplot dvp faust_transpose f_mult greed_omp_chol make_sparse sparse_coeffs UpdateCholesky)
      		configure_file(${FAUST_DEMO_TOOLS_SRC_DIR}/${tools}.m ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools}.m COPYONLY)
      		install(FILES ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools}.m DESTINATION ${FAUST_DEMO_TOOLS_INSTALL_DIR})	
   	endforeach()
