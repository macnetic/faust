
file(GLOB MATLAB_SCRIPTS "*.m")







if(FAUST_USE_MEX)

	configure_file(${FAUST_DEMO_SRC_DIR}/run_all_demo.m ${FAUST_DEMO_BIN_DIR}/run_all_demo.m COPYONLY)
      		install(FILES ${FAUST_DEMO_BIN_DIR}/run_all_demo.m DESTINATION ${FAUST_DEMO_INSTALL_DIR})

	
	############# Brain source localization demo ####################
	file(GLOB BSL_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_BSL_SRC_DIR} "${FAUST_DEMO_BSL_SRC_DIR}/*.m")	
	foreach(matlab_script ${BSL_MATLAB_SCRIPTS})
      		configure_file(${FAUST_DEMO_BSL_SRC_DIR}/${matlab_script} ${FAUST_DEMO_BSL_BIN_DIR}/${matlab_script} COPYONLY)
      		install(FILES ${FAUST_DEMO_BSL_BIN_DIR}/${matlab_script} DESTINATION ${FAUST_DEMO_BSL_INSTALL_DIR})	
   	endforeach()
	
	file(GLOB BSL_DATA_FILE  RELATIVE ${FAUST_DEMO_BSL_DATA_SRC_DIR} "${FAUST_DEMO_BSL_DATA_SRC_DIR}/*.mat")
	foreach(data_BSL ${BSL_DATA_FILE})
      		configure_file(${FAUST_DEMO_BSL_DATA_SRC_DIR}/${data_BSL} ${FAUST_DEMO_BSL_DATA_BIN_DIR}/${data_BSL} COPYONLY)
      		install(FILES ${FAUST_DEMO_BSL_DATA_BIN_DIR}/${data_BSL} DESTINATION ${FAUST_DEMO_BSL_DATA_INSTALL_DIR})	
   	endforeach()
	

	######### Hadamard factorization ###############################
	file(GLOB HADAMARD_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_HADAMARD_SRC_DIR} "${FAUST_DEMO_HADAMARD_SRC_DIR}/*.m")
	foreach(matlab_script ${HADAMARD_MATLAB_SCRIPTS})		
		configure_file(${FAUST_DEMO_HADAMARD_SRC_DIR}/${matlab_script} ${FAUST_DEMO_HADAMARD_BIN_DIR}/${matlab_script} COPYONLY)
      		install(FILES ${FAUST_DEMO_HADAMARD_BIN_DIR}/${matlab_script} DESTINATION ${FAUST_DEMO_HADAMARD_INSTALL_DIR})
	endforeach()

	####### Runtime Comparison ######################################
	file(GLOB TIMECOMP_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_TIMECOMP_SRC_DIR} "${FAUST_DEMO_TIMECOMP_SRC_DIR}/*.m")
	foreach(matlab_script ${TIMECOMP_MATLAB_SCRIPTS})		
		configure_file(${FAUST_DEMO_TIMECOMP_SRC_DIR}/${matlab_script} ${FAUST_DEMO_TIMECOMP_BIN_DIR}/${matlab_script} COPYONLY)
      		install(FILES ${FAUST_DEMO_TIMECOMP_BIN_DIR}/${matlab_script} DESTINATION ${FAUST_DEMO_TIMECOMP_INSTALL_DIR})
	endforeach()
	
	
	#################### Tools directory ############################
	file(GLOB TOOLS_MATLAB_SCRIPTS RELATIVE ${FAUST_DEMO_TOOLS_SRC_DIR}  "${FAUST_DEMO_TOOLS_SRC_DIR}/*.m")
	foreach(tools ${TOOLS_MATLAB_SCRIPTS})
      		configure_file(${FAUST_DEMO_TOOLS_SRC_DIR}/${tools} ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools} COPYONLY)
      		install(FILES ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools} DESTINATION ${FAUST_DEMO_TOOLS_INSTALL_DIR})	
   	endforeach()
		
	


endif()
	
	
