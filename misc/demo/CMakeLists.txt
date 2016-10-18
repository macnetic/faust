
file(GLOB MATLAB_SCRIPTS "*.m")







if(BUILD_MATLAB_MEX_FILES)

	configure_file(${FAUST_DEMO_SRC_DIR}/run_all_demo.m ${FAUST_DEMO_BIN_DIR}/run_all_demo.m COPYONLY)

	
	############# Brain source localization demo ####################
	file(GLOB BSL_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_BSL_SRC_DIR} "${FAUST_DEMO_BSL_SRC_DIR}/*.m")	
	foreach(matlab_script ${BSL_MATLAB_SCRIPTS})
      		configure_file(${FAUST_DEMO_BSL_SRC_DIR}/${matlab_script} ${FAUST_DEMO_BSL_BIN_DIR}/${matlab_script} COPYONLY)
   	endforeach()
	
	file(GLOB BSL_DATA_FILE  RELATIVE ${FAUST_DATA_MAT_DIR} "${FAUST_DATA_MAT_DIR}/*MEG*.mat")
	foreach(data_BSL ${BSL_DATA_FILE})
      		configure_file(${FAUST_DATA_MAT_DIR}/${data_BSL} ${FAUST_DEMO_BSL_DATA_BIN_DIR}/${data_BSL} COPYONLY)
   	endforeach()
	

	######### Hadamard factorization ###############################
	file(GLOB HADAMARD_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_HADAMARD_SRC_DIR} "${FAUST_DEMO_HADAMARD_SRC_DIR}/*.m")
	foreach(matlab_script ${HADAMARD_MATLAB_SCRIPTS})		
		configure_file(${FAUST_DEMO_HADAMARD_SRC_DIR}/${matlab_script} ${FAUST_DEMO_HADAMARD_BIN_DIR}/${matlab_script} COPYONLY)
	endforeach()

	####### Runtime Comparison ######################################
	file(GLOB TIMECOMP_MATLAB_SCRIPTS  RELATIVE ${FAUST_DEMO_TIMECOMP_SRC_DIR} "${FAUST_DEMO_TIMECOMP_SRC_DIR}/*.m")
	foreach(matlab_script ${TIMECOMP_MATLAB_SCRIPTS})		
		configure_file(${FAUST_DEMO_TIMECOMP_SRC_DIR}/${matlab_script} ${FAUST_DEMO_TIMECOMP_BIN_DIR}/${matlab_script} COPYONLY)
	endforeach()


	##### Quick Start ##############################################
	file(GLOB TOOLS_MATLAB_SCRIPTS RELATIVE ${FAUST_DEMO_QUICKSTART_SRC_DIR}  "${FAUST_DEMO_QUICKSTART_SRC_DIR}/*.m")
	foreach(tools ${TOOLS_MATLAB_SCRIPTS})
		configure_file(${FAUST_DEMO_QUICKSTART_SRC_DIR}/${tools} ${FAUST_DEMO_QUICKSTART_BIN_DIR}/${tools} COPYONLY)

	endforeach()
	configure_file(${FAUST_DATA_MAT_DIR}/faust_quick_start.mat ${FAUST_DEMO_QUICKSTART_BIN_DIR}/faust_quick_start.mat COPYONLY)

	
	
	#################### Tools directory ############################
	file(GLOB TOOLS_MATLAB_SCRIPTS RELATIVE ${FAUST_DEMO_TOOLS_SRC_DIR}  "${FAUST_DEMO_TOOLS_SRC_DIR}/*.m")
	foreach(tools ${TOOLS_MATLAB_SCRIPTS})
      		configure_file(${FAUST_DEMO_TOOLS_SRC_DIR}/${tools} ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools} COPYONLY)
   	endforeach()
		
	install(DIRECTORY ${FAUST_DEMO_BIN_DIR} DESTINATION ${CMAKE_INSTALL_MATLAB_PREFIX} FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE DIRECTORY_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE  FILES_MATCHING PATTERN "*.m" PATTERN "*.mat" PATTERN "*CMakeFiles*" EXCLUDE)


endif()
	
	
