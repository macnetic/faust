
file(GLOB MATLAB_SCRIPTS "*.m")


if(BUILD_WRAPPER_MATLAB)

	configure_file(${FAUST_DEMO_SRC_DIR}/run_all_demo.m ${FAUST_DEMO_BIN_DIR}/run_all_demo.m COPYONLY)

	# copy BSL data matrices
	file(GLOB BSL_DATA_FILE  RELATIVE ${FAUST_DATA_MAT_DIR} "${FAUST_DATA_MAT_DIR}/*MEG*.mat")
	foreach(data_BSL ${BSL_DATA_FILE})
		configure_file(${FAUST_DATA_MAT_DIR}/${data_BSL} ${FAUST_MATFAUST_DEMO_DATA_BIN_DIR}/${data_BSL} COPYONLY)
   	endforeach()


	# copy data matrix for quick_start.m
	configure_file(${FAUST_DATA_MAT_DIR}/faust_quick_start.mat ${FAUST_MATFAUST_DEMO_DATA_BIN_DIR}/faust_quick_start.mat COPYONLY)



	#################### Tools directory ############################
	file(GLOB TOOLS_MATLAB_SCRIPTS RELATIVE ${FAUST_DEMO_TOOLS_SRC_DIR}  "${FAUST_DEMO_TOOLS_SRC_DIR}/*.m")
	foreach(tools ${TOOLS_MATLAB_SCRIPTS})
      		configure_file(${FAUST_DEMO_TOOLS_SRC_DIR}/${tools} ${FAUST_DEMO_TOOLS_BIN_DIR}/${tools} COPYONLY)
   	endforeach()




	#the installation is made in wrapper/matlab/CmakeList.txt as the demo is now Matlab code

endif()


