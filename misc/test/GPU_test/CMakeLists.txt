#TEST VERSION GPU (version using CUDA)


include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_SRC_LINEAR_OPERATOR_GPU_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${EIGEN_INC_DIR})

if(BUILD_READ_MAT_FILE)
	include_directories(${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR} ${MATIO_INC_DIR})
endif(BUILD_READ_MAT_FILE)


add_definitions(-DCOMPILE_GPU)


if(BUILD_READ_MAT_FILE AND MATIO_LIB_FILE AND MATIO_INC_DIR ) 
	# generation des executables de algorithme hierarchical_fact en GPU en simple et double precision 
	foreach(TEST_FPP float double)	
		foreach(testin hierarchicalFactorization_gpu ) #Faust_multiplication_gpu
			# copy CPU files cpp.in to the user's ./src/ directory en float et double precision
			configure_file(${FAUST_SRC_TEST_SRC_DIR}/${testin}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
			# Creation des executable en double et en float
			add_executable(${testin}_${TEST_FPP} ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_BIN_TEST_SRC_DIR}/ ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio_mat.cpp)
			target_link_libraries(${testin}_${TEST_FPP} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${CUBLAS_LIB_FILE} ${CUDART_LIB_FILE}  ${CUSPARSE_LIB_FILE})

			install(TARGETS ${testin}_${TEST_FPP} DESTINATION ${FAUST_INSTALL_BIN})	
		endforeach()
	endforeach()
###### IN CASE MATIO IS NOT AVAILABLE YOU CAN RUN THE TEST ######
else(BUILD_READ_MAT_FILE AND MATIO_LIB_FILE AND MATIO_INC_DIR )	
	#message(STATUS "Library matio is not available. Please used .xml input files.")
endif(BUILD_READ_MAT_FILE AND MATIO_LIB_FILE AND MATIO_INC_DIR )


set(TIMEOUT_MEG 14000)
# MATLAB input format
if (BUILD_READ_MAT_FILE AND MATIO_LIB_FILE AND MATIO_INC_DIR)
	foreach(TEST_FPP float double)
#message (STATUS "------------------------------------------------")
#message (STATUS "${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_HADAMARD_32.mat ${FAUST_CONFIG_MAT_DIR}/config_HADAMARD_32.mat 5930 30")

		add_test(NAME MATFILE_FACT_HADAMARD_GPU_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_HADAMARD_32.mat ${FAUST_CONFIG_MAT_DIR}/config_HADAMARD_32.mat 5930 30)
		# add_test(NAME MATFILE_FACT_HADAMARD_GPU_FLOAT COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_float ${FAUST_DATA_MAT_DIR}/matrix_HADAMARD_32.mat ${FAUST_CONFIG_MAT_DIR}/config_HADAMARD_32.mat 5930 30)

	add_test(NAME MATFILE_FACT_HIER_GPU_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_hierarchical_fact.mat ${FAUST_CONFIG_MAT_DIR}/config_hierarchical_fact.mat 9401.5 0.1)
	#add_test(NAME MATFILE_FACT_HIER_GPU_FLOAT COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_float ${FAUST_DATA_MAT_DIR}/matrix_hierarchical_fact.mat ${FAUST_CONFIG_MAT_DIR}/config_hierarchical_fact.mat 9401.5 0.1)
	

	add_test(NAME MATFILE_FACT_MEG_GPU_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG.mat 22450 200)
	#add_test(NAME MATFILE_FACT_MEG_GPU_DOUBLE COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_double ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG.mat 22480 1)
	


	# add_test(NAME MATFILE_MULTIPLICATION_MEG_GPU_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_GPU_DIR}/Faust_multiplication_gpu_${TEST_FPP} ${FAUST_DEMO_BSL_DATA_SRC_DIR}X_meg.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_6.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_8.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_16.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_25.mat)
	# add_test(NAME MATFILE_MULTIPLICATION_MEG_GPU_FLOAT COMMAND ${FAUST_BIN_TEST_GPU_DIR}/Faust_multiplication_gpu_float ${FAUST_DEMO_BSL_DATA_SRC_DIR}X_meg.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_6.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_8.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_16.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_25.mat)

	set_tests_properties(MATFILE_FACT_MEG_GPU_${TEST_FPP} PROPERTIES TIMEOUT ${TIMEOUT_MEG})	
	endforeach()
	

endif(BUILD_READ_MAT_FILE AND MATIO_LIB_FILE AND MATIO_INC_DIR)



