#TEST VERSION GPU (version using CUDA)


include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_SRC_LINEAR_OPERATOR_GPU_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${EIGEN_INC_DIR})

if(BUILD_READ_MAT_FILE)
	include_directories(${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR} ${MATIO_INC_DIR})
endif(BUILD_READ_MAT_FILE)


add_definitions(-DCOMPILE_GPU)
# generation des executables de algorithme hierarchical_fact en GPU en simple et double precision 
foreach(TEST_FPP float double)	
		foreach(testin hierarchicalFactorization_gpu)
			# copy CPU files cpp.in to the user's ./src/ directory en float et double precision
			configure_file(${FAUST_SRC_TEST_SRC_DIR}/${testin}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
			# Creation des executable en double et en float
			add_executable(${testin}_${TEST_FPP} ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_BIN_TEST_SRC_DIR}/ ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio_mat.cpp)
			target_link_libraries(${testin}_${TEST_FPP} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${CUBLAS_LIB_FILE} ${CUDART_LIB_FILE}  ${CUSPARSE_LIB_FILE})

			install(TARGETS ${testin}_${TEST_FPP}  DESTINATION ${FAUST_INSTALL_BIN})	
	endforeach()
endforeach()

# test GPU
add_test(NAME HADAMARD_GPU_FACT_DOUBLE COMMAND  ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_double 	${FAUST_DATA_MAT_DIR}/config_HADAMARD.mat)
add_test(NAME FAUST_HIER_GPU_DOUBLE COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_double 	${FAUST_DATA_MAT_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
add_test(NAME FAUST_HIER_GPU_FLOAT COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_float ${FAUST_DATA_MAT_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
add_test(NAME MEG_GPU_FACT_FLOAT COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_float ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22332 1)
add_test(NAME MEG_GPU_FACT_DOUBLE COMMAND ${FAUST_BIN_TEST_GPU_DIR}/hierarchicalFactorization_gpu_double ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22480 1)
set(TIMEOUT_MEG 14000)
set_tests_properties(MEG_GPU_FACT_FLOAT MEG_GPU_FACT_DOUBLE PROPERTIES TIMEOUT ${TIMEOUT_MEG})



