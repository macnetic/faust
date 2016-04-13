#TEST VERSION GPU (version using CUDA)


include_directories(${FAUST_MATRIX_CU_SRC_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR} ${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_PALM4MSA_CU_SRC_DIR} ${EIGEN_INC_DIR} ${FAUST_EXCEPTION_SRC_DIR} )

if(FAUST_USE_MATIO)
	include_directories(${FAUST_MATIO_SRC_DIR} ${MATIO_INC_DIR})
endif(FAUST_USE_MATIO)


add_definitions(-DCOMPILE_GPU)
# generation des executables de algorithme hierarchical_fact en GPU en simple et double precision 
foreach(TEST_FPP float double)	
		foreach(testin hierarchical_fact_cu)
			# copy CPU files cpp.in to the user's ./src/ directory en float et double precision
			configure_file(${FAUST_TESTSRC_SRC_DIR}/${testin}.in ${FAUST_TESTSRC_BIN_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
			# Creation des executable en double et en float
			add_executable(${testin}_${TEST_FPP} ${FAUST_TESTSRC_BIN_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_TESTSRC_BIN_DIR}/ ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
			target_link_libraries(${testin}_${TEST_FPP} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${CUBLAS_LIB_FILE} ${CUDART_LIB_FILE}  ${CUSPARSE_LIB_FILE})

			install(TARGETS ${testin}_${TEST_FPP}  DESTINATION ${FAUST_INSTALL_TESTING_BIN})	
	endforeach()
endforeach()

# test GPU
add_test(NAME HADAMARD_GPU_FACT_DOUBLE COMMAND  ${FAUST_TEST_GPU_BIN_DIR}/hierarchical_fact_cu_double 	${FAUST_TESTDATA_SRC_DIR}/config_HADAMARD.mat)
add_test(NAME FAUST_HIER_GPU_DOUBLE COMMAND ${FAUST_TEST_GPU_BIN_DIR}/hierarchical_fact_cu_double 	${FAUST_TESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
add_test(NAME FAUST_HIER_GPU_FLOAT COMMAND ${FAUST_TEST_GPU_BIN_DIR}/hierarchical_fact_cu_float ${FAUST_TESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
add_test(NAME MEG_GPU_FACT_FLOAT COMMAND ${FAUST_TEST_GPU_BIN_DIR}/hierarchical_fact_cu_float ${FAUST_TESTDATA_SRC_DIR}/config_MEG.mat 22332 1)
add_test(NAME MEG_GPU_FACT_DOUBLE COMMAND ${FAUST_TEST_GPU_BIN_DIR}/hierarchical_fact_cu_double ${FAUST_TESTDATA_SRC_DIR}/config_MEG.mat 22480 1)
set(TIMEOUT_MEG 14000)
set_tests_properties(MEG_GPU_FACT_FLOAT MEG_GPU_FACT_DOUBLE PROPERTIES TIMEOUT ${TIMEOUT_MEG})



