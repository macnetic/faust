#TEST VERSION GPU (version using CUDA)


include_directories(${FAUST_MATRIX_CU_SRC_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR} ${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_PALM4MSA_CU_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} )

if(FAUST_USE_MATIO)
	include_directories(${FAUST_MATIO_SRC_DIR} ${MATIO_INC_DIR})
endif(FAUST_USE_MATIO)

#rajout OPENBLAS_INC_DIR pour la compatibilite UBUNTU
#if(FAUST_USE_OPENBLAS)
#	include_directories(${OPENBLAS_INC_DIR})
#endif()


# copy GPU files cpp.in to the user's ./src/ directory
foreach(gpucppin hierarchical_fact_test_cu MEG_fact_cu) 
	#configure_file(${FAUST_TESTSRC_SRC_DIR}/${gpucppin}.cpp.in ${FAUST_TESTSRC_SRC_DIR}/${gpucppin}.cpp @ONLY)
	configure_file(${FAUST_TESTSRC_SRC_DIR}/${gpucppin}.cpp.in ${FAUST_TESTSRC_BIN_DIR}/${gpucppin}.cpp @ONLY)
	# USING GPU (cuda)
	#add_definitions(-D__COMPILE_GPU__ -D__COMPILE_SPMAT__)
	add_executable(${gpucppin} ${FAUST_TESTSRC_BIN_DIR}/${gpucppin}.cpp ${FAUST_TESTSRC_BIN_DIR}/ ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)

	target_link_libraries(${gpucppin} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${CUBLAS_LIB_FILE} ${CUDART_LIB_FILE}  ${CUSPARSE_LIB_FILE})

endforeach()


install(TARGETS  hierarchical_fact_test_cu MEG_fact_cu DESTINATION ${FAUST_INSTALL_TESTING_BIN})


