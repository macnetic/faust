
if(MATIO_LIB_FILE AND HDF5_LIB_FILE AND MATIO_INC_DIR)
   include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_MATIO_SRC_DIR} ${MATIO_INC_DIR})

   # palm
   add_executable(palm ${FAUST_PALM4MSATESTSRC_SRC_DIR}/palm4MSA_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(palm ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(palm ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   # hier_fact
   add_executable(faust_hier ${FAUST_PALM4MSATESTSRC_SRC_DIR}/hierarchical_fact_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(faust_hier ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(faust_hier ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   # meg
   add_executable(meg ${FAUST_PALM4MSATESTSRC_SRC_DIR}/MEG_fact.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(meg ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(meg ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   # faust
   add_executable(faust_test ${FAUST_PALM4MSATESTSRC_SRC_DIR}/faust_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_params.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(faust_test ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE} )
   else()
      target_link_libraries(faust_test ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   CONFIGURE_FILE(${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_compared_palm2.mat ${FAUST_PALM4MSATESTDATA_BIN_DIR}/config_compared_palm2.mat COPYONLY)
   CONFIGURE_FILE(${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_MEG.mat ${FAUST_PALM4MSATESTDATA_BIN_DIR}/config_MEG.mat COPYONLY)
   CONFIGURE_FILE(${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat ${FAUST_PALM4MSATESTDATA_BIN_DIR}/config_compared_hierarchical_fact.mat COPYONLY)

   #install(TARGETS palm faust_hier meg faust_test DESTINATION ${FAUST_INSTALL_TESTING_BIN})
   install(TARGETS faust_hier meg faust_test DESTINATION ${FAUST_INSTALL_TESTING_BIN})
   install(FILES ${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_compared_palm2.mat 
                 ${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_MEG.mat 
                 ${FAUST_PALM4MSATESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 
           DESTINATION ${FAUST_INSTALL_TESTING_DATA})
   

   
   #time_comp
	add_executable(time_comp time_comp_spmat_vec.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp)	
	target_link_libraries(time_comp ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	
	#time_comp_mat
	add_executable(time_comp_mat time_comp_spmat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp)	
	target_link_libraries(time_comp_mat ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	
else()	

	#test mkl faust_core multiply
	#message(STATUS "MMMMMMMMMMMMMMMMMMMMMMMMMMM : MKL TEST")	
	# include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR})
	#add_executable(fcore_test faust_core_test.cpp)	
	#target_link_libraries(fcore_test ${FAUST_TARGET})
   
   
   
endif()
