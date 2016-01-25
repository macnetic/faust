

configure_file(${FAUST_MEXINTERFACETOOLS_SRC_DIR}/matlab_faust.m.in ${FAUST_TESTMATLABSCRIPT_BIN_DIR}/matlab_faust.m COPYONLY)


if(FAUST_USE_MATIO)

	if(FAUST_USE_XML)
		include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_MATIO_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${MATIO_INC_DIR} ${FAUST_XML_SRC_DIR} ${XML_INC_DIR})
	else()
		include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_MATIO_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${MATIO_INC_DIR})
	endif()
else()
	if(FAUST_USE_XML)
		include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR} ${FAUST_XML_SRC_DIR} ${XML_INC_DIR})
	else()
		include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${FAUST_EIGEN_SRC_DIR} ${FAUST_EXCEPTION_SRC_DIR})
	endif()
endif()

if(MATIO_LIB_FILE AND HDF5_LIB_FILE AND MATIO_INC_DIR AND FAUST_USE_MATIO)
   

   foreach(cppin faust_test hierarchical_fact_test MEG_fact palm4MSA_test)
      configure_file(${FAUST_TESTSRC_SRC_DIR}/${cppin}.cpp.in ${FAUST_TESTSRC_BIN_DIR}/${cppin}.cpp @ONLY)
   endforeach()
   #CONFIGURE_FILE(${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp COPYONLY)
   #CONFIGURE_FILE(${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp COPYONLY)
   add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp" "${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp")
   add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp" "${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp")
   add_custom_command(OUTPUT ${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m COMMAND cmake -E copy "${FAUST_TESTRUNCOMP_SRC_DIR}/drawComptime1.m" "${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m" DEPENDS "${FAUST_TESTRUNCOMP_SRC_DIR}/drawComptime1.m")
	#CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/main_RCG.m.in ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m @ONLY)


   # palm
   add_executable(palm ${FAUST_TESTSRC_BIN_DIR}/palm4MSA_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(palm ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(palm ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()


   # hier_fact
   add_executable(faust_hier ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(faust_hier ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(faust_hier ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   # meg
   add_executable(meg ${FAUST_TESTSRC_BIN_DIR}/MEG_fact.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(meg ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(meg ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   # faust
   add_executable(faust_test ${FAUST_TESTSRC_BIN_DIR}/faust_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_params.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(faust_test ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE} )
   else()
      target_link_libraries(faust_test ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()
    
	#comptime1
	add_custom_target(drawComp DEPENDS ${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m)
	add_executable(comptime1 ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_core.cpp)
	add_dependencies(comptime1 drawComp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(comptime1 ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE} )
   else()
      target_link_libraries(comptime1 ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()
	
	
	
   # comptime0
   add_executable(comptime0 ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_core.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(comptime0 ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(comptime0 ${MATLAB_ROOT}/bin/glnxa64/libmex.so ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
   endif()

   CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/main_RCG.m.in ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m @ONLY)
   INSTALL(FILES ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m DESTINATION ${FAUST_INSTALL_TESTRUNCOMP})
   FOREACH(matlabscript demo_compTime.m demo_compTime_matrix_draw.m dvp.m figures_compTime.m)
      CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/${matlabscript} ${FAUST_TESTRUNCOMP_BIN_DIR}/${matlabscript} COPYONLY)
      INSTALL(FILES ${FAUST_TESTRUNCOMP_BIN_DIR}/${matlabscript} DESTINATION ${FAUST_INSTALL_TESTRUNCOMP})
   ENDFOREACH()



   install(TARGETS  palm faust_hier meg faust_test DESTINATION ${FAUST_INSTALL_TESTING_BIN})
   #install(TARGETS faust_hier meg faust_test comptime0 DESTINATION ${FAUST_INSTALL_TESTING_BIN})
   #install(FILES ${FAUST_TESTDATA_SRC_DIR}/config_compared_palm2.mat 
   #              ${FAUST_TESTDATA_SRC_DIR}/config_MEG.mat 
   #              ${FAUST_TESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 
   #        DESTINATION ${FAUST_INSTALL_TESTING_DATA})
   

   
   #time_comp
	#add_executable(time_comp time_comp_spmat_vec.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp)	
	#target_link_libraries(time_comp ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	
	#time_comp_mat
	#add_executable(time_comp_mat time_comp_spmat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp)	
	#target_link_libraries(time_comp_mat ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	
else()	
	### IN CASE MATIO IS NOT AVAILABLE YOU CAN RUN THE TEST
	foreach(cppinbis hierarchical_fact_test MEG_fact)
      configure_file(${FAUST_TESTSRC_SRC_DIR}/${cppinbis}.cpp.inbis ${FAUST_TESTSRC_BIN_DIR}/${cppinbis}.cpp @ONLY)
	endforeach()
	foreach(datafile hier_fact_config MEG_config)
		configure_file(${FAUST_TESTDATA_SRC_DIR}/${datafile}.txt.in ${FAUST_TESTDATA_SRC_DIR}/${datafile}.txt @ONLY)
	endforeach()
	
	
	add_executable(faust_hier ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_test.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(faust_hier  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(faust_hier ${FAUST_TARGET} )
   endif()
   
   add_executable(meg ${FAUST_TESTSRC_BIN_DIR}/MEG_fact.cpp)
   if(FAUST_USE_OPENBLAS)
      target_link_libraries(meg  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
   else()
      target_link_libraries(meg  ${FAUST_TARGET})
   endif()
   
   
   
endif()

if (XML_INC_DIR AND XML2_LIB_FILE AND FAUST_USE_XML)	
	foreach(xmldatafile hier_fact_config MEG_config palm4MSA_config)
		configure_file(${FAUST_TESTDATA_SRC_DIR}/${xmldatafile}.xml.in ${FAUST_TESTDATA_SRC_DIR}/${xmldatafile}.xml @ONLY)
	endforeach()
	
	foreach(algofile launch_hierarchical_fact launch_palm4MSA)
		configure_file(${FAUST_TESTSRC_SRC_DIR}/${algofile}.cpp ${FAUST_TESTSRC_BIN_DIR}/${algofile}.cpp COPYONLY)
		add_executable(${algofile} ${FAUST_TESTSRC_BIN_DIR}/${algofile}.cpp ${FAUST_XML_SRC_DIR}/xml_utils.cpp)
		if(FAUST_USE_OPENBLAS)
			target_link_libraries(${algofile}  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE} ${XML2_LIB_FILE})
		else()
			target_link_libraries(${algofile}  ${FAUST_TARGET} ${XML2_LIB_FILE})
		endif()
	endforeach()



   
endif()   
   


	
	# add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/hierarchical_fact_basic.cpp" "${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/hierarchical_fact_basic.cpp")
	
	 # add_executable(faust_hier_basic ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp)
   # if(FAUST_USE_OPENBLAS)
      # target_link_libraries(faust_hier_basic  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
   # else()
      # target_link_libraries(faust_hier_basic ${FAUST_TARGET})
   # endif()
