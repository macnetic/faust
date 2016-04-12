# WRAPPER : TESTS EXAMPLES USING THE LIBRARY FAUST WITH A CPU VERSION

include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} ${EIGEN_INC_DIR} ${FAUST_EXCEPTION_SRC_DIR})

if(FAUST_USE_MATIO)
	include_directories(${FAUST_MATIO_SRC_DIR} ${MATIO_INC_DIR})
endif(FAUST_USE_MATIO)

#rajout OPENBLAS_INC_DIR pour la compatibilite UBUNTU
if(FAUST_USE_OPENBLAS)
	include_directories(${OPENBLAS_INC_DIR})
endif()


####### IF THE MATIO LIBRARY IS AVAILABLE ELSE WE USE TXT OR XML DATA  #######
if(MATIO_LIB_FILE AND MATIO_INC_DIR AND FAUST_USE_MATIO AND HDF5_LIB_FILE)
	########## HIERARCHICAL FACTORIZATION IN CPU MODE (simple and double precision)
	foreach(TEST_FPP float double)	
		foreach(testin hierarchical_fact)
			# copy CPU files cpp.in to the user's ./src/ directory en float et double precision
			configure_file(${FAUST_TESTSRC_SRC_DIR}/${testin}.in ${FAUST_TESTSRC_BIN_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
			# Generation of the binary files in double and float precision
			add_executable(${testin}_${TEST_FPP} ${FAUST_TESTSRC_BIN_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_TESTSRC_BIN_DIR}/ ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
   	
			if(FAUST_USE_OPENBLAS)
		      		target_link_libraries(${testin}_${TEST_FPP}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
			else(FAUST_USE_OPENBLAS)
		      		target_link_libraries(${testin}_${TEST_FPP}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
			endif(FAUST_USE_OPENBLAS)	
		endforeach()
	endforeach()


### test_faust : Evaluate the time difference using faust or classical-dense matrix, in a context of the multiplication between an matrix and a vector
### WARNING les donnée test/data/facts_MEG%d.txt n'existent plus... 
### test_palm4MSA : Run a test of palm4MSA
	foreach(cppin test_palm4MSA) #test_faust
		configure_file(${FAUST_TESTSRC_SRC_DIR}/${cppin}.cpp.in ${FAUST_TESTSRC_BIN_DIR}/${cppin}.cpp @ONLY)
	endforeach()

	# test_palm4MSA
	add_executable(test_palm4MSA ${FAUST_TESTSRC_BIN_DIR}/test_palm4MSA.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
	if(FAUST_USE_OPENBLAS)
		target_link_libraries(test_palm4MSA  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
	else()
		target_link_libraries(test_palm4MSA  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	endif()

	# test_faust 
#	add_executable(test_faust ${FAUST_TESTSRC_BIN_DIR}/test_faust.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_params.cpp)
#	if(FAUST_USE_OPENBLAS)
#		target_link_libraries(test_faust  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE} )
#	else()
#		target_link_libraries(test_faust ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
#	endif()
    

	#CONFIGURE_FILE(${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp COPYONLY)
	#CONFIGURE_FILE(${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp COPYONLY)
	
	
	#### comp_time0 ???????????	
	#add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp" "${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/comp_time0.cpp")
	#### comp_time1 ???????????	
	#add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp" "${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/comp_time1.cpp")
	#### drawComptime1 ???????????	
	#add_custom_command(OUTPUT ${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m COMMAND cmake -E copy "${FAUST_TESTRUNCOMP_SRC_DIR}/drawComptime1.m" "${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m" DEPENDS "${FAUST_TESTRUNCOMP_SRC_DIR}/drawComptime1.m")
	
	
	#CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/main_RCG.m.in ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m @ONLY) 
	#comptime1
	#add_custom_target(drawComp DEPENDS ${FAUST_TESTRUNCOMP_BIN_DIR}/drawComptime1.m)
	#add_executable(comptime1 ${FAUST_TESTSRC_BIN_DIR}/comp_time1.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_core.cpp)
	#add_dependencies(comptime1 drawComp)
	#if(FAUST_USE_OPENBLAS)
	#	target_link_libraries(comptime1  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE} )
	#else()
	#	target_link_libraries(comptime1  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	#endif()	
	
	# comptime0
	#add_executable(comptime0 ${FAUST_TESTSRC_BIN_DIR}/comp_time0.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_core.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
	#if(FAUST_USE_OPENBLAS)
	#	target_link_libraries(comptime0  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
	#else()
	#	target_link_libraries(comptime0  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
	#endif()

	CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/main_RCG.m.in ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m @ONLY)
   	INSTALL(FILES ${FAUST_TESTRUNCOMP_BIN_DIR}/main_RCG.m DESTINATION ${FAUST_INSTALL_TESTRUNCOMP})
	foreach(matlabscript demo_compTime.m demo_compTime_matrix_draw.m dvp.m figures_compTime.m)
      		CONFIGURE_FILE(${FAUST_TESTRUNCOMP_SRC_DIR}/${matlabscript} ${FAUST_TESTRUNCOMP_BIN_DIR}/${matlabscript} COPYONLY)
		INSTALL(FILES ${FAUST_TESTRUNCOMP_BIN_DIR}/${matlabscript} DESTINATION ${FAUST_INSTALL_TESTRUNCOMP})
	endforeach()


	install(TARGETS hierarchical_fact_double hierarchical_fact_float test_palm4MSA DESTINATION ${FAUST_INSTALL_TESTING_BIN})


###### IN CASE MATIO IS NOT AVAILABLE YOU CAN RUN THE TEST ######
else(MATIO_LIB_FILE AND MATIO_INC_DIR AND FAUST_USE_MATIO AND HDF5_LIB_FILE)	

	message(STATUT "Library matio is not available. Please used .txt or .xml input files")
#	foreach(cppinbis hierarchical_fact_test MEG_fact)
#		configure_file(${FAUST_TESTSRC_SRC_DIR}/${cppinbis}.cpp.inbis ${FAUST_TESTSRC_BIN_DIR}/${cppinbis}.cpp @ONLY)
#	endforeach()
#	foreach(datafile hier_fact_config MEG_config)
#		configure_file(${FAUST_TESTDATA_SRC_DIR}/${datafile}.txt.in ${FAUST_TESTDATA_SRC_DIR}/${datafile}.txt @ONLY)
#	endforeach()
#	
#	add_executable(faust_hier ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_test.cpp)
#   	
#   	if(FAUST_USE_OPENBLAS)
# 		target_link_libraries(faust_hier  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
#  	else()
#  		target_link_libraries(faust_hier ${FAUST_TARGET} )
#  	endif()
#   
#  	add_executable(meg ${FAUST_TESTSRC_BIN_DIR}/MEG_fact.cpp)
#   	
#   	if(FAUST_USE_OPENBLAS)
#    		target_link_libraries(meg  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
#   	else()
#		target_link_libraries(meg  ${FAUST_TARGET})
#	endif()
 
endif(MATIO_LIB_FILE AND MATIO_INC_DIR AND FAUST_USE_MATIO AND HDF5_LIB_FILE)



if (FAUST_USE_MEX)
	configure_file(${FAUST_TEST_SRC_TOOL_DIR}/set_path.m.in ${FAUST_TESTTOOLS_BIN_DIR}/set_path.m @ONLY)
	configure_file(${FAUST_TEST_SRC_TOOL_DIR}/launch_hierarchical_fact.m ${FAUST_TESTTOOLS_BIN_DIR}/launch_hierarchical_fact.m COPYONLY)
	
	FOREACH(matlabtest hier_fact.m MEG_fact.m)
		CONFIGURE_FILE(${FAUST_TESTSRC_SRC_DIR}/${matlabtest} ${FAUST_TEST_BIN_DIR}/${matlabtest} COPYONLY)
	ENDFOREACH()
endif(FAUST_USE_MEX)



##############################################################################
# Testing with CTest
# CTest is a testing tool distributed as a part of CMake.
# Running the differents binary created in this project  
##############################################################################
# Critical time, which fixe the limit time of the process
set(TIMEOUT_MEG 7200)


# MATLAB input format
if (FAUST_USE_MATIO)
	add_test(NAME HADAMARD_FACT_DOUBLE COMMAND  ${FAUST_TEST_BIN_DIR}/hierarchical_fact_double 	${FAUST_TESTDATA_SRC_DIR}/config_HADAMARD.mat)
	add_test(NAME FAUST_HIER_DOUBLE COMMAND ${FAUST_TEST_BIN_DIR}/hierarchical_fact_double 	${FAUST_TESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
	add_test(NAME FAUST_HIER_FLOAT COMMAND ${FAUST_TEST_BIN_DIR}/hierarchical_fact_float ${FAUST_TESTDATA_SRC_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
	add_test(NAME MEG_FACT_FLOAT COMMAND ${FAUST_TEST_BIN_DIR}/hierarchical_fact_float ${FAUST_TESTDATA_SRC_DIR}/config_MEG.mat 22332 1)
	add_test(NAME MEG_FACT_DOUBLE COMMAND ${FAUST_TEST_BIN_DIR}/hierarchical_fact_double ${FAUST_TESTDATA_SRC_DIR}/config_MEG.mat 22480 1)
	set_tests_properties(MEG_FACT_FLOAT MEG_FACT_DOUBLE PROPERTIES TIMEOUT ${TIMEOUT_MEG})
endif (FAUST_USE_MATIO)

# XML input format
if (FAUST_USE_XML)
	add_test(NAME FAUST_HIER_CMDLINE COMMAND ${FAUST_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CMDLINE_SRC_CONFIGFILE_DIR}/hier_fact_config.xml ${FAUST_CMDLINE_SRC_DATA_DIR}/hier_fact_data.txt )

	add_test(NAME FAUST_HIER_TRANSPOSE_CMDLINE COMMAND ${FAUST_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CMDLINE_SRC_CONFIGFILE_DIR}/trans_hier_fact_config.xml ${FAUST_CMDLINE_SRC_DATA_DIR}/hier_fact_data.txt  T)

	add_test(NAME MEG_FACT_CMDLINE COMMAND ${FAUST_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CMDLINE_SRC_CONFIGFILE_DIR}/MEG_config.xml ${FAUST_CMDLINE_SRC_DATA_DIR}/MEG_data.txt )

	set_tests_properties(MEG_FACT_CMDLINE PROPERTIES TIMEOUT ${TIMEOUT_MEG})	 
endif(FAUST_USE_XML)


#if(FAUST_USE_MEX)
#	#add_test(NAME FAUST_HIER_MEX COMMAND matlab  -nojvm -nodisplay -r ${FAUST_TEST_BIN_DIR}/hier_fact.m;exit)
#	add_test(NAME FAUST_HIER_MEX COMMAND matlab -nojvm -r "addpath('${FAUST_TEST_BIN_DIR}');hier_fact;exit" )	
#	set_tests_properties(FAUST_HIER_MEX PROPERTIES TIMEOUT 60)
#endif() 

#matlab -nojvm -r "addpath('/home/aleman/WORK/FAUST/faust/trunk/devcpp/Build-AL/wrapper/run_test/bin/');hier_fact;exit"

	
	# add_custom_command(OUTPUT ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp COMMAND cmake -E copy "${FAUST_TESTSRC_SRC_DIR}/hierarchical_fact_basic.cpp" "${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp" DEPENDS "${FAUST_TESTSRC_SRC_DIR}/hierarchical_fact_basic.cpp")
	
	 # add_executable(faust_hier_basic ${FAUST_TESTSRC_BIN_DIR}/hierarchical_fact_basic.cpp)
   # if(FAUST_USE_OPENBLAS)
      # target_link_libraries(faust_hier_basic  ${FAUST_TARGET} ${OPENBLAS_LIB_FILE})
   # else()
      # target_link_libraries(faust_hier_basic ${FAUST_TARGET})
   # endif()
   
   
   
   
