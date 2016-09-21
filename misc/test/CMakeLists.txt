# WRAPPER : TESTS EXAMPLES USING THE LIBRARY FAUST WITH A CPU VERSION

message(STATUS "CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")

#${FAUST_EXCEPTION_SRC_DIR}
include_directories(${FAUST_SRC_LINEAR_OPERATOR_DIR} ${FAUST_LINEAR_OPERATOR_CPU_SRC_DIR} ${FAUST_ALGORITHM_CONSTRAINT_SRC_DIR} ${FAUST_ALGORITHM_FACTORIZATION_SRC_DIR} ${EIGEN_INC_DIR})

if(BUILD_READ_MAT_FILE)
	include_directories(${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR} ${MATIO_INC_DIR})
endif(BUILD_READ_MAT_FILE)

#rajout OPENBLAS_INC_DIR pour la compatibilite UBUNTU
if(BUILD_OPENBLAS)
	include_directories(${OPENBLAS_INC_DIR})
endif()


##############################################################################
# Testing with CTest
# CTest is a testing tool distributed as a part of CMake.
# Running the differents binary created in this project  
##############################################################################


# Critical time, which fixe the limit time of the process
set(TIMEOUT_MEG 8000)
set(TIMEOUT_COMP 100)


# MATLAB TEST
if(BUILD_MATLAB_MEX_FILES)

	if(WIN32 AND NOT UNIX)
		add_test(NAME MATLAB_FACT_HIER_MEX COMMAND matlab -nojvm -r "delete MATLAB_FACT_HIER_MEX.txt;diary MATLAB_FACT_HIER_MEX.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_compared_hierarchical_fact',9402,1,'MEX'); catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FACT_HIER_MEX.txt)
		add_test(NAME MATLAB_FACT_HIER_MATLAB COMMAND matlab -nojvm -r "delete MATLAB_FACT_HIER_MATLAB.txt;diary MATLAB_FACT_HIER_MATLAB.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_compared_hierarchical_fact',9402,1,'MATLAB'); catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FACT_HIER_MATLAB.txt)
		add_test(NAME MATLAB_FACT_MEG_MEX COMMAND matlab -nojvm -r "delete MATLAB_FACT_MEG_MEX.txt;diary MATLAB_FACT_MEG_MEX.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_MEG',22480,1,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FACT_MEG_MEX.txt)
		add_test(NAME MATLAB_FACT_MEG_MATLAB COMMAND matlab -nojvm -r "delete MATLAB_FACT_MEG_MATLAB.txt;diary MATLAB_FACT_MEG_MATLAB.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_MEG',22480,1,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FACT_MEG_MATLAB.txt)
		add_test(NAME MATLAB_FACT_HADAMARD_MEX COMMAND matlab -nojvm -r "delete MATLAB_FACT_HADAMARD_MEX.txt;diary MATLAB_FACT_HADAMARD_MEX.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_HADAMARD',5941,30,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)" ) #&& type MATLAB_FACT_HADAMARD_MEX.txt)
		add_test(NAME MATLAB_FACT_HADAMARD_MATLAB COMMAND matlab -nojvm -r "delete MATLAB_FACT_HADAMARD_MATLAB.txt;diary MATLAB_FACT_HADAMARD_MATLAB.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_HADAMARD',5941,30,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FACT_HADAMARD_MATLAB.txt)
		add_test(NAME MATLAB_DEMO_INSTALL COMMAND matlab -nodesktop -r "delete MATLAB_DEMO_INSTALL.txt;diary MATLAB_DEMO_INSTALL.txt;try;testpass=0;addpath('${CMAKE_INSTALL_MATLAB_PREFIX}');setup_FAUST;run_all_demo;catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_DEMO_INSTALL.txt)
		add_test(NAME MATLAB_DEMO_BUILD COMMAND matlab -nodesktop -r "delete MATLAB_DEMO_BUILD.txt;diary MATLAB_DEMO_BUILD.txt;try;testpass=0;addpath('${FAUST_MATLAB_BIN_DIR}');setup_FAUST;run_all_demo;catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)" ) #&& type MATLAB_DEMO_BUILD.txt)
		add_test(NAME MATLAB_FAUST_1 COMMAND matlab -nojvm -r "delete MATLAB_FAUST_1.txt;diary MATLAB_FAUST_1.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust;catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)" ) #&& type MATLAB_FAUST_1.txt)
		## test the matlab Faust class methods with incorrect parameter
		add_test(NAME MATLAB_FAUST_2 COMMAND matlab -nojvm -r "delete MATLAB_FAUST_2.txt;diary MATLAB_FAUST_2.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust2;catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FAUST_2.txt)
		add_test(NAME MATLAB_FAUST_SPEED COMMAND matlab -nojvm -r "delete MATLAB_FAUST_SPEED.txt;diary MATLAB_FAUST_SPEED.txt;try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust_time;catch ME ;testpass=-1;disp(getReport(ME)); end ; disp (testpass); diary OFF; exit(testpass)") # && type MATLAB_FAUST_SPEED.txt)
		
	else(WIN32 AND NOT UNIX)	
		add_test(NAME MATLAB_FACT_HIER_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_compared_hierarchical_fact',9402,1,'MEX'); catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)" )
		add_test(NAME MATLAB_FACT_HIER_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_compared_hierarchical_fact',9402,1,'MATLAB'); catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)" )
		add_test(NAME MATLAB_FACT_MEG_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_MEG',22480,1,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FACT_MEG_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_MEG',22480,1,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FACT_HADAMARD_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_HADAMARD',5941,30,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FACT_HADAMARD_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('config_HADAMARD',5941,30,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_DEMO_INSTALL COMMAND matlab -nodesktop -r "try;testpass=0;addpath('${CMAKE_INSTALL_MATLAB_PREFIX}');setup_FAUST;run_all_demo;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_DEMO_BUILD COMMAND matlab -nodesktop -r "try;testpass=0;addpath('${FAUST_MATLAB_BIN_DIR}');setup_FAUST;run_all_demo;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FAUST_1 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		# test the matlab Faust class methods with incorrect parameter
		add_test(NAME MATLAB_FAUST_2 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust2;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FAUST_SPEED COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust_time;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")

	endif(WIN32 AND NOT UNIX)
	# limit the time of the process in case of bug... 
	#CTEST_FULL_OUTPUT
	
	#set_tests_properties(MATLAB_FACT_HIER_MEX MATLAB_FACT_HIER_MATLAB   PROPERTIES TIMEOUT 100)
	#set_tests_properties(MATLAB_DEMO_INSTALL MATLAB_DEMO_BUILD PROPERTIES TIMEOUT 6000)
	#set_tests_properties(MATLAB_FAUST_1 MATLAB_FAUST_2  PROPERTIES TIMEOUT 30)
	#set_tests_properties(MATLAB_FACT_MEG_MATLAB MATLAB_FACT_MEG_MEX PROPERTIES TIMEOUT ${TIMEOUT_MEG})
endif()












####### IF THE MATIO LIBRARY IS AVAILABLE ELSE WE USE TXT OR XML DATA  #######
if(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE) # AND HDF5_LIB_FILE)

	### list of executable generated in double and simple precision
	# hierarchical_fact : test the hierarchical_fact algorithm
	# test_palm4MSA : test the palm4MSA algorithm
	# multiply_compare_time : time comparison between Faust-vector product and Dense matrix-vector product
 
	foreach(TEST_FPP float double)	
		foreach(testin hierarchicalFactorization test_palm4MSA Faust_multiplication)
			# copy CPU files cpp.in to the user's ./src/ directory in float et double precision
			configure_file(${FAUST_SRC_TEST_SRC_DIR}/${testin}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
			# Generation of the binary files in double and float precision
			add_executable(${testin}_${TEST_FPP} ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_BIN_TEST_SRC_DIR}/ ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio_mat.cpp)
   	
			if(BUILD_OPENBLAS)
		      		target_link_libraries(${testin}_${TEST_FPP}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
			else(BUILD_OPENBLAS)
		      		target_link_libraries(${testin}_${TEST_FPP}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
			endif(BUILD_OPENBLAS)	
		endforeach()
	endforeach()

#message(STATUS "ALALALLLLLLLLLLL FAUST_INSTALL_BIN=${FAUST_INSTALL_BIN}")

	install(TARGETS hierarchicalFactorization_double hierarchicalFactorization_float test_palm4MSA_double test_palm4MSA_float DESTINATION ${FAUST_INSTALL_BIN})


	# matlab script useful for displaying the result of the multiply_compare_time executable ( 
	configure_file(${FAUST_SRC_TEST_TIME_COMP_DIR}/display_multiply_comptime.m.in ${FAUST_BIN_TEST_TIME_COMPARE_DIR}/display_multiply_comptime.m @ONLY)
	install(FILES ${FAUST_BIN_TEST_TIME_COMPARE_DIR}/display_multiply_comptime.m DESTINATION ${FAUST_INSTALL_BIN_TIME_COMP})


###### IN CASE MATIO IS NOT AVAILABLE YOU CAN RUN THE TEST ######
else(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE)	
	#message(STATUS "Library matio is not available. Please used .xml input files.")
endif(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE)



if (BUILD_MATLAB_MEX_FILES)
	configure_file(${FAUST_SRC_TEST_TOOL_DIR}/set_path.m.in ${FAUST_BIN_TEST_TOOLS_DIR}/set_path.m @ONLY)
	
	
	foreach(matlabtest hier_fact_test.m test_matlab_faust.m test_matlab_faust2.m test_matlab_faust_time.m)
		configure_file(${FAUST_SRC_TEST_SRC_DIR}/${matlabtest} ${FAUST_BIN_TEST_BIN_DIR}/${matlabtest} @ONLY)
	endforeach()
endif(BUILD_MATLAB_MEX_FILES)







# MATLAB input format
if (BUILD_READ_MAT_FILE)
	
	# test for the algorithm (palm4MSA and hierarchical_fact) for different configuration (MEG, Hadamard,...) in simple and 	double precision
	
	add_test(NAME MATFILE_FACT_MEG_float COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_float ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22450 200)
	add_test(NAME MATFILE_FACT_MEG_double COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_double ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22480 1)
	set_tests_properties(MATFILE_FACT_MEG_float MATFILE_FACT_MEG_double PROPERTIES TIMEOUT ${TIMEOUT_MEG})

	foreach(TEST_FPP float double)
		add_test(NAME MATFILE_FACT_HADAMARD_${TEST_FPP} COMMAND  ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_${TEST_FPP} 	${FAUST_DATA_MAT_DIR}/config_HADAMARD.mat)
		add_test(NAME MATFILE_FACT_HIER_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_${TEST_FPP} 	${FAUST_DATA_MAT_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
	
		add_test(NAME MATFILE_MULTIPLICATION_MEG_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/Faust_multiplication_${TEST_FPP} ${FAUST_DEMO_BSL_DATA_SRC_DIR}X_meg.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_6.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_8.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_16.mat ${FAUST_DEMO_BSL_DATA_SRC_DIR}M_25.mat)

		add_test(NAME MATFILE_PALM4MSA_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/test_palm4MSA_${TEST_FPP} )

		# limit the time of the process in case of bug... 
		set_tests_properties(MATFILE_FACT_HADAMARD_${TEST_FPP} MATFILE_FACT_HIER_${TEST_FPP} PROPERTIES TIMEOUT ${TIMEOUT_COMP})
	
endforeach()


	

	
	


	## time comparison between Faust-vector product and Dense matrix-vector product for various configuration ( dimension, RCG , number of factor)
	#add_test(NAME MULTIPLY_COMPTIME COMMAND ${FAUST_BIN_TEST_BIN_DIR}/multiply_compare_time_double)
	## display the result of MULTIPLY_COMPTIME and save it in folder FAUST_BIN_TEST_FIG_DIR as jpeg image (WARNING must be launched after MULTIPLY_COMPTIME)
	#add_test(NAME FIG_MULTIPLY_COMPTIME COMMAND	matlab -nodesktop -r "try;addpath('${FAUST_BIN_TEST_TIME_COMPARE_DIR}');display_multiply_comptime('${FAUST_BIN_TEST_OUTPUT_DIR}/multiply_compare_time.mat','${FAUST_BIN_TEST_FIG_DIR}') ; catch ME ; disp('ERROR IN TEST display_multiply_comptime.m Please chek the output to visualize the error. ') ; end ; exit")
	
	#set_tests_properties(FIG_MULTIPLY_COMPTIME PROPERTIES TIMEOUT 300)

endif (BUILD_READ_MAT_FILE)


# XML input format
if (BUILD_READ_XML_FILE)
	add_test(NAME XML_FACT_HIER COMMAND ${FAUST_BIN_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CONFIG_XML_DIR}/hier_fact_config.xml ${FAUST_DATA_TXT_DIR}/hier_fact_data.txt )

	add_test(NAME XML_FACT_HIER_TRANSPOSE COMMAND ${FAUST_BIN_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CONFIG_XML_DIR}/trans_hier_fact_config.xml ${FAUST_DATA_TXT_DIR}/hier_fact_data.txt  T)

	add_test(NAME XML_FACT_MEG COMMAND ${FAUST_BIN_CMDLINE_BIN_DIR}/launch_hierarchical_fact ${FAUST_CONFIG_XML_DIR}/MEG_config.xml ${FAUST_DATA_TXT_DIR}/MEG_data.txt )
	# limit the time of the process in case of bug... 
	set_tests_properties(XML_FACT_HIER XML_FACT_HIER_TRANSPOSE PROPERTIES TIMEOUT ${TIMEOUT_COMP})
	set_tests_properties(XML_FACT_MEG PROPERTIES TIMEOUT ${TIMEOUT_MEG})	 
endif(BUILD_READ_XML_FILE)



 



	



########################################################################
##################  TEST VERSION GPU (version using CUDA)
#if (BUILD_USE_GPU)
#	#${FAUST_EXCEPTION_SRC_DIR} 
#	include_directories(${FAUST_MATRIX_CU_SRC_DIR} ${FAUST_FAUSTCORE_CU_SRC_DIR} ${FAUST_ALGORITHM_CU_SRC_DIR} ${CUDA_INC_DIR} ${CUDA_RUNTIME_API_INC_DIR} ${CUSPARSE_INC_DIR} ${CUBLAS_V2_INC_DIR}  ${EIGEN_INC_DIR})

#	if(BUILD_READ_MAT_FILE)
#		include_directories(${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR} ${MATIO_INC_DIR})
#	endif(BUILD_READ_MAT_FILE)


#	add_definitions(-DCOMPILE_GPU)
#	# generation des executables de algorithme hierarchical_fact en GPU en simple et double precision 
#	foreach(TEST_FPP float double)	
#			foreach(testin hierarchical_fact_cu)
#				# copy CPU files cpp.in to the user's ./src/ directory en float et double precision
#				configure_file(${FAUST_SRC_TEST_SRC_DIR}/${testin}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp @ONLY)
#				# Creation des executable en double et en float
#				add_executable(${testin}_${TEST_FPP} ${FAUST_BIN_TEST_SRC_DIR}/${testin}_${TEST_FPP}.cpp ${FAUST_BIN_TEST_SRC_DIR}/ ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio_mat.cpp)
#				target_link_libraries(${testin}_${TEST_FPP} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${CUBLAS_LIB_FILE} ${CUDART_LIB_FILE}  ${CUSPARSE_LIB_FILE})

#				install(TARGETS ${testin}_${TEST_FPP}  DESTINATION ${FAUST_INSTALL_BIN})	
#		endforeach()
#	endforeach()

#	# test GPU
#	add_test(NAME HADAMARD_GPU_FACT_DOUBLE COMMAND  ${FAUST_BIN_TEST_BIN_DIR}/hierarchical_fact_cu_double 	${FAUST_DATA_MAT_DIR}/config_HADAMARD.mat)
#	add_test(NAME FAUST_HIER_GPU_DOUBLE COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchical_fact_cu_double 	${FAUST_DATA_MAT_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
#	add_test(NAME FAUST_HIER_GPU_FLOAT COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchical_fact_cu_float ${FAUST_DATA_MAT_DIR}/config_compared_hierarchical_fact.mat 9401.5 0.1)
#	add_test(NAME MEG_GPU_FACT_FLOAT COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchical_fact_cu_float ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22332 1)
#	add_test(NAME MEG_GPU_FACT_DOUBLE COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchical_fact_cu_double ${FAUST_DATA_MAT_DIR}/config_MEG.mat 22480 1)
#	set(TIMEOUT_MEG 14000)
#	set_tests_properties(MEG_GPU_FACT_FLOAT MEG_GPU_FACT_DOUBLE PROPERTIES TIMEOUT ${TIMEOUT_MEG})
#endif(BUILD_USE_GPU)
########################################################################


