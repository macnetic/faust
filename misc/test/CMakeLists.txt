
if(BUILD_TESTING)
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
set(TIMEOUT_MEG_TRANSPOSED 16000)
set(TIMEOUT_COMP 100)


# MATLAB TEST
if(BUILD_WRAPPER_MATLAB)

	if(WIN32 AND NOT UNIX)
		#WINDOWS :
		# all the matlab test have the same structure, they launch a batch file  ${FAUST_SRC_TEST_TOOL_DIR}/test_matlab.bat
		# This script takes 2 arguments :
		# 	1st arg : the list of matlab script, function that will be executed
		#		WARNING : due to batch file issue,
		#		in the Matlab command passed as a string to the batch file,  
		#		"@" means "," and "$" means ";"
		#   2nd arg : the ouputfile that stores all the things that will be display in the matlab command window
		
		add_test(NAME MATLAB_FACT_HIER_MEX COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_hierarchical_fact'@'config_hierarchical_fact'@9402@1@'MEX')$" "${FAUST_BIN_TEST_OUTPUT_DIR}/MATLAB_FACT_HIER_MEX") 
		add_test(NAME MATLAB_FACT_HIER_MATLAB COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_hierarchical_fact'@'config_hierarchical_fact'@9402@1@'MATLAB')$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FACT_HIER_MATLAB")
		
		
		if(${SLOW_TESTS}) 
			add_test(NAME MATLAB_FACT_MEG_MEX COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_MEG'@'config_MEG'@22480@1@'MEX')$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FACT_MEG_MEX")
			add_test(NAME MATLAB_FACT_MEG_MATLAB COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_MEG'@'config_MEG'@22480@1@'MATLAB')$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FACT_MEG_MATLAB")
		endif()
		add_test(NAME MATLAB_FACT_HADAMARD_MEX COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_HADAMARD_32'@'config_HADAMARD_32'@5941@30@'MEX')$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FACT_HADAMARD_MEX")
		add_test(NAME MATLAB_FACT_HADAMARD_MATLAB COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$hier_fact_test('matrix_HADAMARD_32'@'config_HADAMARD_32'@5941@30@'MATLAB')$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FACT_HADAMARD_MATLAB")

		if(${SLOW_TESTS})
			add_test(NAME MATLAB_DEMO_INSTALL COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${CMAKE_INSTALL_MATLAB_PREFIX}')$setup_FAUST$import matfaust.demo.runall$runall$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_DEMO_INSTALL")
			add_test(NAME MATLAB_DEMO_BUILD COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_MATLAB_BIN_DIR}')$setup_FAUST$import matfaust.demo.runall$runall$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_DEMO_BUILD")
		endif()
		add_test(NAME MATLAB_FAUST_CONFIG1 COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$test_matlab_faust(5@4@10@3)$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FAUST_1")
		add_test(NAME MATLAB_FAUST_CONFIG2 COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$test_matlab_faust(2@15@1@1)$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FAUST_1")
		# test the matlab Faust class methods with incorrect parameter
		add_test(NAME MATLAB_FAUST_2 COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$test_matlab_faust2$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FAUST_2")
		add_test(NAME MATLAB_FAUST_SPEED COMMAND ${FAUST_SRC_TEST_TOOL_DIR}/windows_test_matlab.bat "addpath('${FAUST_BIN_TEST_BIN_DIR}'@'${FAUST_BIN_TEST_TOOLS_DIR}')$set_path$test_matlab_faust_time$" "${FAUST_BIN_TEST_BIN_DIR}/MATLAB_FAUST_SPEED")
		
	
	
	
	else(WIN32 AND NOT UNIX)	
		add_test(NAME MATLAB_FACT_HIER_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_hierarchical_fact','config_hierarchical_fact',9402,1,'MEX'); catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)" )
	
		
		add_test(NAME MATLAB_FACT_HIER_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_hierarchical_fact','config_hierarchical_fact',9402,1,'MATLAB'); catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)" )
	
		if(${SLOW_TESTS})
			add_test(NAME MATLAB_FACT_MEG_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_MEG','config_MEG',22480,1,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
			add_test(NAME MATLAB_FACT_MEG_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_MEG','config_MEG',22480,1,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		endif()

		add_test(NAME MATLAB_FACT_HADAMARD_MEX COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_HADAMARD_32','config_HADAMARD_32',5941,30,'MEX');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		add_test(NAME MATLAB_FACT_HADAMARD_MATLAB COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;hier_fact_test('matrix_HADAMARD_32','config_HADAMARD_32',5941,30,'MATLAB');catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		
		if(${SLOW_TESTS})
			add_test(NAME MATLAB_DEMO_INSTALL COMMAND matlab -nodesktop -r "try;testpass=0;addpath('${CMAKE_INSTALL_MATLAB_PREFIX}');setup_FAUST;import matfaust.demo.runall;runall;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
			add_test(NAME MATLAB_DEMO_BUILD COMMAND matlab -nodesktop -r "try;testpass=0;addpath('${FAUST_MATLAB_BIN_DIR}');setup_FAUST;import matfaust.demo.runall;runall;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		endif(${SLOW_TESTS})
		
		#test Faust with 3 real factors : 
		# 1st factor sparse size 10x3, 
		# 2nd factor dense size 3x2, 
		# 3rd factor sparse size 2x5		
		add_test(NAME MATLAB_FAUST_CONFIG_REAL1 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([10,3,2,5],{'real.sparse','real.dense','real.sparse'},0.5);test_matlab_faust(factors,prod,30,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")	
		
		#test Faust with 3 real factors : 
		# 1st factor sparse size 5x2, 
		# 2nd factor dense size 2x3, 
		# 3rd factor sparse size 3x10		
		add_test(NAME MATLAB_FAUST_CONFIG_REAL2 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([5,2,3,10],{'real.dense','real.sparse','real.dense'},0.5);test_matlab_faust(factors,prod,30,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		#test Faust with 1 real factor : 
		# 1st factor dense size 10x3, 		
		add_test(NAME MATLAB_FAUST_CONFIG_REAL3 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([100,10],{'real.dense'},1.0);test_matlab_faust(factors,prod,30,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")		
		
		#test Faust with 1 real factor : 
		# 1st factor sparse size 10x3, 		
		add_test(NAME MATLAB_FAUST_CONFIG_REAL4 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([100,10],{'real.sparse'},0.1);test_matlab_faust(factors,prod,1000,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")



		#test Faust with 3 factors : 
		# 1st factor real dense size 10x3, 
		# 2nd factor complex dense size 3x2, 
		# 3rd factor complex dense size 2x5		
		add_test(NAME MATLAB_FAUST_CONFIG_CPLX1 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([10,3,2,5],{'real.dense','complex.dense','complex.dense'},0.5);test_matlab_faust(factors,prod,30,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")
		
		#test Faust with 3 factors : 
		# 1st factor real sparse size 10x3, 
		# 2nd factor complex dense size 3x2, 
		# 3rd factor complex sparse size 2x5		
		add_test(NAME MATLAB_FAUST_CONFIG_CPLX2 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;[factors,prod]=generate_Faust_config([10,3,2,5],{'real.sparse','complex.dense','complex.sparse'},0.5);test_matlab_faust(factors,prod,30,0);catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")


		## test the matlab Faust class methods with incorrect parameter
		add_test(NAME MATLAB_FAUST_2 COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust2;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")

		## make some speed-up comparison
		add_test(NAME MATLAB_FAUST_SPEED COMMAND matlab -nojvm -r "try;testpass=0;addpath('${FAUST_BIN_TEST_BIN_DIR}','${FAUST_BIN_TEST_TOOLS_DIR}');set_path;test_matlab_faust_time;catch ME ;testpass=-1;disp(getReport(ME)); end ; exit(testpass)")


		add_test(NAME MATLAB_FAUST_UNIT_TESTS COMMAND matlab -nodisplay -nojvm -r "try;testpass=0;addpath('${PROJECT_BINARY_DIR}/wrapper/matlab:${PROJECT_BINARY_DIR}/run_test/tools:${PROJECT_BINARY_DIR}/run_test/bin');set_path;cd ${FAUST_SRC_TEST_SRC_MATLAB_DIR};FaustTest ${PROJECT_BINARY_DIR}/wrapper/matlab/ ${PROJECT_BINARY_DIR}/run_test/bin ${PROJECT_BINARY_DIR}/run_test/tools;catch ME;testpass=-1;displ(getReport(ME));end;exit(testpass)")

		add_test(NAME MATLAB_FAUST_FACTORY_UNIT_TESTS COMMAND matlab -nodisplay -nojvm -r "try;testpass=0;addpath('${PROJECT_BINARY_DIR}/wrapper/matlab:${PROJECT_BINARY_DIR}/run_test/tools:${PROJECT_BINARY_DIR}/run_test/bin');set_path;cd ${FAUST_SRC_TEST_SRC_MATLAB_DIR};FaustFactoryTest ${PROJECT_BINARY_DIR}/wrapper/matlab/ ${PROJECT_BINARY_DIR}/run_test/bin ${PROJECT_BINARY_DIR}/run_test/tools;catch ME;testpass=-1;displ(getReport(ME));end;exit(testpass)")


	#FIN MODIF NB VM MAC XCode

	endif(WIN32 AND NOT UNIX)
	# limit the time of the process in case of bug... 
	#CTEST_FULL_OUTPUT
	
	#DEBUT MODIF NB VM MAC XCode
	set_tests_properties(MATLAB_FACT_HIER_MEX MATLAB_FACT_HIER_MATLAB   PROPERTIES TIMEOUT 100)
	if(${SLOW_TESTS})
		set_tests_properties(MATLAB_DEMO_INSTALL MATLAB_DEMO_BUILD PROPERTIES TIMEOUT 6000)
	endif()
	if(NOT WIN32)
		set_tests_properties(MATLAB_FAUST_CONFIG_REAL1 MATLAB_FAUST_CONFIG_REAL2 MATLAB_FAUST_CONFIG_REAL3 MATLAB_FAUST_CONFIG_REAL4 MATLAB_FAUST_CONFIG_CPLX1 MATLAB_FAUST_CONFIG_CPLX2 MATLAB_FAUST_2  PROPERTIES TIMEOUT 30)
	endif()
	if(${SLOW_TESTS})
		set_tests_properties(MATLAB_FACT_MEG_MATLAB MATLAB_FACT_MEG_MEX PROPERTIES TIMEOUT ${TIMEOUT_MEG})
	endif()
	#FIN MODIF NB VM MAC XCode
endif()


####### UNIT TEST ########



foreach(TEST_FPP float double)
	foreach(FILE faust_mult faust_mult_cplx test_Vect_min test_MatDense_min test_MatDense_get_row test_MatDense_lower_upper_tri test_MatDense_nonzeros_indices test_Transform_move)
		set(TEST_BIN_FILE ${FILE}_${TEST_FPP})
		set(TEST_FILE_CPP ${TEST_BIN_FILE}.cpp)
		message(STATUS ${TEST_FILE_CPP})
		configure_file(${FAUST_SRC_TEST_SRC_UNIT_DIR}/${FILE}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${TEST_FILE_CPP} @ONLY)

	# Generation of the binary files in double and float precision
	add_executable(${TEST_BIN_FILE} ${FAUST_BIN_TEST_SRC_DIR}/${TEST_FILE_CPP} ${FAUST_BIN_TEST_SRC_DIR}/)
   	
	target_link_libraries(${TEST_BIN_FILE} ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
			
	add_test(NAME UNIT_${TEST_BIN_FILE} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/${TEST_BIN_FILE})
				
	endforeach(FILE)
endforeach(TEST_FPP)
	



####### FIN UNIT TEST #####







####### IF THE MATIO LIBRARY IS AVAILABLE ELSE WE USE TXT OR XML DATA  #######
if(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE) # AND HDF5_LIB_FILE)

	### list of executable generated in double and simple precision
	# hierarchical_fact : test the hierarchical_fact algorithm
	# test_palm4MSA : test the palm4MSA algorithm
	# faust_multiplication : time comparison between Faust-vector product and Dense matrix-vector product
 
	foreach(TEST_FPP float double complex<float> complex<double>)	
		foreach(testin hierarchicalFactorization hierarchicalFactorizationFFT test_palm4MSA test_palm4MSAFFT faust_multiplication faust_matdense_conjugate GivensFGFT GivensFGFTParallel)
		if(${TEST_FPP} MATCHES complex AND ${testin} MATCHES GivensFGFT)
			continue()
			# GivensFGFT doesn't handle complex matrices
		endif()
		string(REGEX REPLACE "<|>" "" TEST_FILE_CPP ${TEST_FPP})
		set(TEST_CPP_NAME ${testin}_${TEST_FILE_CPP})
		set(TEST_FILE_CPP ${TEST_CPP_NAME}.cpp)
		message(STATUS ${TEST_FILE_CPP})			# copy CPU files cpp.in to the user's ./src/ directory in float et double precision
		if(TEST_FPP MATCHES float)
			set(TEST_FPP2 "float")
		else()
			set(TEST_FPP2 "double")
		endif()
		configure_file(${FAUST_SRC_TEST_SRC_CPP_DIR}/${testin}.cpp.in ${FAUST_BIN_TEST_SRC_DIR}/${TEST_FILE_CPP} @ONLY)
			# Generation of the binary files in double and float precision
			add_executable(${TEST_CPP_NAME} ${FAUST_BIN_TEST_SRC_DIR}/${TEST_FILE_CPP} ${FAUST_BIN_TEST_SRC_DIR}/ ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_CMDLINE_TYPE_FORMAT_MAT_SRC_DIR}/faust_init_from_matio_mat.cpp)
   	
			if(BUILD_OPENBLAS)
				target_link_libraries(${TEST_CPP_NAME}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE} ${OPENBLAS_LIB_FILE})
			else(BUILD_OPENBLAS)
				target_link_libraries(${TEST_CPP_NAME}  ${FAUST_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})
			endif(BUILD_OPENBLAS)	
		endforeach()
	endforeach()






###### IN CASE MATIO IS NOT AVAILABLE YOU CAN RUN THE TEST ######
else(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE)	
	#message(STATUS "Library matio is not available. Please used .xml input files.")
endif(MATIO_LIB_FILE AND MATIO_INC_DIR AND BUILD_READ_MAT_FILE)



if (BUILD_WRAPPER_MATLAB)
	configure_file(${FAUST_SRC_TEST_TOOL_DIR}/set_path.m.in ${FAUST_BIN_TEST_TOOLS_DIR}/set_path.m @ONLY)
	
	
	foreach(matlabtest hier_fact_test.m generate_Faust_config.m test_matlab_faust.m test_matlab_faust2.m test_matlab_faust_mult.m test_matlab_faust_time.m)
		configure_file(${FAUST_SRC_TEST_SRC_MATLAB_DIR}/${matlabtest} ${FAUST_BIN_TEST_BIN_DIR}/${matlabtest} @ONLY)
	endforeach()
endif(BUILD_WRAPPER_MATLAB)




#PYTHON WRAPPER
if (BUILD_WRAPPER_PYTHON)
	foreach(PY_VER RANGE 2 3)
		math(EXPR PY_INDEX "${PY_VER}-2")
		list(LENGTH PYTHON_EXES NUM_PY_EXES)
		if(${NUM_PY_EXES} GREATER ${PY_INDEX})
			list(GET PYTHON_EXES ${PY_INDEX} PYTHON_EXE)
			add_test(NAME PYTHON${PY_VER}_FAUST COMMAND ${PYTHON_EXE} ${FAUST_SRC_TEST_SRC_PYTHON_DIR}/test_pyFaust.py ${FAUST_PYTHON_BIN_DIR})


			if(PYTHON_MODULE_SCIPY)
				add_test(NAME PYTHON${PY_VER}_FAUST_TIME COMMAND ${PYTHON_EXE} ${FAUST_SRC_TEST_SRC_PYTHON_DIR}/test_pyFaust_time.py ${FAUST_PYTHON_BIN_DIR} ${FAUST_BIN_TEST_FIG_DIR})
				add_test(NAME PYTHON${PY_VER}_FAUST_DEMO_INSTALL COMMAND ${PYTHON_EXE} -c "import sys; sys.path += ['${FAUST_PYTHON_BIN_DIR}'];from pyfaust.demo import quickstart; quickstart.run()")
				add_test(NAME PYTHON${PY_VER}_FAUST_UNIT_TESTS COMMAND ${PYTHON_EXE} ${FAUST_SRC_TEST_SRC_PYTHON_DIR}/test_FaustPy.py ${FAUST_PYTHON_BIN_DIR}) # second arg. is the FaustPy's dir. to add to PYTHONPATH
			endif(PYTHON_MODULE_SCIPY)
		endif()
	endforeach()	
endif(BUILD_WRAPPER_PYTHON)


# MATLAB input format
if (BUILD_READ_MAT_FILE)
	
	foreach(TEST_FPP float double)
		add_test(NAME MATFILE_FACT_HADAMARD_${TEST_FPP} COMMAND  ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_HADAMARD_32.mat	${FAUST_CONFIG_MAT_DIR}/config_HADAMARD_32.mat 5930 30)
		add_test(NAME MATFILE_FACT_HIER_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_hierarchical_fact.mat	${FAUST_CONFIG_MAT_DIR}/config_hierarchical_fact.mat 9401.5 0.1)
	
		add_test(NAME MATFILE_MULTIPLICATION_MEG_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat T ${FAUST_DATA_MAT_DIR}/faust_MEG_rcg_6.mat ${FAUST_DATA_MAT_DIR}/faust_MEG_rcg_8.mat ${FAUST_DATA_MAT_DIR}/faust_MEG_rcg_16.mat ${FAUST_DATA_MAT_DIR}/faust_MEG_rcg_25.mat)

		if (EXISTS ${FAUST_DATA_MAT_DIR}/visage/)
			# data equipe visage test
			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G1_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G1.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G1_RCG_10_factor_3.mat)
	
			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G2_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G2.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G2_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G3_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G3.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G3_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G4_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G4.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G4_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G5_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G5.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G5_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G6_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G6.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G6_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G7_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G7.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G7_RCG_10_factor_3.mat)

			add_test(NAME MATFILE_MULTIPLICATION_VISAGE_G8_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/faust_multiplication_${TEST_FPP} ${FAUST_DATA_MAT_DIR}/visage/matrix_G8.mat N ${FAUST_DATA_MAT_DIR}/visage/faust_G8_RCG_10_factor_3.mat)

		else ()
			message (STATUS "Warning : The data for demo Visage are not available. ")		
		endif()

		add_test(NAME MATFILE_PALM4MSA_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/test_palm4MSA_${TEST_FPP} )
		add_test(NAME MATFILE_PALM4MSA_FFT_${TEST_FPP} COMMAND ${FAUST_BIN_TEST_BIN_DIR}/test_palm4MSAFFT_${TEST_FPP} )

		# limit the time of the process in case of bug... 
		set_tests_properties(MATFILE_FACT_HADAMARD_${TEST_FPP} MATFILE_FACT_HIER_${TEST_FPP} PROPERTIES TIMEOUT ${TIMEOUT_COMP})

	endforeach()

	if(${SLOW_TESTS})	
		# test for the algorithm (palm4MSA and hierarchical_fact) for different configuration (MEG, Hadamard,...) in simple and 	double precision
		
		add_test(NAME MATFILE_FACT_MEG_float COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_float ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG.mat 22450 500)
		add_test(NAME MATFILE_FACT_MEG_double COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_double ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG.mat 22480 1)
		
		add_test(NAME MATFILE_FACT_MEG_TRANSPOSED_float COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_float ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG_transposed.mat 22450 400 T)
		add_test(NAME MATFILE_FACT_MEG_TRANSPOSED_double COMMAND ${FAUST_BIN_TEST_BIN_DIR}/hierarchicalFactorization_double ${FAUST_DATA_MAT_DIR}/matrix_MEG.mat ${FAUST_CONFIG_MAT_DIR}/config_MEG_transposed.mat 22480 1 T)

		set_tests_properties(MATFILE_FACT_MEG_float MATFILE_FACT_MEG_double PROPERTIES TIMEOUT ${TIMEOUT_MEG})
		set_tests_properties(MATFILE_FACT_MEG_TRANSPOSED_float MATFILE_FACT_MEG_TRANSPOSED_double PROPERTIES TIMEOUT ${TIMEOUT_MEG_TRANSPOSED})
		
	endif(${SLOW_TESTS})

	


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

endif(BUILD_TESTING)

 
