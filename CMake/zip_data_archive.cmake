set(DATA_ARCHIVE_FILES;config_compared_palm2.mat;faust_MEG_rcg_16.mat;faust_MEG_rcg_25.mat;faust_MEG_rcg_6.mat;faust_MEG_rcg_8.mat;faust_quick_start.mat;F_DYNPROG.mat;F_GREEDY.mat;HierarchicalFactFFT_test_U_L_params.mat;Laplacian_1024_community.mat;Laplacian_1024_erdos_renyi.mat;Laplacian_1024_path.mat;Laplacian_1024_random_ring.mat;Laplacian_1024_ring.mat;Laplacian_1024_sensor.mat;Laplacian_128_community.mat;Laplacian_128_erdos_renyi.mat;Laplacian_128_path.mat;Laplacian_128_random_ring.mat;Laplacian_128_ring.mat;Laplacian_128_sensor.mat;Laplacian_256_community.mat;Laplacian_256_erdos_renyi.mat;Laplacian_256_path.mat;Laplacian_256_random_ring.mat;Laplacian_256_ring.mat;Laplacian_256_sensor.mat;Laplacian_512_community.mat;Laplacian_512_erdos_renyi.mat;Laplacian_512_path.mat;Laplacian_512_random_ring.mat;Laplacian_512_ring.mat;Laplacian_512_sensor.mat;matrix_HADAMARD_32.mat;matrix_hierarchical_fact.mat;matrix_MEG.mat;M_DYNPROG.mat;M_GREEDY.mat;ref_test_PALM4SMA_FFT2.mat;supports.mat;test_GivensDiag_Lap_U_J.mat;test_GivensDiagParallel_Lap_U_J_choices.mat)
add_custom_target(zip_data_archive zip faust_data.zip -O ${PROJECT_BINARY_DIR}/doc/html/faust_data.zip ${DATA_ARCHIVE_FILES} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/misc/data/mat COMMENT "Archiving FAµST data into faust_data.zip")
