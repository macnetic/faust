
include_directories(${FAUST_MATRIX_SRC_DIR} ${FAUST_FAUSTCORE_SRC_DIR} ${FAUST_PALM4MSA_SRC_DIR} $ENV{EIGEN_DIR} ${FAUST_MATIO_SRC_DIR} $ENV{MATIODIR}/include)

# palm
add_executable(palm palm4MSA_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
target_link_libraries(palm ${FAUST_PALM4MSA_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})

# hier_fact
add_executable(faust_hier hierarchical_fact_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
target_link_libraries(faust_hier ${FAUST_PALM4MSA_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})

# hier_fact2
add_executable(faust_hier2 hierarchical_fact_test2.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
target_link_libraries(faust_hier2 ${FAUST_PALM4MSA_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})

# meg
add_executable(meg MEG_fact.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp)
target_link_libraries(meg ${FAUST_PALM4MSA_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})

# faust
add_executable(faust faust_test.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_mat.cpp ${FAUST_MATIO_SRC_DIR}/faust_init_from_matio_params.cpp)
target_link_libraries(faust ${FAUST_PALM4MSA_TARGET} ${MATIO_LIB_FILE} ${HDF5_LIB_FILE})

