foreach(REAL_TYPE IN ITEMS float double)
  if(REAL_TYPE MATCHES double)
    set(SVD_COMPLEX 1)
  else()
    set(SVD_COMPLEX 0)
  endif()
  configure_file(${FAUST_ALGORITHM_SRC_DIR}/factorization/faust_SVDTJ.cpp.in ${FAUST_ALGORITHM_SRC_DIR}/factorization/faust_SVDTJ_${REAL_TYPE}.cpp)
endforeach()

