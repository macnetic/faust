#ifndef GIVENS_TEST_UTILITY_@GIVENS_SCALAR@
#define GIVENS_TEST_UTILITY_@GIVENS_SCALAR@
/**
 * Tests iteration errors by comparing matlab and C++ impl. errors.
 */
template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_ite_errors(const Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file, int ref_ite_period = 0);

/**
 * Tests eigenvalues getting them in unordered or ordered sequence, comparing them to matlab ref. ones.
 */
template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_eigenvalues(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file);

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_eigentransform(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file);

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_err_against_Laplacian(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file);

/**
 * Tests that pivot choices made by the C++ impl. are sufficiently close to the ones made by matlab ref.
 */
template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_pivot_choices(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file, const float same_pivot_target_rate=.99f, const int stop_count_ite = 57);

#include "EigTJUtil@GIVENS_SCALAR@.hpp"
#endif
