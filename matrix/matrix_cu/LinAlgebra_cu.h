#ifndef __LINALGEBRA_CU_H__
#define __LINALGEBRA_CU_H__

#include "faust_constant.h"
#ifdef __COMPILE_SPMAT__
   #include "cusparse.h"
  template <typename faust_real> class faust_cu_spmat;
#endif
#include "cublas_v2.h"

template <typename faust_real> class faust_cu_mat;
template <typename faust_real> class faust_cu_vec;



//----------------------------------------------//
//------------ FUNCTION DECLARATIONS -----------//
//----------------------------------------------//


////// FUNCTIONS with faust_cu_mat matrices //////

// opA1=nb_rows(cu_A^opA) ; opA2=nb_cols(cu_A^opA)
template <typename faust_real>
void setOp(const faust_cu_mat<faust_real>& cu_A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);

// C = A + B;
template <typename faust_real>
void add(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C);

// compute the biggest eigenvalue of A, A must be semi-definite positive 
template<typename T>	
T power_iteration(const  faust_cu_mat<T> & cu_A, const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag, cublasHandle_t cublasHandle);


// y = alpha*op(A)*x + beta*y
template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const faust_real alpha, const faust_real beta, char typeA, cublasHandle_t cublasHandle);
// y = op(A) * x 
template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, char typeA, cublasHandle_t cublasHandle);
// y = A * x 
template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cublasHandle_t cublasHandle);
// y = A * x 
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cublasHandle_t cublasHandle); 

// C = alpha*op(A)*op(B) + beta*C
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C,const faust_real alpha, const faust_real beta, char  typeA, char  typeB, cublasHandle_t cublasHandle);
// C = op(A) * op(B) 
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, char typeA, char typeB, cublasHandle_t cublasHandle);
// C = A * B 
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle);
// C = A * B 
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle); 

//////////////////////////////////////////////////



#ifdef __COMPILE_SPMAT__
///// FUNCTIONS with faust_cu_spmat matrices /////

// opA1=nb_rows(cu_S^opA) ; opA2=nb_cols(cu_S^opA)
template <typename faust_real>
void setOp(const faust_cu_spmat<faust_real>& cu_S, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);

// y = alpha*op(A)*x + beta*y
template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const faust_real alpha, const faust_real beta, const char opA, cusparseHandle_t cusparseHandle);
// y = op(A) * x 
template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const char opA, cusparseHandle_t cusparseHandle);
// y = A * x 
template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cusparseHandle_t cusparseHandle);
// y = A * x 
template <typename faust_real>
void multiply(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle); 

// C = alpha*op(A)*op(B) + beta*C ; with A sparse
template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const  faust_real alpha, const faust_real beta, const char opA, const char opB, cusparseHandle_t cusparseHandle);
// C = alpha*op(A)*op(B) + beta*C ; with B sparse
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const  faust_real alpha, const faust_real beta, const char opA, const char opB, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle);
// C = op(A)*op(B) ; with A sparse
template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const char opA, const char opB, cusparseHandle_t cusparseHandle);
// C = op(A)*op(B) ; with B sparse
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const char opA, const char opB, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle);
// C = A*B ; with A sparse
template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cusparseHandle_t cusparseHandle);
// C = A*B ; with B sparse
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle);
// C = A*B ; with A sparse
template <typename faust_real>
void multiply(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cusparseHandle_t cusparseHandle); 
// C = A*B ; with B sparse
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle); 

//////////////////////////////////////////////////
#endif // __COMPILE_SPMAT__








//----------------------------------------------//
//------------- INLINE DEFINITIONS -------------//
//----------------------------------------------//

// because of header include issues, the definition of inline functions have been moved to the end LinAlgebra_cu.hpp
//

#include "LinAlgebra_cu.hpp"

#endif
