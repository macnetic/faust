#ifndef __LINALGEBRA_GPU_H__
#define __LINALGEBRA_GPU_H__


//modif ALAL
#define __COMPILE_SPMAT__

#include "faust_constant_gpu.h"

#ifdef __COMPILE_SPMAT__
   #include "cusparse.h"
   #include "faust_SpBlasHandle_gpu.h"
#endif
#include "cublas_v2.h"
#include "faust_Vect_gpu.h"



#include <iostream>
#include "faust_MatDense_gpu.h"
#include "faust_cuda.h"
#include "kernels.h"
#ifdef __COMPILE_SPMAT__
    #include "faust_MatSparse_gpu.h"
#endif
#include "faust_exception.h"
#include "faust_BlasHandle_gpu.h"
#ifdef __COMPILE_TIMERS__
    #include "faust_Timer.h"
#endif

namespace Faust
{


    template<typename FPP,Device DEVICE> class MatDense;
    template<Device DEVICE> class BlasHandle;

    template<typename FPP,Device DEVICE> class Vect;



    //----------------------------------------------//
    //------------ FUNCTION DECLARATIONS -----------//
    //----------------------------------------------//
    // opA1=nb_rows(cu_A^opA) ; opA2=nb_cols(cu_A^opA)
    template <typename FPP>
    void setOp(const Faust::MatDense<FPP,Gpu>& cu_A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);

    //! \fn add
    //! \brief C = A + B;
    template <typename FPP>
    void add(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C);

    //! \fn power_iteration
    //! \brief compute the biggest eigenvalue of A, A must be semi-definite positive
    template<typename FPP>
    FPP power_iteration(const  Faust::MatDense<FPP,Gpu> & cu_A, const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void Faust::gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const FPP alpha, const FPP beta, char typeA, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs - y = alpha*op(A)*x + beta*y
    template <typename FPP>
    void gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const FPP alpha, const FPP beta, char typeA, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, char typeA, Faust::BlasHandle<Gpu> blasHandle);
    //! y = op(A) * x
    template <typename FPP>
    void gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, char typeA, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs - y = A * x
    template <typename FPP>
    void gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs y = A * x
    template <typename FPP>
    void multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C,const FPP alpha, const FPP beta, char  typeA, char  typeB, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs C = alpha*op(A)*op(B) + beta*C
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C,const FPP alpha, const FPP beta, char  typeA, char  typeB, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, char typeA, char typeB, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs - C = op(A) * op(B)
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, char typeA, char typeB, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle);
    //! Performs C = A * B
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle);

    //! \fn multiply
    //! Performs C = A * B
    template <typename FPP>
    void multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle);

    //////////////////////////////////////////////////

#ifdef __COMPILE_SPMAT__
    ///// FUNCTIONS with faust_cu_spmat matrices /////
    template<typename FPP,Device DEVICE> class MatSparse;
    template<Device DEVICE> class SpBlasHandle;

    //! \fn setOp
    //! Performs opA1=nb_rows(cu_S^opA) ; opA2=nb_cols(cu_S^opA)
    template <typename FPP>
    void setOp(const Faust::MatSparse<FPP,Gpu>& cu_S, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2);

    //! \fn gemv
    //! Performs y = alpha*op(A)*x + beta*y
    template <typename FPP>
    void gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const FPP alpha, const FPP beta, const char opA, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemv
    //! Performs y = op(A) * x
    template <typename FPP>
    void gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const char opA, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemv
    //! Performs y = A * x
    template <typename FPP>
    void gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn multiply
    //! Performs y = A * x
    template <typename FPP>
    void multiply(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = alpha*op(A)*op(B) + beta*C ; with A sparse
    template <typename FPP>
    void gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const  FPP alpha, const FPP beta, const char opA, const char opB, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = alpha*op(A)*op(B) + beta*C ; with B sparse
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const  FPP alpha, const FPP beta, const char opA, const char opB, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = op(A)*op(B) ; with A sparse
    template <typename FPP>
    void gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const char opA, const char opB, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = op(A)*op(B) ; with B sparse
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const char opA, const char opB, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = A*B ; with A sparse
    template <typename FPP>
    void gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn gemm
    //! Performs C = A*B ; with B sparse
    template <typename FPP>
    void gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn multiply
    //! Performs C = A*B ; with A sparse
    template <typename FPP>
    void multiply(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::SpBlasHandle<Gpu> spblasHandle);

    //! \fn multiply
    //! Performs C = A*B ; with B sparse
    template <typename FPP>
    void multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle);

    //////////////////////////////////////////////////
#endif // __COMPILE_SPMAT__




}



//----------------------------------------------//
//------------- INLINE DEFINITIONS -------------//
//----------------------------------------------//

// because of header include issues, the definition of inline functions have been moved to the end LinAlgebra_cu.hpp
//

#include "faust_linear_algebra_gpu.hpp"

#endif
