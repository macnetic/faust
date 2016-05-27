#ifndef __FAUST_LINEAR_OPERATOR_H__
#define __FAUST_LINEAR_OPERATOR_H__

#include "faust_constant.h"


// modif AL AL
//template<typename T,Device DEVICE> class MatDense;
//template<typename FPP,Device DEVICE> class MatDense;


/**
 * \class Faust::LinearOperator faust_LinearOperator.h
 * \brief This virtual template class represent all the different type of matrix that are used<br> (Dense matrix (Faust::MatDense), Sparse matrix (Faust::MatSparse),FAUST (Faust::Transform)). <br>
 *
 * \tparam T scalar numeric type, e.g float or double
*/


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{
    // modif AL AL
    template<typename FPP,Device DEVICE>
    class MatDense;

    template<typename FPP,Device DEVICE>
    class LinearOperator
    {
        public:

        virtual faust_unsigned_int getNbRow() const=0;
        virtual faust_unsigned_int getNbCol() const=0;


        /*!
        *  \brief Transpose the Faust::MatDense
        */
        virtual void transpose()=0;

        /*! \brief faust_gemm — Performs a matrix calculation of the form alpha . (*this).B +beta.C, with B and C 2 faust_matrix and two real constants you supply. <br>

        * C = alpha .op(this).op(B) + beta.C; <br>
        * op(A) = A if typeA='N', op(A) = transpose(A) if typeA='T'<br>
        * op(B) = B if typeB='N', op(B) = transpose(B) if typeB='T'<br>
        * The Object C must be different of B. <br>
        * \param B : Dense Matrix
        * \param C : Dense Matrix
        * \param alpha : Template scalar coefficient
        * \param beta : Template scalar coefficient
        */
        virtual void faust_gemm(const Faust::MatDense<FPP,DEVICE> & B, Faust::MatDense<FPP,DEVICE> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const=0;

    };
}

#endif
