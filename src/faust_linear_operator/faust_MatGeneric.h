#ifndef __FAUST_MAT_GENERIC_H__
#define __FAUST_MAT_GENERIC_H__

#include "faust_constant.h"
#include "faust_LinearOperator.h"
/**
 * \class MatGeneric faust_MatGeneric.h
 * \brief This MatGeneric class serves as a base class for the derived class Faust::MatDense and Faust::MatSparse .
 * 	Member variable dim1 and dim2 correspond to the dimension of the matrix
 *
*/


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    // modif AL AL
    template<typename FPP,Device DEVICE>
    class LinearOperator;

    template<typename FPP,Device DEVICE>
    class MatGeneric : public Faust::LinearOperator<FPP,DEVICE>
    {
        public:
        MatGeneric() : dim1(0), dim2(0) {}
        MatGeneric(faust_unsigned_int dim1_, faust_unsigned_int dim2_) : dim1(dim1_), dim2(dim2_){}

        faust_unsigned_int getNbRow() const {return dim1;}
        faust_unsigned_int getNbCol() const {return dim2;}

        void resize(const faust_unsigned_int dim1_,const faust_unsigned_int dim2_){dim1=dim1_;dim2=dim2_;}

        protected:
        faust_unsigned_int dim1;
        faust_unsigned_int dim2;
    };

}
#endif
