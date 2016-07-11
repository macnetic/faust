#include "faust_exception.h"
template<typename FPP,Device DEVICE>
void Faust::MatGeneric<FPP,DEVICE>::setOp(const char op, faust_unsigned_int& nbRowOp, faust_unsigned_int& nbColOp)const
{
    if(op == 'N')
    {
        nbRowOp=getNbRow();
        nbColOp=getNbCol();
    }
    else if(op == 'T')
    {
        nbRowOp=getNbCol();
        nbColOp=getNbRow();
    }
    else
        handleError("Faust::MatGeneric::","setOp : invalid character");
}
