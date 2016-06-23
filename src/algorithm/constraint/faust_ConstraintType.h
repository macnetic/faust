#ifndef __FAUST_CONSTRAINT_TYPE_H__
#define __FAUST_CONSTRAINT_TYPE_H__


template<typename FPP,Device DEVICE> class ConstraintInt;
template<typename FPP,Device DEVICE> class ConstraintFPP;
template<typename FPP,Device DEVICE> class ConstraintMat;


//! \struct ConstraintType
//! \brief This structure defined the type of constraint. See following table for more precision about the type of constraint. <br>
 //! <img src="../../doc/html/constraint.png" alt="constraint parameters" width=800px />
template<typename FPP,Device DEVICE>
struct ConstraintType
{

	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   ConstraintTypeNormcol;
	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   ConstraintTypeNormlin;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   ConstraintTypeSupp;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   ConstraintTypeConst;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSp;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSpcol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSplin;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSplincol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeSpPos;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   ConstraintTypeBlkdiag;

};


#endif
