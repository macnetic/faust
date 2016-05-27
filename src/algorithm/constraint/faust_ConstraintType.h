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

	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   constraint_type_normcol;
	typedef const  Faust::ConstraintFPP<FPP,DEVICE>   constraint_type_normlin;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   constraint_type_supp;
	typedef const  Faust::ConstraintMat<FPP,DEVICE>   constraint_type_const;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_sp;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_spcol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_splin;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_splincol;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_sp_pos;
	typedef const  Faust::ConstraintInt<FPP,DEVICE>   constraint_type_blkdiag;

};


#endif
