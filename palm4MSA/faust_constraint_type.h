#ifndef __FAUST_CONSTRAINT_TYPE__
#define __FAUST_CONSTRAINT_TYPE__

class faust_constraint_int;
template<typename T> class faust_constraint_real;
template<typename T> class faust_constraint_mat;




/*!
\brief This structure defined the type of constraint. See following table for more precision about the type of constraint. <br>
 <img src="../../doc/html/constraint.png" alt="constraint parameters" width=800px />
*/
template<typename T>
struct constraint_type
{
	typedef const  faust_constraint_real<T>  constraint_type_normcol;
	typedef const  faust_constraint_real<T>  constraint_type_normlin;
	typedef const  faust_constraint_mat<T>   constraint_type_supp;
	typedef const  faust_constraint_mat<T>   constraint_type_const;
	typedef const  faust_constraint_int   constraint_type_sp;
	typedef const  faust_constraint_int   constraint_type_spcol;
	typedef const  faust_constraint_int   constraint_type_splin;
	typedef const  faust_constraint_int   constraint_type_splincol;
	typedef const  faust_constraint_int   constraint_type_l0pen;
	typedef const  faust_constraint_int   constraint_type_l1pen;
	typedef const  faust_constraint_int   constraint_type_wav;
	typedef const  faust_constraint_int   constraint_type_sp_pos;
	typedef const  faust_constraint_int   constraint_type_blkdiag;
	typedef const  faust_constraint_int   constraint_type_splin_test;

	typedef const  faust_constraint_int   constraint_type_toeplitz;
};





#endif
