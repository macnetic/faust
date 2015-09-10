#ifndef __FAUST_CONSTRAINT_TYPE__
#define __FAUST_CONSTRAINT_TYPE__

class faust_constraint_int;
class faust_constraint_real;
class faust_constraint_mat;


typedef const  faust_constraint_int   constraint_type_sp;
typedef const  faust_constraint_int   constraint_type_spcol;
typedef const  faust_constraint_int   constraint_type_splin;
typedef const  faust_constraint_real  constraint_type_normcol;
typedef const  faust_constraint_int   constraint_type_splincol;
typedef const  faust_constraint_int   constraint_type_l0pen;
typedef const  faust_constraint_int   constraint_type_l1pen;
typedef const  faust_constraint_mat   constraint_type_const;
typedef const  faust_constraint_int   constraint_type_wav;
typedef const  faust_constraint_int   constraint_type_sp_pos;
typedef const  faust_constraint_int   constraint_type_blkdiag;
typedef const  faust_constraint_int   constraint_type_splin_test;
typedef const  faust_constraint_mat   constraint_type_supp;
typedef const  faust_constraint_real  constraint_type_normlin;
typedef const  faust_constraint_int   constraint_type_toeplitz;

#endif
