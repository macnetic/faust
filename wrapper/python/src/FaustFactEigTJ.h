#ifndef FAUST_FACT_EIGTJ_H
#define FAUST_FACT_EIGTJ_H
#include "faust_EigTJ.h"

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP* D, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false, const int err_period=100);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj_sparse(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false, const int err_period=100);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj_generic(EigTJ<FPP, Cpu, FPP2>* algo, FPP* D, const int order);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj_cplx(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP2* D, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false, const int err_period=100);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj_sparse_cplx(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t, FPP2* D, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false, const int err_period=100);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_eigtj_generic_cplx(EigTJComplex<FPP, Cpu, FPP2>* algo, FPP2* D, const int order);

template<typename FPP, typename FPP2 = float>
void svdtj( FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S /** end of output parameters */, const FPP* M, unsigned int num_rows, unsigned int num_cols, unsigned int J1, unsigned int J2, unsigned int t1, unsigned int t2, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false /* end of input parameters*/, const int err_period=100);

template<typename FPP, typename FPP2 = float>
void svdtj_sparse(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J1, unsigned int J2, unsigned int t1, unsigned int t2, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false, const int err_period=100);

template<typename FPP, typename FPP2 = float>
void svdtj_cplx( FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S /** end of output parameters */, const FPP* M, unsigned int num_rows, unsigned int num_cols, unsigned int J1, unsigned int J2, unsigned int t1, unsigned int t2, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false /* end of input parameters*/, const int err_period=100);

template<typename FPP, typename FPP2 = float>
void svdtj_sparse_cplx(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J1, unsigned int J2, unsigned int t1, unsigned int t2, unsigned int verbosity = 0, const FPP2 stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false, const int err_period=100);

template<typename FPP, typename FPP2 = float>
void create_svdtj_output(TransformHelper<FPP,Cpu> *U_, TransformHelper<FPP,Cpu>  *V_, FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, Faust::Vect<FPP,Cpu> * S_);

#include "FaustFactEigTJ.hpp"

#endif
