#ifndef FAUST_FACT_GIVENS_FGFT_H
#define FAUST_FACT_GIVENS_FGFT_H
#include "faust_GivensFGFT.h"

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_generic(GivensFGFT<FPP, Cpu, FPP2>* algo, FPP* D, const int order, const bool enable_large_Faust=false);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_cplx(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP2* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse_cplx(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t, FPP2* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */, const bool enable_large_Faust=false);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_generic_cplx(GivensFGFTComplex<FPP, Cpu, FPP2>* algo, FPP2* D, const int order, const bool enable_large_Faust=false);

template<typename FPP, typename FPP2 = float>
void svdtj( FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S /** end of output parameters */, const FPP* M, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false /* end of input parameters*/);

template<typename FPP, typename FPP2 = float>
void svdtj_sparse(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J, unsigned int t, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false);

template<typename FPP, typename FPP2 = float>
void svdtj_cplx( FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S /** end of output parameters */, const FPP* M, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false /* end of input parameters*/);

template<typename FPP, typename FPP2 = float>
void svdtj_sparse_cplx(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J, unsigned int t, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const bool enable_large_Faust = false);

template<typename FPP, typename FPP2 = float>
void create_svdtj_output(TransformHelper<FPP,Cpu> *U_, TransformHelper<FPP,Cpu>  *V_, FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, Faust::Vect<FPP,Cpu> * S_);

#include "FaustFactGivensFGFT.hpp"

#endif
