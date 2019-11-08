#ifndef FAUST_FACT_GIVENS_FGFT_H
#define FAUST_FACT_GIVENS_FGFT_H
#include "faust_GivensFGFT.h"

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_generic(GivensFGFT<FPP, Cpu, FPP2>* algo, FPP* D, const int order);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_cplx(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse_cplx(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true, const int order = 1 /* ascendant order for eigenvalues */);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_generic_cplx(GivensFGFTComplex<FPP, Cpu, FPP2>* algo, FPP* D, const int order);

template<typename FPP, typename FPP2 = float>
void svdtj( FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S /** end of output parameters */, const FPP* M, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true /* end of input parameters*/);

#include "FaustFactGivensFGFT.hpp"

#endif
