#ifndef FAUST_FACT_GIVENS_FGFT_H
#define FAUST_FACT_GIVENS_FGFT_H
#include <faust_GivensFGFT.h>

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true);

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse(const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity = 0, const double stoppingError = 0.0, const bool errIsRel = true);


template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft_generic(GivensFGFT<FPP, Cpu, FPP2>* algo, FPP* D);


#include "FaustFactGivensFGFT.hpp"

#endif
