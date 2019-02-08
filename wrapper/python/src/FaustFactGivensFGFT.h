#ifndef FAUST_FACT_GIVENS_FGFT_H
#define FAUST_FACT_GIVENS_FGFT_H

template<typename FPP, typename FPP2 = float>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D);

#include "FaustFactGivensFGFT.hpp"

#endif
