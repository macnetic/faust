#include "faust_MatDense.h"
#include "faust_GivensFGFTParallel.h"
#include "faust_GivensFGFT.h"
#include "faust_GivensFGFTComplex.h"
using namespace Faust;

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const double stoppingError, const bool errIsRel, const int order)
{
      Faust::MatSparse<FPP, Cpu> mat_Lap(nnz, nrows, ncols, data, id_col, row_ptr);
      GivensFGFT<FPP, Cpu, FPP2>* algo;
      if(t <= 1)
      {
          algo = new GivensFGFT<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel);
      }
      else
      {
          algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel);
      }
      return fact_givens_fgft_generic(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const double stoppingError, const bool errIsRel, const int order)
{
    //TODO: optimization possible here by avoiding Lap copy in MatDense (by
    //just using the data in Lap as underlying pointer of MatDense)
    Faust::MatDense<FPP,Cpu> mat_Lap(Lap, num_rows, num_cols);

    GivensFGFT<FPP, Cpu, FPP2>* algo;
    if(t <= 1)
    {
        algo = new GivensFGFT<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel);
    }
    else
    {
        algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel);
    }
    return fact_givens_fgft_generic(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_generic(GivensFGFT<FPP, Cpu, FPP2>* algo, FPP* D, const int order)
{


    algo->compute_facts();

    Faust::Transform<FPP,Cpu> trans = std::move(algo->get_transform(order));
    TransformHelper<FPP,Cpu> *th = new TransformHelper<FPP,Cpu>(trans, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


    //copy ordered diagonal in buffer D (allocated from the outside)
    algo->get_D(D, order);

    delete algo;

    return new FaustCoreCpp<FPP>(th);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse_cplx(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const double stoppingError, const bool errIsRel, const int order)
{
      Faust::MatSparse<FPP, Cpu> mat_Lap(nnz, nrows, ncols, data, id_col, row_ptr);
      GivensFGFTComplex<FPP, Cpu, FPP2>* algo;
      if(t <= 1)
      {
          algo = new GivensFGFTComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel);
      }
      else
      {
//          algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel);
          return nullptr;
      }
      return fact_givens_fgft_generic_cplx(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_cplx(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const double stoppingError, const bool errIsRel, const int order)
{
    //TODO: optimization possible here by avoiding Lap copy in MatDense (by
    //just using the data in Lap as underlying pointer of MatDense)
    Faust::MatDense<FPP,Cpu> mat_Lap(Lap, num_rows, num_cols);

    GivensFGFTComplex<FPP, Cpu, FPP2>* algo;
    if(t <= 1)
    {
        algo = new GivensFGFTComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel);
    }
    else
    {
//        algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel);
        return nullptr;
    }
    return fact_givens_fgft_generic_cplx(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_generic_cplx(GivensFGFTComplex<FPP, Cpu, FPP2>* algo, FPP* D, const int order)
{


    algo->compute_facts();

    Faust::Transform<FPP,Cpu> trans = std::move(algo->get_transform(order));
    TransformHelper<FPP,Cpu> *th = new TransformHelper<FPP,Cpu>(trans, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


    //copy ordered diagonal in buffer D (allocated from the outside)
    algo->get_D(D, order);

    delete algo;

    return new FaustCoreCpp<FPP>(th);
}
