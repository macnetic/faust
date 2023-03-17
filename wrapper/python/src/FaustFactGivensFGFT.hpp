#include "faust_MatDense.h"
#include "faust_GivensFGFTParallel.h"
#include "faust_GivensFGFT.h"
#include "faust_GivensFGFTComplex.h"
#include "faust_GivensFGFTParallelComplex.h"
#include "faust_SVDTJ.h"
using namespace Faust;

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const int order, const bool enable_large_Faust)
{
      Faust::MatSparse<FPP, Cpu> mat_Lap(nnz, nrows, ncols, data, id_col, row_ptr);
      GivensFGFT<FPP, Cpu, FPP2>* algo;
      if(t <= 1)
      {
          algo = new GivensFGFT<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel, enable_large_Faust);
      }
      else
      {
          algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel, enable_large_Faust);
      }
      return fact_givens_fgft_generic(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const int order, const bool enable_large_Faust)
{
    //TODO: optimization possible here by avoiding Lap copy in MatDense (by
    //just using the data in Lap as underlying pointer of MatDense)
    Faust::MatDense<FPP,Cpu> mat_Lap(Lap, num_rows, num_cols);

    GivensFGFT<FPP, Cpu, FPP2>* algo;
    if(t <= 1)
    {
        algo = new GivensFGFT<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel, enable_large_Faust);
    }
    else
    {
        algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel, enable_large_Faust);
    }
    return fact_givens_fgft_generic(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_generic(GivensFGFT<FPP, Cpu, FPP2>* algo, FPP* D, const int order, const bool enable_large_Faust)
{

    FaustCoreCpp<FPP>* fc = nullptr;

    algo->compute_facts();

    try {
        Faust::Transform<FPP,Cpu> trans = std::move(algo->get_transform(order));
        TransformHelper<FPP,Cpu> *th = new TransformHelper<FPP,Cpu>(trans, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


        //copy ordered diagonal in buffer D (allocated from the outside)
        algo->get_D(D, order);

        fc = new FaustCoreCpp<FPP>(th);
    }
    catch(out_of_range e)
    {

    }
    delete algo;
    return fc;
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_sparse_cplx(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols,
        unsigned int J, unsigned int t /* end of input parameters*/, FPP2* D, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const int order, const bool enable_large_Faust)
{
      Faust::MatSparse<FPP, Cpu> mat_Lap(nnz, nrows, ncols, data, id_col, row_ptr);
      GivensFGFTComplex<FPP, Cpu, FPP2>* algo;
      if(t <= 1)
      {
          algo = new GivensFGFTComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel, enable_large_Faust);
      }
      else
      {
          algo = new GivensFGFTParallelComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel, enable_large_Faust);
      }
      return fact_givens_fgft_generic_cplx(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_cplx(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP2* D, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const int order, const bool enable_large_Faust)
{
    //TODO: optimization possible here by avoiding Lap copy in MatDense (by
    //just using the data in Lap as underlying pointer of MatDense)
    Faust::MatDense<FPP,Cpu> mat_Lap(Lap, num_rows, num_cols);

    GivensFGFTComplex<FPP, Cpu, FPP2>* algo;
    if(t <= 1)
    {
        algo = new GivensFGFTComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, verbosity, stoppingError, errIsRel, enable_large_Faust);
    }
    else
    {
        algo = new GivensFGFTParallelComplex<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t, verbosity, stoppingError, errIsRel, enable_large_Faust);
    }
    return fact_givens_fgft_generic_cplx(algo, D, order);
}

template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft_generic_cplx(GivensFGFTComplex<FPP, Cpu, FPP2>* algo, FPP2* D, const int order, const bool enable_large_Faust)
{

    FaustCoreCpp<FPP>* fc = nullptr;
    algo->compute_facts();

    try
    {
        Faust::Transform<FPP,Cpu> trans = std::move(algo->get_transform(order));
        TransformHelper<FPP,Cpu> *th = new TransformHelper<FPP,Cpu>(trans, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


        //copy ordered diagonal in buffer D (allocated from the outside)
        algo->get_D(D, order);
        fc =new FaustCoreCpp<FPP>(th);
    }
    catch(out_of_range e)
    {
    }
    delete algo;

    return fc;
}

template<typename FPP, typename FPP2>
void svdtj(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* M_data, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const bool enable_large_Faust)
{
    Faust::MatDense<FPP,Cpu> M(M_data, (faust_unsigned_int) num_rows, (faust_unsigned_int) num_cols);
    TransformHelper<FPP,Cpu> *U_ = nullptr,  *V_ = nullptr;
    Faust::Vect<FPP,Cpu> * S_ = nullptr;
    svdtj(M, J, t, stoppingError, verbosity, errIsRel, -1 /* descending order */, enable_large_Faust ,&U_, &V_, &S_);
    create_svdtj_output(U_, V_, U, V, S, S_);
}

template<typename FPP, typename FPP2>
void svdtj_sparse(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J, unsigned int t, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const bool enable_large_Faust)
{
    Faust::MatSparse<FPP, Cpu> M(nnz, nrows, ncols, data, id_col, row_ptr);
    TransformHelper<FPP,Cpu> *U_ = nullptr,  *V_ = nullptr;
    Faust::Vect<FPP,Cpu> * S_ = nullptr;
    svdtj(M, J, t, stoppingError, verbosity, errIsRel, -1 /* descending order */, enable_large_Faust, &U_, &V_, &S_);
    create_svdtj_output(U_, V_, U, V, S, S_);
}

template<typename FPP, typename FPP2>
void svdtj_cplx(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* M_data, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const bool enable_large_Faust)
{
    Faust::MatDense<FPP,Cpu> M(M_data, (faust_unsigned_int) num_rows, (faust_unsigned_int) num_cols);
    TransformHelper<FPP,Cpu> *U_ = nullptr,  *V_ = nullptr;
    Faust::Vect<FPP,Cpu> * S_ = nullptr;
    svdtj_cplx(M, J, t, stoppingError, verbosity, errIsRel, -1 /* descending order */, enable_large_Faust, &U_, &V_, &S_);
    create_svdtj_output(U_, V_, U, V, S, S_);
}

template<typename FPP, typename FPP2>
void svdtj_sparse_cplx(FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, /*start of input parameters*/ const FPP* data, int* row_ptr, int* id_col, int nnz, int nrows, int ncols, unsigned int J, unsigned int t, unsigned int verbosity, const FPP2 stoppingError, const bool errIsRel, const bool enable_large_Faust)
{
    Faust::MatSparse<FPP, Cpu> M(nnz, nrows, ncols, data, id_col, row_ptr);
    TransformHelper<FPP,Cpu> *U_ = nullptr,  *V_ = nullptr;
    Faust::Vect<FPP,Cpu> * S_ = nullptr;
    svdtj_cplx(M, J, t, stoppingError, verbosity, errIsRel, -1 /* descending order */, enable_large_Faust, &U_, &V_, &S_);
    create_svdtj_output(U_, V_, U, V, S, S_);
}

template<typename FPP, typename FPP2>
void create_svdtj_output(TransformHelper<FPP,Cpu> *U_, TransformHelper<FPP,Cpu>  *V_, FaustCoreCpp<FPP>** U, FaustCoreCpp<FPP> **V, FPP* S, Faust::Vect<FPP,Cpu> * S_)
{
    if(U_ != nullptr && V_ != nullptr)
    {
        *U = new FaustCoreCpp<FPP>(U_);
        *V = new FaustCoreCpp<FPP>(V_);
        memcpy(S, S_->getData(), sizeof(FPP)* S_->size());
    }
    if(S_ != nullptr)
        delete S_;
}
