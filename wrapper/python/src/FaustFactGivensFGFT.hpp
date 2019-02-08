#include "faust_MatDense.h"
#include "faust_GivensFGFTParallel.h"
#include "faust_GivensFGFT.h"

using namespace Faust;


template<typename FPP, typename FPP2>
FaustCoreCpp<FPP>* fact_givens_fgft(const FPP* Lap, unsigned int num_rows, unsigned int num_cols, unsigned int J, unsigned int t /* end of input parameters*/, FPP* D)
{
    //TODO: optimization possible here by avoiding Lap copy in MatDense (by
    //just using the data in Lap as underlying pointer of MatDense)
    Faust::MatDense<FPP,Cpu> mat_Lap(Lap, num_rows, num_cols);

    GivensFGFT<FPP, Cpu, FPP2>* algo;


    if(t <= 0)
    {
        algo = new GivensFGFT<FPP, Cpu, FPP2>(mat_Lap, (int)J);
    }
    else
    {
        algo = new GivensFGFTParallel<FPP, Cpu, FPP2>(mat_Lap, (int)J, (int) t);
    }

    algo->compute_facts();


    Faust::Transform<FPP,Cpu> trans = std::move(algo->get_transform(true));
    TransformHelper<FPP,Cpu> *th = new TransformHelper<FPP,Cpu>(trans, true); // 1st true is for ordering (according to ascendant eigenvalues order), 2nd true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)

    //copy ordered diagonal in buffer D (allocated from the outside)
    algo->get_D(D, true);

    delete algo;

    return new FaustCoreCpp<FPP>(th);
}


