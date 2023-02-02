#include "faust_GivensFGFT.h"
#include "faust_GivensFGFTParallel.h"
#include "faust_GivensFGFTParallelComplex.h"
#include "faust_GivensFGFTComplex.h"
#include "faust_constant.h"
#include "faust_GivensFGFTGen.h"
#include <type_traits>
using namespace Faust;

	template<typename FPP, FDevice DEVICE, typename FPP2>
void Faust::svdtj(MatDense<FPP, DEVICE> & dM, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'


	gemm(dM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	gemm(dM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &dM;

	if(verbosity)
	{
		std::cout << "Faust::svdtj input conf:" << std::endl;
		std::cout << " J: " << J << std::endl;
		std::cout << " t: " << t << std::endl;
		std::cout << " tol: " << tol << std::endl;
		std::cout << " relErr: " << relErr << std::endl;
		std::cout << " order: " << order << std::endl;
		std::cout << " enable_large_Faust: " << enable_large_Faust << std::endl;
		std::cout << " matrix norm: " << dM.norm() << std::endl;
	}

	svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
}

	template<typename FPP, FDevice DEVICE, typename FPP2>
void Faust::svdtj(MatSparse<FPP, DEVICE> & sM, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP,DEVICE> dM;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
//	MatDense<FPP, DEVICE> sM_M, sMM_; // M'*M, M*M'

	GivensFGFT<FPP, DEVICE, FPP2>* algoW1;
	GivensFGFT<FPP, DEVICE, FPP2>* algoW2;

	dM = sM;
	spgemm(sM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	spgemm(sM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &sM;
	svdtj_core_gen<FPP,DEVICE, FPP2>(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
}

	template<typename FPP, FDevice DEVICE, typename FPP2>
void Faust::svdtj_cplx(MatDense<FPP, DEVICE> & dM, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'


	gemm(dM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	gemm(dM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &dM;

//	svdtj_core_cplx(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
	svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);

}

	template<typename FPP, FDevice DEVICE, typename FPP2>
void Faust::svdtj_cplx(MatSparse<FPP, DEVICE> & sM, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP,DEVICE> dM;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
//	MatDense<FPP, DEVICE> sM_M, sMM_; // M'*M, M*M'


	dM = sM;
	spgemm(sM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	spgemm(sM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &sM;
//	svdtj_core_cplx(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
	svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
}

template<typename FPP, FDevice DEVICE, typename FPP2>
void Faust::svdtj_core_gen(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order /* ignored anyway, kept here just in case */, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW1;
	GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW2;

	instantiate_algos(&algoW1, &algoW2, dM_M, dMM_, J, t, verbosity, tol, relErr, enable_large_Faust);

	try {

		//TODO: parallelize with OpenMP ?
		algoW1->compute_facts();
		algoW2->compute_facts();

		auto m = M->getNbRow();
		auto n = M->getNbCol();
		auto min_mn = m > n?n:m;

		Vect<FPP,DEVICE> S(min_mn); // zero-initialized


		Transform<FPP,DEVICE> transW1 = std::move(algoW1->get_transform(/*order*/ -1)); // 0 for no order, -1 descending order


		Transform<FPP,DEVICE> transW2 = std::move(algoW2->get_transform(/* order*/ -1));


		// compute S = W1'*M*W2 = W1'*(W2'*M')'
		dM.adjoint();
		MatDense<FPP, DEVICE> MW2 = transW2.multiply(dM, 'H');
		MW2.adjoint();
		Faust::MatDense<FPP,DEVICE> W1_MW2 = transW1.multiply(MW2, 'H');


		TransformHelper<FPP,DEVICE> *thW1 = new TransformHelper<FPP,DEVICE>(transW1, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)

		TransformHelper<FPP,DEVICE> *thW2 = new TransformHelper<FPP,DEVICE>(transW2, true); // the same

#if DEBUG_SVDTJ
		thW1->save_mat_file("/tmp/W1_cpp.mat");
		thW2->save_mat_file("/tmp/W2_cpp.mat");

		MW2.save_to_mat_file("/tmp/MW2_cpp.mat", "MW2_cpp");

		W1_MW2.save_to_mat_file("/tmp/W1_MW2_cpp.mat", "W1_MW2_cpp");
#endif


		// set diagonal vector
		for(int i=0;i<S.size();i++) //TODO: OpenMP?
			S.getData()[i] = W1_MW2(i,i);

#if DEBUG_SVDTJ
		MatDense<FPP, Cpu> S1_mat(S.size(), 1, S.getData());
		S1_mat.save_to_mat_file("/tmp/S1_cpp.mat", "S1_cpp");
#endif

		// order D descendingly according to the abs values
		// and change the sign when the value is negative
		// it gives a signed permutation matrix PS to multiply to W1, and a permutation P is multiplied to W2
		vector<int> ord_indices, indices;
		Vect<FPP,DEVICE>* ordered_S = new Vect<FPP,DEVICE>(S.size());
		vector<FPP> values(m, 0);
		vector<FPP> values2(n, 0);
		vector<int> col_ids(m);
		vector<int> col_ids2(n);
		std::iota(col_ids.begin(), col_ids.end(), 0);
		std::iota(col_ids2.begin(), col_ids2.end(), 0);
		ord_indices.resize(S.size());
		std::iota(ord_indices.begin(), ord_indices.end(), 0);
		sort(ord_indices.begin(), ord_indices.end(), [S, &order](int i, int j) {
				return Faust::fabs(S.getData()[i]) > Faust::fabs(S.getData()[j]);
				});
		vector<int> row_ids(m);
		vector<int> row_ids2(n);
		std::copy(ord_indices.begin(), ord_indices.end(), row_ids.begin());
		std::copy(ord_indices.begin(), ord_indices.end(), row_ids2.begin());
		std::iota(row_ids.begin()+min_mn, row_ids.end(), min_mn); // min_mn == S.size()
		std::iota(row_ids2.begin()+min_mn, row_ids2.end(), min_mn);
		for(int i=0;i<ord_indices.size();i++)
		{
			ordered_S->getData()[i] = Faust::fabs(S.getData()[ord_indices[i]]);
			if(complex<float>(S.getData()[ord_indices[i]]).real() < 0)
			{
				values[i] = -1;
			}
			else
				values[i] = 1;
			values2[i] = 1;
		}
		MatSparse<FPP, DEVICE>* PS = new MatSparse<FPP, DEVICE>(row_ids, col_ids, values, m, m);
		MatGeneric<FPP,Cpu>* lf = (MatGeneric<FPP,Cpu>*)(thW1->get_fact_addr(thW1->size()-1));
		lf->multiplyRight(*PS);
//		thW1->push_back(PS); // multiplying is more efficient than pushing a new factor

		MatSparse<FPP, DEVICE>* P = new MatSparse<FPP, DEVICE>(row_ids2, col_ids2, values2, n, n);
		lf = (MatGeneric<FPP,Cpu>*)(thW2->get_fact_addr(thW2->size()-1));
		lf->multiplyRight(*P);
		// thW2->push_back(P);
		delete PS;
		delete P;
		*U = thW1;
		*V = thW2;
		*S_ = ordered_S;

#if DEBUG_SVDTJ
		MatDense<FPP, Cpu> So_mat(ordered_S->size(), 1, ordered_S->getData());
		So_mat.save_to_mat_file("/tmp/So_cpp.mat", "So_cpp");
#endif
	}
	catch(out_of_range e)
	{
		// algorithm stopped before the first iteration (nGivens too large)
	}

	delete algoW1;
	delete algoW2;

}
