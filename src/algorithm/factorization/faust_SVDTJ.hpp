#include "faust_GivensFGFT.h"
#include "faust_GivensFGFTParallel.h"

using namespace Faust;

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj(MatDense<FPP, DEVICE> & dM, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
	BlasHandle<DEVICE> blas_handle;
	SpBlasHandle<DEVICE> spblas_handle;


	gemm(dM, dM, dM_M, 1.0, 0.0, 'T', 'N');
	gemm(dM, dM, dMM_, 1.0, 0.0, 'N', 'T');
	M = &dM;

	svdtj_core(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, U, V, S_);
}

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj(MatSparse<FPP, DEVICE> & sM, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP,DEVICE> dM;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
//	MatDense<FPP, DEVICE> sM_M, sMM_; // M'*M, M*M'
	BlasHandle<DEVICE> blas_handle;
	SpBlasHandle<DEVICE> spblas_handle;

	GivensFGFT<FPP, DEVICE, FPP2>* algoW1;
	GivensFGFT<FPP, DEVICE, FPP2>* algoW2;

	dM = sM;
	spgemm(sM, dM, dM_M, 1.0, 0.0, 'T', 'N');
	spgemm(sM, dM, dMM_, 1.0, 0.0, 'N', 'T');
	M = &sM;
	svdtj_core(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, U, V, S_);
}

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj_core(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, double tol, unsigned int verbosity, bool relErr, int order /* ignored anyway, kept here just in case */, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	GivensFGFT<FPP, DEVICE, FPP2>* algoW1;
	GivensFGFT<FPP, DEVICE, FPP2>* algoW2;

	if(t <= 1)
	{
		algoW1 = new GivensFGFT<FPP,DEVICE,FPP2>(dMM_, J, verbosity, tol, relErr);
		algoW2 = new GivensFGFT<FPP,DEVICE,FPP2>(dM_M, J, verbosity, tol, relErr);
	}
	else
	{
		algoW1 = new GivensFGFTParallel<FPP,DEVICE,FPP2>(dMM_, J, t, verbosity, tol, relErr);
		algoW2 = new GivensFGFTParallel<FPP,DEVICE,FPP2>(dM_M, J, t, verbosity, tol, relErr);
	}

	//TODO: parallelize with OpenMP
	algoW1->compute_facts();
	algoW2->compute_facts();

	Faust::Vect<FPP,DEVICE> S(M->getNbRow());


	Faust::Transform<FPP,DEVICE> transW1 = std::move(algoW1->get_transform(/*order*/ 0)); // 0 for no order
	TransformHelper<FPP,DEVICE> *thW1 = new TransformHelper<FPP,DEVICE>(transW1, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


	Faust::Transform<FPP,DEVICE> transW2 = std::move(algoW2->get_transform(/* order*/ 0)); //0 for no order
	TransformHelper<FPP,DEVICE> *thW2 = new TransformHelper<FPP,DEVICE>(transW2, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


	// compute S = W1'*M*W2 = W1'*(W2^T*M)^T
	dM.transpose();
	Faust::MatDense<FPP,DEVICE> MW2 = thW2->multiply(dM, /* transpose */ true);
	MW2.transpose();
	Faust::MatDense<FPP,DEVICE> W1_MW2 = thW1->multiply(MW2, /* transpose */ true);

	// create diagonal vector
	for(int i=0;i<S.size();i++){
		S.getData()[i] = W1_MW2(i,i);
	}

	// order D descendantly according to the abs value
	// and change the sign when the value is negative
	// it gives a signed permutation matrix P to append to W1, abs(P2) is append to W2
	vector<int> ord_indices;
	Faust::Vect<FPP,DEVICE>* ordered_S = new Faust::Vect<FPP,DEVICE>(S.size());
	vector<FPP> values(S.size());
	vector<FPP> values2(S.size());
	vector<int> col_ids(S.size());
	ord_indices.resize(0);
	order = 1;
	for(int i=0;i<S.size();i++)
		ord_indices.push_back(i);
	sort(ord_indices.begin(), ord_indices.end(), [S, &order](int i, int j) {
			return Faust::fabs(S.getData()[i]) > Faust::fabs(S.getData()[j]);
			});
	for(int i=0;i<ord_indices.size();i++)
	{
		col_ids[i] = i;
		ordered_S->getData()[i] = Faust::fabs(S.getData()[ord_indices[i]]);
		if(S.getData()[ord_indices[i]] < 0)
			values[i] = -1;
		else
			values[i] = 1;
		values2[i] = 1;
	}
	Faust::MatSparse<FPP, DEVICE>* PS = new Faust::MatSparse<FPP, DEVICE>(ord_indices, col_ids, values, M->getNbRow(), M->getNbCol());
	MatGeneric<FPP,Cpu>* lf = (MatGeneric<FPP,Cpu>*)(thW1->get_fact_addr(thW1->size()-1));
	lf->multiplyRight(*PS);
	// thW1->push_back(PS);

	Faust::MatSparse<FPP, DEVICE>* P = new Faust::MatSparse<FPP, DEVICE>(ord_indices, col_ids, values2, M->getNbRow(), M->getNbCol());
	lf = (MatGeneric<FPP,Cpu>*)(thW2->get_fact_addr(thW2->size()-1));
	lf->multiplyRight(*P);
	// thW2->push_back(P);
	delete algoW1;
	delete algoW2;

	*U = thW1;
	*V = thW2;
	*S_ = ordered_S;
}

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj_cplx(MatDense<FPP, DEVICE> & dM, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
	BlasHandle<DEVICE> blas_handle;
	SpBlasHandle<DEVICE> spblas_handle;


	gemm(dM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	gemm(dM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &dM;

	svdtj_core_cplx(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, U, V, S_);
}

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj_cplx(MatSparse<FPP, DEVICE> & sM, int J, int t, double tol, unsigned int verbosity, bool relErr, int order, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	MatGeneric<FPP,DEVICE>* M;
	MatDense<FPP,DEVICE> dM;
	MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
//	MatDense<FPP, DEVICE> sM_M, sMM_; // M'*M, M*M'
	BlasHandle<DEVICE> blas_handle;
	SpBlasHandle<DEVICE> spblas_handle;

	GivensFGFT<FPP, DEVICE, FPP2>* algoW1;
	GivensFGFT<FPP, DEVICE, FPP2>* algoW2;

	dM = sM;
	spgemm(sM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
	spgemm(sM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
	M = &sM;
	svdtj_core_cplx(M, dM, dM_M, dMM_, J, t, tol, verbosity, relErr, order, U, V, S_);
}

	template<typename FPP, Device DEVICE, typename FPP2>
void Faust::svdtj_core_cplx(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, double tol, unsigned int verbosity, bool relErr, int order /* ignored anyway, kept here just in case */, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Faust::Vect<FPP,DEVICE> ** S_)
{
	GivensFGFTComplex<FPP, DEVICE, FPP2>* algoW1;
	GivensFGFTComplex<FPP, DEVICE, FPP2>* algoW2;

	if(t <= 1)
	{
		algoW1 = new GivensFGFTComplex<FPP,DEVICE,FPP2>(dMM_, J, verbosity, tol, relErr);
		algoW2 = new GivensFGFTComplex<FPP,DEVICE,FPP2>(dM_M, J, verbosity, tol, relErr);
	}
	else
	{
		algoW1 = new GivensFGFTParallelComplex<FPP,DEVICE,FPP2>(dMM_, J, t, verbosity, tol, relErr);
		algoW2 = new GivensFGFTParallelComplex<FPP,DEVICE,FPP2>(dM_M, J, t, verbosity, tol, relErr);
	}

	//TODO: parallelize with OpenMP
	algoW1->compute_facts();
	algoW2->compute_facts();

	Faust::Vect<FPP,DEVICE> S(M->getNbRow());


	Faust::Transform<FPP,DEVICE> transW1 = std::move(algoW1->get_transform(/*order*/ 0)); // 0 for no order
	TransformHelper<FPP,DEVICE> *thW1 = new TransformHelper<FPP,DEVICE>(transW1, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


	Faust::Transform<FPP,DEVICE> transW2 = std::move(algoW2->get_transform(/* order*/ 0)); //0 for no order
	TransformHelper<FPP,DEVICE> *thW2 = new TransformHelper<FPP,DEVICE>(transW2, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)


	// compute S = W1'*M*W2 = W1'*(W2^T*M)^T
	dM.transpose();
	Faust::MatDense<FPP,DEVICE> MW2 = thW2->multiply(dM, /* transpose */ true);
	MW2.transpose();
	Faust::MatDense<FPP,DEVICE> W1_MW2 = thW1->multiply(MW2, /* transpose */ true);

	// create diagonal vector
	for(int i=0;i<S.size();i++){
		S.getData()[i] = W1_MW2(i,i);
	}

	// order D descendantly according to the abs value
	// and change the sign when the value is negative
	// it gives a signed permutation matrix P to append to W1, abs(P2) is append to W2
	vector<int> ord_indices;
	Faust::Vect<FPP,DEVICE>* ordered_S = new Faust::Vect<FPP,DEVICE>(S.size());
	vector<FPP> values(S.size());
	vector<FPP> values2(S.size());
	vector<int> col_ids(S.size());
	ord_indices.resize(0);
	order = 1;
	for(int i=0;i<S.size();i++)
		ord_indices.push_back(i);
	sort(ord_indices.begin(), ord_indices.end(), [S, &order](int i, int j) {
			return Faust::fabs(S.getData()[i]) > Faust::fabs(S.getData()[j]);
			});
	for(int i=0;i<ord_indices.size();i++)
	{
		col_ids[i] = i;
		ordered_S->getData()[i] = Faust::fabs(S.getData()[ord_indices[i]]);
		if(complex<float>(S.getData()[ord_indices[i]]).real() < 0)
			values[i] = -1;
		else
			values[i] = 1;
		values2[i] = 1;
	}
	Faust::MatSparse<FPP, DEVICE>* PS = new Faust::MatSparse<FPP, DEVICE>(ord_indices, col_ids, values, M->getNbRow(), M->getNbCol());
	MatGeneric<FPP,Cpu>* lf = (MatGeneric<FPP,Cpu>*)(thW1->get_fact_addr(thW1->size()-1));
	lf->multiplyRight(*PS);
	// thW1->push_back(PS);

	Faust::MatSparse<FPP, DEVICE>* P = new Faust::MatSparse<FPP, DEVICE>(ord_indices, col_ids, values2, M->getNbRow(), M->getNbCol());
	lf = (MatGeneric<FPP,Cpu>*)(thW2->get_fact_addr(thW2->size()-1));
	lf->multiplyRight(*P);
	// thW2->push_back(P);
	delete algoW1;
	delete algoW2;

	*U = thW1;
	*V = thW2;
	*S_ = ordered_S;
}
