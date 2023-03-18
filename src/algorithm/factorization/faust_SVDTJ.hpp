#include "faust_GivensFGFT.h"
#include "faust_GivensFGFTParallel.h"
#include "faust_GivensFGFTParallelComplex.h"
#include "faust_GivensFGFTComplex.h"
#include "faust_constant.h"
#include "faust_GivensFGFTGen.h"
#include <type_traits>

namespace Faust
{

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj(MatDense<FPP, DEVICE> & dM, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
		{
			MatGeneric<FPP,DEVICE>* M;
			MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'


			gemm(dM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
			gemm(dM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
			M = &dM;

			if(verbosity)
			{
				std::cout << "svdtj input conf:" << std::endl;
				std::cout << " J1: " << J1 << std::endl;
				std::cout << " J2: " << J2 << std::endl;
				std::cout << " t1: " << t1 << std::endl;
				std::cout << " t2: " << t2 << std::endl;
				std::cout << " tol: " << tol << std::endl;
				std::cout << " relErr: " << relErr << std::endl;
				std::cout << " order: " << order << std::endl;
				std::cout << " enable_large_Faust: " << enable_large_Faust << std::endl;
				std::cout << " matrix norm: " << dM.norm() << std::endl;
			}

			svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj(MatSparse<FPP, DEVICE> & sM, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
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
			svdtj_core_gen<FPP,DEVICE, FPP2>(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj_cplx(MatDense<FPP, DEVICE> & dM, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
		{
			MatGeneric<FPP,DEVICE>* M;
			MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'


			gemm(dM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
			gemm(dM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
			M = &dM;

			//	svdtj_core_cplx(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
			svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);

		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj_cplx(MatSparse<FPP, DEVICE> & sM, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
		{
			MatGeneric<FPP,DEVICE>* M;
			MatDense<FPP,DEVICE> dM;
			MatDense<FPP, DEVICE> dM_M, dMM_; // M'*M, M*M'
											  //	MatDense<FPP, DEVICE> sM_M, sMM_; // M'*M, M*M'


			dM = sM;
			spgemm(sM, dM, dM_M, FPP(1.0), FPP(0.0), 'H', 'N');
			spgemm(sM, dM, dMM_, FPP(1.0), FPP(0.0), 'N', 'H');
			M = &sM;
			//	svdtj_core_cplx(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
			svdtj_core_gen<FPP,DEVICE,FPP2>(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj_core_gen_blind(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
		{

			assert(order != 0);  // we need to use a specific order of singular values/vectors approximates to know where to get the singular values in W1' M W2 matrix

			GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW1;
			GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW2;

			if (t1 <= 0)
				//automatic t1
				t1 = dMM_.getNbRow() / 2;
			if (t2 <= 0)
				t2 = dM_M.getNbRow() / 2;

			instantiate_algos(&algoW1, &algoW2, dM_M, dMM_, J1, J2, t1, t2, verbosity, tol, relErr, enable_large_Faust);

			try {

				//TODO: parallelize with OpenMP ?
				algoW1->compute_facts();
				algoW2->compute_facts();

				auto m = M->getNbRow();
				auto n = M->getNbCol();
				auto min_mn = m > n?n:m;

				Vect<FPP,DEVICE> S(min_mn); // zero-initialized


				Transform<FPP,DEVICE> transW1 = std::move(algoW1->get_transform(order)); // 0 for no order, -1 descending order


				Transform<FPP,DEVICE> transW2 = std::move(algoW2->get_transform(order));


				// compute S = W1'*M*W2 = W1'*(W2'*M')'
				dM.adjoint();
				MatDense<FPP, DEVICE> MW2 = transW2.multiply(dM, 'H');
				MW2.adjoint();
				MatDense<FPP,DEVICE> W1_MW2 = transW1.multiply(MW2, 'H');

				// set singular values
				if(order < 0)
					for(int i=0;i<min_mn;i++) //TODO: OpenMP?
						S[i] = W1_MW2(i,i);
				else if (order >= 0)
					for(int i=0;i<min_mn;i++) //TODO: OpenMP?
						S[i] = W1_MW2(m - min_mn + i,n - min_mn + i);

#if DEBUG_SVDTJ
				MatDense<FPP, Cpu> S1_mat(S.size(), 1, S.getData());
				S1_mat.save_to_mat_file("/tmp/S1_cpp.mat", "S1_cpp");
#endif


				TransformHelper<FPP,DEVICE> *thW1 = new TransformHelper<FPP,DEVICE>(transW1, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)

				TransformHelper<FPP,DEVICE> *thW2 = new TransformHelper<FPP,DEVICE>(transW2, true); // the same

#if DEBUG_SVDTJ
				dM.save_to_mat_file("/tmp/M.mat", "M");
				thW1->save_mat_file("/tmp/W1_cpp.mat");
				thW2->save_mat_file("/tmp/W2_cpp.mat");

				MW2.save_to_mat_file("/tmp/MW2_cpp.mat", "MW2_cpp");

				W1_MW2.save_to_mat_file("/tmp/W1_MW2_cpp.mat", "W1_MW2_cpp");
#endif
				svdtj_sign_W1_S(m, n, order, S, S_, *thW1);

				*U = thW1;
				*V = thW2;

			}
			catch(out_of_range e)
			{
				// algorithm stopped before the first iteration (nGivens too large)
			}

			delete algoW1;
			delete algoW2;

		}

	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth1(MatDense<FPP,DEVICE> &dM, bool &dM_adj, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2)
		{

			// compute S = W1'*M*W2 = W1'*(W2'*M')'
			// do it the simplest way: a) compute MW2 = (W2' * M')' as a dense matrix then W1'* MW2 in a second matrix W1_MW2 (the one returned)
			// COST:
			// memory cost: more expensive than meth2 and meth3, //TODO: verify precisely
			// computation cost: about the same cost than meth2 and slower than meth3
			if(! dM_adj)
			{
				dM.adjoint();
				dM_adj = true;
			}
			MatDense<FPP, DEVICE> MW2 = tW2.multiply(dM, 'H');
			MW2.adjoint();
			MatDense<FPP,DEVICE> W1_MW2 = tW1.multiply(MW2, 'H');
			return W1_MW2;
		}

	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth2(MatDense<FPP,DEVICE> &dM, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2)
		{
			// compute S = W1'*M*W2 = W1'*(W2'*M')'
			// the main idea is to build two "Fausts" (TransformHelper) one for the expression F =(W2'*M')', the other for  G = W1'
			// then W1'*M*W2 is compute in two steps M = F.get_product() then G * M' gives S
			// NOTE: we don't use a single Faust for the whole S expression because some factors are in adjoint state, the other not, so it can't be set in one Faust whose the adjoint state is uniform (not set by factor)
			// NOTE: this method has no reason to be faster than *meth1 but is a first step forward *meth3.
			// However, note that it avoids the first dM.adjoint()
			std::vector<MatGeneric<FPP, DEVICE>*> W1_facts(tW1.begin(), tW1.end());
			std::vector<MatGeneric<FPP, DEVICE>*> W2_facts(tW2.begin(), tW2.end());
			TransformHelper<FPP, DEVICE> thW1(W1_facts, FPP(1.0), /* no copy*/ false, false);
			auto thW1_H = thW1.adjoint(/*inplace*/true); // thW1 same obj as thW1_H
			TransformHelper<FPP, DEVICE> thM_W2(W2_facts, FPP(1.0), /* no copy*/ false, false);
			thM_W2.push_first(&dM, /* no copy */ false, false);
			auto thM_W2_H = thM_W2.adjoint(/*inplace*/true); // thM_W2_H same obj as thM_W2

			// disable dtors because we use directly tW1 and tW2 data, they are responsible to free the memoory
			thW1.disable_dtor();
			thW1_H->disable_dtor();
			thM_W2.disable_dtor();
			thM_W2_H->disable_dtor();

			auto W1_MW2 = thM_W2_H->get_product();
			W1_MW2.adjoint();
			W1_MW2 = thW1_H->multiply(W1_MW2);

			return W1_MW2;
		}


	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth3(const MatDense<FPP,DEVICE> &dM, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2, MatDense<FPP, DEVICE> & prev_W1H_M_W2, const int err_step_period, const int k1, const int k2, const int t1, const int t2, const bool new_W1, const bool new_W2)
		{
			int id_shift1, id_shift2;

			bool first_call = prev_W1H_M_W2.getNbRow() == 0;
			if(first_call)
			{
				// this is the first time we compute W1' M W2
				prev_W1H_M_W2 = dM;
				id_shift2 = tW2.size();
				id_shift1 = tW1.size();
			}
			else
			{
				// W1' M W2 has already been computed once (at least)
				// take into account the last permuted factor of W1 (resp. W2)
				// it has been computed in the previous call of this function and must be retrieved here again because it won't be the one permuted here
				// (it's the reason why there is a +1 in shifts)
				auto n1 = k1 / t1 % err_step_period;
				auto n2 = k2 / t2 % err_step_period;
				if (n1)
				{
					// this function call is not due to a periodic error verification
					// it might be because of the stopping criterion (e.g. J1 reaching)
					// compute the number factors to consider since last call
					id_shift1 = n1 + 1;
				}
				else
					id_shift1 = err_step_period + 1;
				// do the same for W2
				if(n2)
					id_shift2 = n2 + 1;
				else
					id_shift2 = err_step_period + 1;

			}

			assert(new_W1 || new_W2 || first_call); // it does'nt make sense to compute W1H_M_W2 again if neither W1 or W2 have changed
																	   // except if this is the first time this function is called

			// NOTE1: disable dtors below because we use directly tW1 and tW2 data
			// they are responsible to free the memoory themselves

			TransformHelper<FPP, DEVICE>* thW1_H = nullptr;
			TransformHelper<FPP, DEVICE>* thM_W2_H = nullptr;

			if(new_W1 || first_call)
			{
				std::vector<MatGeneric<FPP, DEVICE>*> W1_facts(tW1.end() - id_shift1, tW1.end() - 1); // works only if err_step_period >= 3
				TransformHelper<FPP, DEVICE> thW1(W1_facts, FPP(1.0), /* no copy*/ false, false);
				thW1.disable_dtor(); // cf. NOTE1
				thW1_H = thW1.adjoint(/*inplace*/ false); // no inplace because thW1 is in local scope
				thW1_H->disable_dtor();
			}

			if(new_W2 || first_call)
			{
				std::vector<MatGeneric<FPP, DEVICE>*> W2_facts(tW2.end() - id_shift2, tW2.end() - 1);
				TransformHelper<FPP, DEVICE> thM_W2(W2_facts, FPP(1.0), /* no copy*/ false, false);
				thM_W2.disable_dtor(); //cf. NOTE1
				thM_W2.push_first(&prev_W1H_M_W2, /* no copy */ false, false);
				thM_W2_H = thM_W2.adjoint(/*inplace*/ true);
				thM_W2_H->disable_dtor();
				prev_W1H_M_W2 = thM_W2_H->get_product();
				prev_W1H_M_W2.adjoint();
			}

			if(new_W1 || first_call) // thW1_H not nullptr
			{
				prev_W1H_M_W2 = thW1_H->multiply(prev_W1H_M_W2);
				delete thW1_H;
			}

			MatDense<FPP, DEVICE> W1_MW2 = prev_W1H_M_W2;

			// now multiply by W1 and W2 last permuted factors to get full W1_MW2
			auto lastW1 = dynamic_cast<MatSparse<FPP, DEVICE>*>(*(tW1.end() - 1));
			lastW1->multiply(W1_MW2, 'H');

			auto lastW2 = dynamic_cast<MatSparse<FPP, DEVICE>*>(*(tW2.end() - 1));
			W1_MW2.multiplyRight(*lastW2);

			return W1_MW2;
		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj_core_gen_step(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_)
		{

			// algorithm main idea: run two eigtj for U and V step by step until the target error is verified or the limit number of Givens is reached

			assert(order != 0);  // we need to use a specific order of singular values/vectors approximates to know where to get the singular values in W1' M W2 matrix
			GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW1;
			GivensFGFTGen<Real<FPP>, DEVICE, FPP2, FPP>* algoW2;

			if (t1 <= 0)
				//automatic t1
				t1 = dMM_.getNbRow() / 2;

			if (t2 <= 0)
				t2 = dM_M.getNbRow() / 2;

			instantiate_algos(&algoW1, &algoW2, dM_M, dMM_, /*J*/ 1, 1, t1, t2, verbosity, /*tol*/ FPP2(0), relErr, enable_large_Faust);
			// J == 1 to prevent exception raising in eigtj when J == 0 and tol == 0
			// tol == 0 because the error criterion is handled in this function if it is the stopping criterion (otherwise this is the iteration number determined by J and t)

#define continue_loop() \
			((J1 == 0 || k1 < J1) && (J2 == 0 || k2 < J2)) && (tol == FPP2(0) || tol < err) && (new_W1 || new_W2)
			// three scenarii:
			// 1) stopping only on error criterion : J == 0, 0 < tol (<= 1), stop if err <= tol
			// 2) stopping only when a certain number of Givens are built: J > 0, t Givens built per iteration/factor for both U/W1 and V/W2, stop if k (the number of Givens built on each side) is greater or equal to J
			// 3) Concurrent 1) and 2), the most restrictive criterion applies

			try {

				auto m = M->getNbRow();
				auto n = M->getNbCol();
				auto min_mn = m > n?n:m;

				Vect<FPP,DEVICE> S(min_mn); // zero-initialized
				int k1 = 0, k2 = 0; // k counts the number of Givens built (t per factor) on each side (for U (W1) and V (W2))
				auto M_norm = dM.norm();
				bool dM_adj = false; //cf. svdtj_compute_W1H_M_W2_meth1
				Real<FPP> err = 1; // relative error of the SVDTJ on norms
				Transform<FPP, DEVICE> tW1, tW2;
				bool loop = true;
				auto err_step_period = 100; //TODO: it should be a user parameter one day, (cf. svdtj_compute_W1H_M_W2_meth)
				MatDense<FPP, DEVICE> W1_MW2_; // previous computation of W1' M W2 (initialization at M; no Givens yet, in svdtj_compute_W1H_M_W2_meth3)
				MatDense<FPP, DEVICE> W1_MW2;

				bool better_W1 = true, better_W2 = true;

				auto W1_last_err = FPP2(-1);
				auto W2_last_err = FPP2(-1);

				auto W1_err = FPP2(-1);
				auto W2_err = FPP2(-1);

				auto W1_max_size = std::numeric_limits<int>::max();
				auto W2_max_size = std::numeric_limits<int>::max();
				bool new_W1 = true, new_W2 = true; // true if W1 (resp. W2) has grown up during the iteration

				bool W1_too_long = false, W2_too_long = false; // used only if enable_large_Faust is false
				if(! enable_large_Faust)
				{
					W1_max_size = m * m / 4;
					W2_max_size = n * n / 4;
				}

				assert(tol > FPP2(0) || J1 > 0 && J2 > 0);
				assert(J1 == 0 || J1 > t1);
				assert(J2 == 0 || J2 > t2);
//				assert(J1 == 0 || J1 > t1);
//				assert(J2 == 0 || J2 > t2);
				while(loop) // if the error is not a stop crit then tol is normally 0
				{

					// make eigtj steps for W1/U W2/V
					// but don't continue the eigtj stopped to enhance the eigenvalues accuracy
					if(k1 == 0 || new_W1 && better_W1 && ! W1_too_long) // k1 == 0°<=> k2 == 0
						algoW1->next_step();
					else
						new_W1 = false;

					if(k1 == 0 || new_W2 && better_W2 && ! W2_too_long)
						algoW2->next_step();
					else
						new_W2 = false;

					if(! enable_large_Faust && ! W1_too_long)
					{
						W1_too_long = algoW1->nfacts() >= W1_max_size;
						if(W1_too_long)
							std::cerr << "Warning: eigtj stopped on U approximate which is too long to be worth it, set enable_large_Faust to true if you want to continue anyway..." << std::endl;
					}

					if(! enable_large_Faust && ! W2_too_long)
					{
						W2_too_long = algoW2->nfacts() >= W2_max_size;
						if(W2_too_long)
							std::cerr << "Warning: eigtj stopped on V approximate which is too long to be worth it, set enable_large_Faust to true if you want to continue anyway..." << std::endl;
					}

					k1 += t1;
					k2 += t2;
					loop = continue_loop();

					bool verif_err_ite = ! (k1 % (err_step_period * t1));

					if(verif_err_ite)
					{
						if(better_W1)
							W1_err = algoW1->calc_err();

						if(better_W2)
							W2_err = algoW2->calc_err();

						// check the errors obtained so far
						if(better_W1 && W1_last_err >= 0)
						{
							better_W1 = W1_last_err > W1_err;
							if(verbosity > 0 && ! better_W1)
								std::cout << "Warning: U approximate accuracy is not increasing anymore, stopped eigtj tunning for U." << std::endl;
						}

						if(better_W2 && W2_last_err >= 0)
						{
							better_W2 = W2_last_err > W2_err;
							if(verbosity > 0 && ! better_W2)
								std::cout << "Warning: V approximate accuracy is not increasing anymore, stopped eigtj tunning for V." << std::endl;
						}

						if(better_W1)
							W1_last_err = W1_err;

						if(better_W2)
							W2_last_err = W2_err;

						if(new_W1)
							new_W1 = better_W1;
						if(new_W2)
							new_W2 = better_W2;

						// if both U and V are not enhancing, it's time to stop the loop
						loop = continue_loop();
					}

					if(tol != FPP2(0) && verif_err_ite || ! loop)
					{
						// we have either to verify the error of approximate or to finish the algorithm (compute S based on W1 and W2)
						tW1 = std::move(algoW1->get_transform(/*order*/ order, /* copy */ false, /*nfacts*/ better_W1?-1:algoW1->nfacts() - err_step_period));
						tW2 = std::move(algoW2->get_transform(/* order*/ order, /* copy */ false, /*nfacts*/ better_W2?-1:algoW2->nfacts() - err_step_period));
						// compute S = W1'*M*W2 = W1'*(W2'*M')'

						//				//	method 1
						//				auto W1_MW2 = svdtj_compute_W1H_M_W2_meth1(dM, dM_adj, tW1, tW2);

						//				// method 2
						//				auto W1_MW2 = svdtj_compute_W1H_M_W2_meth2(dM, tW1, tW2);

						// method 3
						if(new_W1 || new_W2 || ! W1_MW2_.getNbRow())
						{
							W1_MW2 = svdtj_compute_W1H_M_W2_meth3(dM, tW1, tW2, W1_MW2_, err_step_period, k1, k2, t1, t2, new_W1, new_W2);

							if(order < 0)
								for(int i=0;i<min_mn;i++) //TODO: OpenMP?
									S[i] = W1_MW2(i,i);
							else if (order >= 0)
								for(int i=0;i<min_mn;i++) //TODO: OpenMP?
									S[i] = W1_MW2(m - min_mn + i,n - min_mn + i);

							// compute error
							auto Sd_norm = S.norm();
							err = (M_norm - Sd_norm);
							if(relErr) err /= M_norm;
							// erase permuted factors because algoW1/W2 doesn't keep account for them (cf. GivensFGFTGen::get_transform with ord == true and copy == false
							if(verbosity)
							{
								std::cout << "computing the error at iteration k: " << int(k1 / t1) << " t1, t2: " << t1 << ", " << t2 << std::endl;
								std::cout << "norm err: " << err << std::endl;
							}

						}
						if(order)
						{
							tW1.erase(tW1.size()-1);
							tW2.erase(tW2.size()-1);
						}
						// re-evaluate the stopping criterion (because of the error update)
						loop = continue_loop();
					}
				}

				// don't use tW1 and tW2 here because we want an independent copy (not linked to algoW1 and algoW2 internal data)
				Transform<FPP,DEVICE> transW1 = std::move(algoW1->get_transform(/*order*/ order, /*copy*/ true, /*nfacts*/ better_W1?-1:algoW1->nfacts() - err_step_period));

				Transform<FPP,DEVICE> transW2 = std::move(algoW2->get_transform(/*order*/ order, /*copy*/ true, /*nfacts*/ better_W2?-1:algoW2->nfacts() - err_step_period));

				TransformHelper<FPP,DEVICE> *thW1 = new TransformHelper<FPP,DEVICE>(transW1, true); // true is for moving and not copying the Transform object into TransformHelper (optimization possible cause we know the original object won't be used later)

				TransformHelper<FPP,DEVICE> *thW2 = new TransformHelper<FPP,DEVICE>(transW2, true); // the same

#if DEBUG_SVDTJ
				dM.save_to_mat_file("/tmp/M.mat", "M");
				thW1->save_mat_file("/tmp/W1_cpp.mat");
				thW2->save_mat_file("/tmp/W2_cpp.mat");

				W1_MW2.save_to_mat_file("/tmp/W1_MW2_cpp.mat", "W1_MW2_cpp");

				MatDense<FPP, Cpu> S1_mat(S.size(), 1, S.getData());
				S1_mat.save_to_mat_file("/tmp/S1_cpp.mat", "S1_cpp");
#endif

				svdtj_sign_W1_S(m, n, order, S, S_, *thW1);

				*U = thW1;
				*V = thW2;

			}
			catch(out_of_range e)
			{
				// algorithm stopped before the first iteration (nGivens too large)
			}

			delete algoW1;
			delete algoW2;

		}

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void svdtj_core_gen(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J1, int J2, int t1, int t2, FPP2 tol, unsigned int verbosity, bool relErr, int order /* ignored anyway, kept here just in case */, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_, const bool by_step/*=true*/)
		{
			if(by_step)
				svdtj_core_gen_step(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
			else
				svdtj_core_gen_blind(M, dM, dM_M, dMM_, J1, J2, t1, t2, tol, verbosity, relErr, order, enable_large_Faust, U, V, S_);
		}
}
