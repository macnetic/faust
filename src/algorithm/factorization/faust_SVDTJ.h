#ifndef __FAUST_SVDTJ__
#define __FAUST_SVDTJ__
#include "faust_MatDense.h"
#include "faust_MatSparse.h"

#include "faust_constant.h"
namespace Faust
{

	//TODO document the prototypes
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj(MatDense<FPP, DEVICE> & M, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S);
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj(MatSparse<FPP, DEVICE> & M, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S);
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj_cplx(MatDense<FPP, DEVICE> & M, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S);
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj_cplx(MatSparse<FPP, DEVICE> & M, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust,TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S);


	/**
	 * \brief Wrapper for blind and step versions of SVDTJ.
	 */
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj_core_gen(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_, const bool by_step=true);

	/**
	 * This version runs two eigtjs step by step until the tol error or the number of Givens J is reached.
	 */
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj_core_gen_step(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_);

	/**
	 * This version runs two eigtjs blindly until they reach the tol error for the eigen decomposition or the number of Givens is reached.
	 */
	template<typename FPP, FDevice DEVICE, typename FPP2 = float>
		void svdtj_core_gen_blind(MatGeneric<FPP,DEVICE>* M, MatDense<FPP,DEVICE> &dM, MatDense<FPP,DEVICE> &dM_M, MatDense<FPP,DEVICE> &dMM_, int J, int t, FPP2 tol, unsigned int verbosity, bool relErr, int order, const bool enable_large_Faust, TransformHelper<FPP,DEVICE> ** U, TransformHelper<FPP,DEVICE> **V, Vect<FPP,DEVICE> ** S_);

	template<typename FPP, FDevice DEVICE, typename FPP2>
		void instantiate_algos(GivensFGFTGen<Real<FPP>, Cpu, FPP2, FPP>** algoW1, GivensFGFTGen<Real<FPP>, Cpu, FPP2, FPP>** algoW2, Faust::MatDense<FPP,DEVICE> &dM_M, Faust::MatDense<FPP,DEVICE> &dMM_, int J1, int J2, int t1, int t2, unsigned int verbosity, FPP2 tol, bool relErr, bool enable_large_Faust);

	/**
	 * \brief Move signs of S in matching singular vector of W1 (U).
	 *
	 * \param m: the number of rows of the matrix we compute SVD.
	 * \param n: the number of columns of the matrix we compute SVD.
	 * \param order: the order according to S has already been sorted (before call). The affected W1 "vectors" depend of that order.
	 * \param S: the ordered singular values.
	 * \param S_: the output for absolute singular values.
	 * \param thW1: the input/output of the left singular vectors (W1 / U).
	 */
	template<typename FPP, FDevice DEVICE>
		void svdtj_sign_W1_S(const faust_unsigned_int &m, const faust_unsigned_int &n, const int order, const Vect<FPP, DEVICE> &S, Vect<FPP, DEVICE> **S_, TransformHelper<FPP, DEVICE> &thW1);

	/**
	 * Utility function of svdtj_core_gen_step for computing U'MV (or W1' M W2), method 1 (the simplest way).
	 *
	 * This method is not used but is kept as reference because the code goes more an more complicated for methods 2 and 3.
	 *
	 * \param dM: (in out) matrix of which we calculate the SVD approximate as a MatDense.
	 * \param dM_adj: (in out) true if dM is already the adjoint of the original matrix, false otherwise (in which case the adjoint is computed internvally and dM is modified). The goal is to avoid computing the adjoint M many times in svdtj_core_gen_step.
	 * \param tW1: the "Faust" of Givens matrices for U (W1).
	 * \param tW2: the "Faust" of Givens matrices for V (W2).
	 *
	 * \return the MatDense for W1'*M*W2.
	 */
	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth1(MatDense<FPP,DEVICE> &dM, bool &dM_adj, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2);

	/**
	 * Utility function of svdtj_core_gen_step for computing U'MV (or W1' M W2), method 2 (not faster than method 1, just a sketching of method 3).
	 *
	 * This method is not used but is kept as reference because the code goes more an more complicated in method 3.
	 *
	 * \param dM: matrix of which we calculate the SVD approximate as a MatDense.
	 * \param tW1: the "Faust" of Givens matrices for U (W1).
	 * \param tW2: the "Faust" of Givens matrices for V (W2).
	 *
	 * \return the MatDense for W1'*M*W2.
	 */
	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth2(const MatDense<FPP,DEVICE> &dM, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2);


	/**
	 * Utility function of svdtj_core_gen_step for computing U'MV (or W1' M W2), method 3 (fastest method).
	 *
	 * It is advised to refer to *meth1 and *meth2 function in order to understand this one.
	 *
	 * The principle of this method is to compute W1H_M_W2 recursively reusing previous computations.
	 * if W1 is a sequence of factors W1,0 W1,1 ... W1,N,
	 * W2 a sequence of factors W2,0 W2,1 ... W2,N
	 * S_j the successive products U' M V along the SVDTJ algorithm execution.
	 * Then S_{n} = W1,N' * S_{n-1} * W2,N.
	 *
	 *
	 * \param dM: matrix of which we calculate the SVD approximate as a MatDense.
	 * \param tW1: the "Faust" of Givens matrices for U (W1).
	 * \param tW2: the "Faust" of Givens matrices for V (W2).
	 * \param prev_W1H_M_W2: previous result of W1'*M*W2 to base the computation of this one.
	 * \param err_step_period: the period in number of Givens according which the SVDJT computes the error. It matters to determine the order of the recursive calculation of W1' M W2.
	 * \param k1, k2: number of Givens used for W1 and W2 until the latest iteration of SVDTJ.
	 * \param t1, t2: the number of Givens matrix in each W1 W2 factors.
	 *
	 * \return the MatDense for W1'*M*W2.
	 */
	template<typename FPP, FDevice DEVICE>
		MatDense<FPP,DEVICE> svdtj_compute_W1H_M_W2_meth3(const MatDense<FPP,DEVICE> &dM, const Transform<FPP, DEVICE> &tW1, const Transform<FPP, DEVICE> &tW2, MatDense<FPP, DEVICE> & prev_W1H_M_W2, const int err_step_period, const int k1, const int k2, const int t1, const int t2, const bool new_W1, const bool new_W2);

}

#include "faust_SVDTJ.hpp"
#endif
