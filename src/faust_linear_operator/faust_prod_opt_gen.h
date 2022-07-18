#ifndef __FAUST_PROD_OPT_GEN__
#define __FAUST_PROD_OPT_GEN__
namespace Faust
{
	/**
	 * Multipy the matrix chain "factors" using the dynamic programming method to determine the order (parenthesis positions) in the product (reducing the computation cost).
	 */
	template<typename FPP>
		MatDense<FPP, Cpu> dynprog_multiply(std::vector<MatGeneric<FPP, Cpu>*>& factors, const char op='N', const MatGeneric<FPP, Cpu>* A=nullptr);


	int cost_bsr_dense(int A_bnnz, int A_bm, int A_bn, int B_nrows, int B_ncols, bool B_transp);


	int cost_bsr_sparse(int A_bnnz, int A_bm, int A_bn, bool A_transp, int B_nnz);


	int cost_bsr_bsr(int A_bnnz, int A_bm, int A_bn, bool A_transp, int B_nnz);


	int cost_sparse_dense(int A_nnz, int B_nrows, int B_ncols, bool B_transp);


	int cost_sparse_sparse(int A_nrows, int A_ncols, int A_nnz, bool A_transp, int B_nnz);


	int cost_sparse_bsr(int A_nnz, int B_bnnz, int B_bm, int B_bn, bool B_transp);


	int cost_dense_dense(int A_nrows, int A_ncols, int B_nrows, int B_ncols, bool B_transp);


	int cost_dense_sparse(int A_nrows, int A_ncols, bool A_transp, int B_nnz);


	int cost_dense_bsr(int A_nrows, int A_ncols, bool A_transp, int B_bnrows, int B_bncols, int B_bnnz);
}
#include "faust_prod_opt_gen.hpp"
#endif
