#include "faust_prod_opt.h"

namespace Faust
{

	int cost_bsr_dense(int A_bnnz, int A_bm, int A_bn, int B_nrows, int B_ncols, bool B_transp)
	{
		auto c = A_bnnz * A_bm * A_bn;
		if(B_transp)
			return  c * B_nrows;
		else
			return c * B_ncols;
	}


	int cost_bsr_sparse(int A_bnnz, int A_bm, int A_bn, bool A_transp, int B_nnz)
	{
		if(A_transp)
			return A_bnnz*A_bn*B_nnz;
		else
			return A_bnnz*A_bm*B_nnz;
	}


	int cost_bsr_bsr(int A_bnnz, int A_bm, int A_bn, bool A_transp, int B_nnz)
	{
		auto cost_tosparse = A_bnnz * A_bm * A_bn; //TODO: very unsure (maybe too high compared to the mul cost)
		return cost_tosparse + cost_bsr_sparse(A_bnnz, A_bm, A_bn, A_transp, B_nnz);
	}


	int cost_sparse_dense(int A_nnz, int B_nrows, int B_ncols, bool B_transp)
	{
		if(B_transp)
			return A_nnz*B_nrows;
		else
			return A_nnz*B_ncols;
	}


	int cost_sparse_sparse(int A_nrows, int A_ncols, int A_nnz, bool A_transp, int B_nnz)
	{
		if(A_transp)
			return A_nnz*B_nnz/A_ncols;
		else
			return A_nnz*B_nnz/A_nrows;
	}


	int cost_sparse_bsr(int A_nnz, int B_bnnz, int B_bm, int B_bn, bool B_transp)
	{
		if(B_transp)
			return A_nnz * B_bnnz * B_bn;
		else
			return A_nnz * B_bnnz * B_bm;
	}


	int cost_dense_dense(int A_nrows, int A_ncols, int B_nrows, int B_ncols, bool B_transp)
	{
		int cost = A_nrows*A_ncols; // whether A is transp or not
		if(B_transp)
			cost *= B_nrows;
		else
			cost *= B_ncols;
		return cost;
	}


	int cost_dense_sparse(int A_nrows, int A_ncols, bool A_transp, int B_nnz)
	{
		if(A_transp)
			return A_ncols*B_nnz;
		return A_nrows*B_nnz;
	}


	int cost_dense_bsr(int A_nrows, int A_ncols, bool A_transp, int B_bnrows, int B_bncols, int B_bnnz)
	{
		int cost = B_bnnz * B_bnrows * B_bncols;
		if(A_transp)
			cost *= A_nrows;
		else
			cost *= A_ncols;
		return cost;
	}
}
