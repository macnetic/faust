#ifdef USE_GPU_MOD
#include "faust_linear_algebra_gpu.h"
#endif
namespace Faust
{
	template<typename FPP, FDevice DEV>
		MatGeneric<FPP, DEV>* dynprog_multiply_rec(const std::vector<MatGeneric<FPP, DEV>*>& factors, int** inds, int i, int j, const char op='N', const char last_op='N')
		{
			int j_minus_i = j-i;
			int p_nrows, p_ncols; // prod numbers rows, cols
			char op_left, op_right;
			if(j_minus_i == 0)
				return factors[i];
			else if(j_minus_i == 1)
			{
				op_left = op; // i < factors.size()-1 because j_minus_i == 1
				if(j < factors.size()-1)
					op_right = op;
				else // j == last factor
					op_right = last_op;
				if(op_left == 'N')
					p_nrows = factors[i]->getNbRow();

				else
					p_nrows = factors[i]->getNbCol();
				if(op_right == 'N')
					p_ncols = factors[j]->getNbCol();
				else
					p_ncols = factors[j]->getNbRow();
				auto prod = new MatDense<FPP, DEV>(p_nrows, p_ncols);
//				std::cout << factors[i]->getNbRow() << "x" << factors[i]->getNbCol() << "*" << factors[j]->getNbRow() << "x" << factors[j]->getNbCol() << std::endl;
				gemm_gen(*factors[i], *factors[j], *prod, (FPP)1.0, (FPP)0.0, op_left, op_right);
				return prod;
			}
			else
			{
				auto k = inds[i][j];


				auto prod1 = dynprog_multiply_rec(factors, inds, i, k, op, last_op);
				auto prod2 = dynprog_multiply_rec(factors, inds, k+1, j, op, last_op);

				if(i == k)
					op_left = op;
				else // prod1 was computed here, it's not one of factors[i], so no need to take care of op
					op_left = 'N';
				if(k+1 == j)
					if(j == factors.size()-1)
						op_right = last_op;
					else
						op_right = op;
				else // prod2 was computed here, it's not one of factors[i], so no need to take care of op
					op_right = 'N';

				if(op_left == 'N')
					p_nrows = prod1->getNbRow();
				else
					p_nrows = prod1->getNbCol();

				if(op_right == 'N')
					p_ncols = prod2->getNbCol();
				else
					p_ncols = prod2->getNbRow();

				auto prod12 = new MatDense<FPP, DEV>(p_nrows, p_ncols);
//				std::cout << prod1->getNbRow() << "x" << prod1->getNbCol() << "*" << prod2->getNbRow() << "x" << prod2->getNbCol() << std::endl;
				gemm_gen(*prod1, *prod2, *prod12, FPP(1.0), FPP(0.0), op_left, op_right);
				//TODO/ verify if a memory leak exists
				// delete matrices allocated in this function's recursive calls
				if(k-i > 1)
					delete prod1;
				if(j-k-1 > 1)
					delete prod2;
				return prod12;
			}
		}


	template<typename FPP, FDevice DEV>
		MatDense<FPP, DEV> dynprog_multiply(std::vector<MatGeneric<FPP, DEV>*>& factors, const char op/*='N'*/, const MatGeneric<FPP, DEV>* A/*=nullptr*/)
		{
			// manage useless cases when factors.size is too small
			if(factors.size() == 1)
			{
				MatDense<FPP, DEV> P;
				if(A != nullptr)
				{
					gemm_gen(*factors[0], *A, P, (FPP)1.0, (FPP)0.0, op, 'N');
					return P;
				}
				else
				{
					auto sp_mat = dynamic_cast<MatSparse<FPP, DEV>*>(factors[0]);
					if(sp_mat)
						return MatDense<FPP,DEV>(*sp_mat);
					else
						*dynamic_cast<MatDense<FPP, DEV>*>(factors[0]);
				}
			}
			char last_op = op;
			if(A != nullptr)
			{
				factors.push_back(const_cast<MatGeneric<FPP, DEV>*>(A)); // A won't be touched
				last_op = 'N';
			}
			// this function initializes a triplet of boolean depending on fac concrete type (MatDense, MatSparse or MatBSR)
			auto init_fac_type_bools = [](const MatGeneric<FPP, DEV>* fac, bool& fac_dense, bool &fac_sparse, bool &fac_bsr)
			{
				fac_dense = fac_sparse = fac_bsr = false;
				std::runtime_error et("dynprog_multiply error: non-supported matrix type (only MatDense, MatSparse, MatBSR are)");
				if(! (fac_dense = dynamic_cast<const MatDense<FPP, DEV>*>(fac)))
					if(! (fac_sparse = dynamic_cast<const MatSparse<FPP, DEV>*>(fac)))
						if(! (fac_bsr = dynamic_cast<const MatBSR<FPP, DEV>*>(fac)))
							throw et;
			};
			const int n = factors.size();
			int** c = new int*[n]; // TODO: reduce the memory used because only the upper triangle of the array is used
			int** inds = new int*[n]; // TODO: idem
			int j, k, cost;
			int c_i, c_j;
			MatBSR<FPP, DEV> *lf_bsr_mat, *rf_bsr_mat;
			for(int i=0;i<n;i++)
			{
				c[i] = new int[n];
				inds[i] = new int[n];
				c[i][i] = 0;
			}
			for(int p=1;p<n;p++)
				for(int i=0;i<n-p;i++)
				{
					j = p + i;
					// left factor is dense, sparse or bsr
					bool lf_dense, lf_sparse, lf_bsr;
					// right factor is dense, sparse or bsr
					bool rf_dense, rf_sparse, rf_bsr;
					init_fac_type_bools(factors[i], lf_dense, lf_sparse, lf_bsr);
					init_fac_type_bools(factors[j], rf_dense, rf_sparse, rf_bsr);

					k = i;
					c[i][j] = std::numeric_limits<int>::max();
					auto lf_nrows = factors[i]->getNbRow();
					auto lf_ncols = factors[i]->getNbCol();
					auto rf_nrows = factors[j]->getNbRow();
					auto rf_ncols = factors[j]->getNbCol();
					while(k < j)
					{
						auto lf_init = k == i; // if true : left factor and middle factor are the same matrix which is an initial matrix (it doesn't result from intermediate product)
						auto rf_init = k + 1 == j; // true means that the right factor is composed of only one initial matrix (it doesn't result from intermediate product)
						cost = c[i][k] + c[k+1][j];
						// middle factor nrows/ncols
						auto mf_nrows = factors[k]->getNbRow();
						auto mf_ncols = factors[k]->getNbCol();

						if(lf_bsr) lf_bsr_mat = dynamic_cast<MatBSR<FPP,DEV>*>(factors[i]);
						if(rf_bsr) rf_bsr_mat = dynamic_cast<MatBSR<FPP,DEV>*>(factors[j]);
						if(lf_dense && rf_dense)
							// left and right factors are dense matrices
							cost += cost_dense_dense(lf_nrows, mf_ncols, rf_nrows, rf_ncols, op != 'N' && rf_init); // no need to consider op for product resulting dense matrices
						else if(lf_dense && rf_sparse)
							cost += cost_dense_sparse(lf_nrows, mf_ncols, op != 'N' && lf_init, factors[j]->getNonZeros());
						else if(lf_dense && rf_bsr)
							cost += cost_dense_bsr(lf_nrows, mf_ncols, op != 'N' && lf_init, rf_bsr_mat->getNbBlockRow(), rf_bsr_mat->getNbBlockCol(), rf_bsr_mat->getNBlocks());
						else if(lf_sparse && rf_sparse)
							// lf_sparse => lf_init
							cost += cost_sparse_sparse(lf_nrows, lf_ncols, factors[i]->getNonZeros(), op != 'N', factors[j]->getNonZeros());
						else if(lf_sparse && rf_dense /*! rf_sparse*/)
							cost += cost_sparse_dense(factors[i]->getNonZeros(), rf_nrows, rf_ncols, op != 'N' && rf_init);
						else if(lf_sparse && rf_bsr)
							// rf_bsr => rf_fac_initial
							cost += cost_sparse_bsr(factors[i]->getNonZeros(), rf_bsr_mat->getNBlocks(), rf_bsr_mat->getNbBlockRow(), rf_bsr_mat->getNbBlockCol(), op != 'N');
						else if(lf_bsr && rf_dense)
							// lf_bsr => lf_init
							cost += cost_bsr_dense(lf_bsr_mat->getNBlocks(), lf_bsr_mat->getNbBlockRow(), lf_bsr_mat->getNbBlockCol(), rf_nrows, rf_ncols, op != 'N');
						else if(lf_bsr && rf_sparse)
							// lf_bsr => lf_init
							cost += cost_bsr_sparse(lf_bsr_mat->getNBlocks(), lf_bsr_mat->getNbBlockRow(), lf_bsr_mat->getNbBlockCol(), op != 'N', factors[j]->getNonZeros());
						else if(lf_bsr && rf_bsr)
							// lf_bsr => lf_init
							cost += cost_bsr_bsr(lf_bsr_mat->getNBlocks(), lf_bsr_mat->getNbBlockRow(), lf_bsr_mat->getNbBlockCol(), op != 'N', factors[j]->getNonZeros());
						if(cost < c[i][j])
						{
							c[i][j] = cost;
							inds[i][j] = k;
						}
						k++;
					}
				}
			auto prod = dynamic_cast<MatDense<FPP, DEV>*>(dynprog_multiply_rec(factors, inds, 0, n-1, op, last_op));
			for(int i=0;i<n;i++)
			{
				delete[] c[i];
				delete[] inds[i];
			}
			delete[] c;
			delete[] inds;
			auto M = std::move(*prod);
			delete prod;
			if(A != nullptr)
				factors.erase(factors.end()-1);
			return std::move(M);
		}

}
