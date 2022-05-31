/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2021):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/

#include "faust_FFT.h"
#include "faust_WHT.h"
#include "faust_linear_algebra.h"
#include "faust_prod_opt.h"
#include <chrono>
#include <cstdlib>

namespace Faust {


	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
				const FPP lambda_, const bool optimizedCopy, const bool cloning_fact,
				const bool internal_call) : TransformHelper<FPP,Cpu>()
	{
		if(lambda_ != FPP(1.0) && ! internal_call)
			cerr << "WARNING: the constructor argument for multiplying the Faust by a scalar is DEPRECATED and might not be supported in next versions of FAµST." << endl;
		this->transform = make_shared<Transform<FPP,Cpu>>(facts, lambda_, optimizedCopy, cloning_fact);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : TransformHelperGen<FPP,Cpu>()
	{
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(Transform<FPP,Cpu> &t, const bool moving /* default to false */) : TransformHelper<FPP,Cpu>()
	{
		if(moving)
			this->transform = make_shared<Transform<FPP,Cpu>>(std::move(t));
		else
			this->transform = make_shared<Transform<FPP,Cpu>>(t);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const TransformHelper<FPP,Cpu>* th_left, const TransformHelper<FPP,Cpu>* th_right)
		: TransformHelper<FPP,Cpu>()
		{
			this->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(th_left)->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(th_left)->eval_fancy_idx_Transform();
			bool right_left_transposed = th_left->is_transposed && th_right->is_transposed;
			bool right_left_conjugate = th_left->is_conjugate && th_right->is_conjugate;
			if(right_left_transposed)
			{
				auto tmp = th_left;
				th_left = th_right;
				th_right = tmp;
			}
			this->transform = make_shared<Transform<FPP,Cpu>>(th_left->transform.get(), th_left->is_transposed && ! right_left_transposed, th_left->is_conjugate && ! right_left_conjugate,
					th_right->transform.get(), th_right->is_transposed && ! right_left_transposed, th_right->is_conjugate && ! right_left_conjugate);
			// if the both are transposed, the factors won't be transposed in the Transform underlying object,
			// for optimization just set the transpose flag here
			this->is_transposed = right_left_transposed;
			//likewise for the conjugate
			this->is_conjugate = right_left_conjugate;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate) : Faust::TransformHelper<FPP,Cpu>()
	{
		this->transform = th->transform;
		this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
		this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		this->copy_slice_state(*th, transpose);
		this->copy_fancy_idx_state(*th, transpose);
		copy_mul_mode_state(*th);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th): TransformHelper<FPP,Cpu>()
	{
		this->copy_state(*th);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>& th): TransformHelper<FPP,Cpu>()
	{
		*this = th;
	}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::operator=(TransformHelper<FPP,Cpu>& th)
		{
			this->copy_state(th);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]): TransformHelper<FPP,Cpu>()
	{
		this->init_sliced_transform(th, s);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols): TransformHelper<FPP,Cpu>()
	{
		this->init_fancy_idx_transform(th, row_ids, num_rows, col_ids, num_cols);
	}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::copy_mul_mode_state(const TransformHelper<FPP,Cpu>& th)
		{
			TransformHelperGen<FPP,Cpu>::copy_mul_mode_state(th);
		}

#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
	template<typename FPP>
		template<typename ... GList>
		TransformHelper<FPP,Cpu>::TransformHelper(GList& ... t) : TransformHelper<FPP,Cpu>()
		{
			this->push_back_(t...);
		}
#endif
	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatSparse<FPP,Cpu> &A)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			//TODO: refactor with multiply(MatDense)
			Faust::MatDense<FPP,Cpu> M;
#ifdef FAUST_TORCH
			if(this->mul_order_opt_mode >= TORCH_CPU_L2R && this->mul_order_opt_mode <= TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R && !tensor_data.size())
			{
				// init tensor data cache
				convMatGenListToTensorList(this->transform->data, tensor_data, at::kCPU, /* clone */ false, /* transpose */ ! this->is_transposed);
//				display_TensorList(tensor_data);
			}
#endif

			switch(this->mul_order_opt_mode)
			{
				case GREEDY_ALL_ENDS:
				case GREEDY_1ST_BEST:
				case GREEDY_ALL_BEST_CONVDENSE:
				case GREEDY_ALL_BEST_GENMAT:
					{
						std::vector<Faust::MatGeneric<FPP,Cpu>*> data(this->transform->data);
						data.resize(data.size()+1); // for the matrix multiplying Faust
						data[data.size()-1] = const_cast<Faust::MatSparse<FPP,Cpu>*>(&A);
						std::vector<char> transconj_flags = {'N'}; // 'N' for all Faust factors and the matrix A
						if(this->is_transposed)
						{
							transconj_flags = std::vector<char>(data.size(), this->isTransposed2char()); // TODO: use assign in C++ 20
							*(transconj_flags.end()-1) = 'N'; // the matrix multiplyed is not transposed
							std::reverse(data.begin(), data.end()-1); // reversing the Faust factors (because they are transposed or transconjugate)
						}
						Faust::multiply_order_opt(this->mul_order_opt_mode, data, M, /*alpha */ FPP(1.0), /* beta */ FPP(0.0), transconj_flags);
					}
					break;
				case DYNPROG:
					{
						std::vector<Faust::MatGeneric<FPP,Cpu>*> data = this->transform->data;
						if(this->is_transposed)
							std::reverse(data.begin(), data.end());
						M = std::move(dynprog_multiply(data, this->isTransposed2char(), &A));
					}
					break;
				case CPP_PROD_PAR_REDUC:
				case OMP_PROD_PAR_REDUC:
					throw std::runtime_error("CPP_PROD_PAR_REDUC and OMP_PROD_PAR_REDUC are not capable to handle Faust-MatSparse product (only Faust-MatDense product is available).");
					break;
#ifdef FAUST_TORCH
				case TORCH_CPU_L2R:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_GREEDY:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false,  /*clone */ false,/* chain_opt */ true, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ true, !this->is_transposed);
					break;
#endif
				default:
					M = this->transform->multiply(A, this->isTransposed2char());
					break;
			}
			return std::move(M);
		}


	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> &x)
		{
			// the prototypes below use this function
			Vect<FPP,Cpu> v;
			if(this->is_sliced)
				return this->multiply(x.getData()); // this->is_sliced will be checked twice but it doesn't really matter (it is more important to factorize the code)
			else
				v = std::move(this->transform->multiply(x,
							this->transposed2char(this->is_transposed, this->is_conjugate)));
				//NOTE: the function below depends on this one for the non-sliced case
			return v;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const FPP *x)
		{
			int x_size;
			// assuming that x size is valid, infer it from this size
			if(this->is_transposed)
				x_size = this->transform->getNbRow();
			else
				x_size = this->transform->getNbCol();
			if(this->is_sliced)
			{
				Vect<FPP, Cpu> v;
//				int v_size;
//				if(this->is_transposed)
//					v_size = this->transform->getNbCol();
//				else
//					v_size = this->transform->getNbRow();
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
				v.resize(this->getNbRow());
				sliceMultiply(this->slices, x, v.getData(), 1);
				return v;
#else
				this->eval_sliced_Transform();
				return this->multiply(x);
#endif
			}
			else if(this->is_fancy_indexed && this->is_all_dense()) // benchmarks have shown that indexMultiply worths it only if the Faust is dense
			{
				size_t id_lens[2] = {this->fancy_num_rows, this->fancy_num_cols};
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
				auto v = indexMultiply(this->fancy_indices, id_lens, x);
				return v;
#else
				this->eval_fancy_idx_Transform();
				return this->multiply(x);
#endif
			}
			else
			{
				this->eval_fancy_idx_Transform();
				Vect<FPP, Cpu> vx(x_size, x);
				return std::move(this->multiply(vx));
				// do not use the prototype below because it results in fact in a larger number of copies
			}
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::multiply(const FPP *x, FPP* y)
		{
			if(this->is_sliced)
			{
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
				sliceMultiply(this->slices, x, y, 1);
#else
				this->eval_sliced_Transform();
				return this->multiply(x, y);
#endif
			}
			else if(this->is_fancy_indexed && this->is_all_dense()) // benchmarks have shown that indexMultiply worths it only if the Faust is dense
			{
				size_t id_lens[2] = {this->fancy_num_rows, this->fancy_num_cols};
				auto y_vec = indexMultiply(this->fancy_indices, id_lens, x);
				memcpy(y, y_vec.getData(), sizeof(FPP)*y_vec.size()); // TODO: avoid this copy
			}
			else
			{
				this->eval_fancy_idx_Transform();
				int x_size;
				// assuming that x size is valid, infer it from this size
				if(this->is_transposed)
					x_size = this->transform->getNbRow();
				else
					x_size = this->transform->getNbCol();
				Vect<FPP, Cpu> vx(x_size, x);
				auto y_vec = std::move(this->multiply(vx));
				memcpy(y, y_vec.getData(), sizeof(FPP)*y_vec.size());
				// this alternative call is commented out because even if it avoids all copies made above, it is slower because the method above keeps only one Eigen vector to compute the whole product contrary to the following function that uses two vectors (output and operand vectors), which is slower
//				this->transform->multiply(x, y, this->transposed2char(this->is_transposed, this->is_conjugate));
			}
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> &A)
		{
			Faust::MatDense<FPP,Cpu> M;
#ifdef FAUST_TORCH
			if(this->mul_order_opt_mode >= TORCH_CPU_L2R && this->mul_order_opt_mode <= TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R && !tensor_data.size())
			{
				// init tensor data cache
				convMatGenListToTensorList(this->transform->data, tensor_data, at::kCPU, /* clone */ false, /* transpose */ ! this->is_transposed);
//				display_TensorList(tensor_data);
			}
#endif

			switch(this->mul_order_opt_mode)
			{
				case GREEDY_ALL_ENDS:
				case GREEDY_1ST_BEST:
				case GREEDY_ALL_BEST_CONVDENSE:
				case GREEDY_ALL_BEST_GENMAT:
					{
						std::vector<Faust::MatGeneric<FPP,Cpu>*> data(this->transform->data);
						data.resize(data.size()+1); // for the matrix multiplying Faust
						data[data.size()-1] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
						std::vector<char> transconj_flags = {'N'}; // 'N' for all Faust factors and the matrix A
						if(this->is_transposed)
						{
							transconj_flags = std::vector<char>(data.size(), this->isTransposed2char()); // TODO: use assign in C++ 20
							*(transconj_flags.end()-1) = 'N'; // the matrix multiplying is not transposed
							std::reverse(data.begin(), data.end()-1); // reversing the Faust factors (because they are transposed or transconjugate)
						}
						Faust::multiply_order_opt(this->mul_order_opt_mode, data, M, /*alpha */ FPP(1.0), /* beta */ FPP(0.0), transconj_flags);
					}
					break;
				case DYNPROG:
					{
						std::vector<Faust::MatGeneric<FPP,Cpu>*> data = this->transform->data;
						if(this->is_transposed)
							std::reverse(data.begin(), data.end());
						M = std::move(dynprog_multiply(data, this->isTransposed2char(), &A));
					}
					break;
				case CPP_PROD_PAR_REDUC:
					M = Faust::multiply_par(this->transform->data, A, this->isTransposed2char());
					break;
				case OMP_PROD_PAR_REDUC:
					M = Faust::multiply_omp(this->transform->data, A, this->isTransposed2char());
					break;
#ifdef FAUST_TORCH
				case TORCH_CPU_L2R:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_GREEDY:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false,  /*clone */ false,/* chain_opt */ true, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ true, !this->is_transposed);
					break;
#endif
				default:
					M = this->transform->multiply(A, this->isTransposed2char());
					break;
			}
			return std::move(M);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::multiply(const FPP* A, int A_ncols, FPP* C)
		{
			if(this->is_sliced && (A_ncols == 1 || this->size() > 1)) // benchmarks have shown that a single factor Faust is less efficient to multiply a marix (A_ncols > 1) with sliceMultiply than using eval_sliced_Transform and multiply
			{
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
				this->sliceMultiply(this->slices, A, C, A_ncols);
#else
				this->eval_sliced_Transform();
				return multiply(A, A_ncols, C);
#endif
			}
			else if(this->is_fancy_indexed && this->is_all_dense()) // benchmarks have shown that indexMultiply worths it only if the Faust is dense
			{
				size_t id_lens[2] = {this->fancy_num_rows, this->fancy_num_cols};
				auto mat = this->indexMultiply(this->fancy_indices, id_lens, A, A_ncols);
				memcpy(C, mat.getData(), sizeof(FPP)*mat.getNbRow()*mat.getNbCol());
			}
			else
			{
				this->eval_sliced_Transform();
				this->eval_fancy_idx_Transform();
				this->transform->multiply(A, A_ncols, C, this->isTransposed2char());
			}
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize(const bool transp /* deft to false */)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			//TODO: need a nsamples argument to feed optimize_time*
			Faust::TransformHelper<FPP,Cpu> *th = this->pruneout(/*nnz_tres=*/0), *th2;
			th2 = th->optimize_storage(false);
			delete th;
			th = th2;
			th->optimize_time(transp, true);
			return th;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_time(const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			// choose the quickest method for the Faust "toarray"
			auto t = this->optimize_time_full(transp, inplace, nsamples);
			return t;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_time_full(const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			return this->optimize_multiply([this](){this->get_product();}, transp, inplace, nsamples, "Faust-toarray");
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_time_prod(const MatGeneric<FPP, Cpu>* test_mat, const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			std::function<void(void)> benchmark_func;
			auto md = dynamic_cast<const MatDense<FPP,Cpu>*>(test_mat);
			auto ms = dynamic_cast<const MatSparse<FPP,Cpu>*>(test_mat);
			if(! md && ! ms)
				throw std::runtime_error("optimize_time_prod supports only MatDense or MatSparse benchmarking.");
			return this->optimize_multiply([this, ms, md]()
					{
					if(md) this->multiply(*md);
					else /* ms != nullptr */ this->multiply(*ms);
					}, transp, inplace, nsamples, "Faust-matrix product");
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_multiply(std::function<void()> f, const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples, const char* op_name)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			std::vector<string> meth_names = {"DEFAULT_L2R", "GREEDY_ALL_ENDS", "GREEDY_1ST_BEST", "GREEDY_ALL_BEST_CONVDENSE", "GREEDY", "DYNPROG", "CPP_PROD_PAR_REDUC", "OMP_PROD_PAR_REDUC", "TORCH_CPU_L2R", "TORCH_CPU_GREEDY","TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R" }; //TODO: it should be a function of faust_prod_opt module
			// GREEDY_ALL_BEST_GENMAT is printed out as GREEDY to follow the wrappers name
			TransformHelper<FPP,Cpu>* t_opt = nullptr;
			int NMETS = 11;
			std::chrono::duration<double> * times = new std::chrono::duration<double>[NMETS]; //use heap because of VS14 (error C3863)
//			MatDense<FPP,Cpu>* M = MatDense<FPP,Cpu>::randMat(transp?this->getNbRow():this->getNbCol(), 2048);
			int old_meth = this->get_mul_order_opt_mode();
			int nmuls = nsamples, opt_meth=0;
			std::vector<int> disabled_meths = {CPP_PROD_PAR_REDUC, OMP_PROD_PAR_REDUC}; // skip openmp/C++ threads methods because they are unfruitful when Eigen is multithreaded
			// disable experimental (mostly less efficient) greedy methods
			disabled_meths.push_back(GREEDY_ALL_ENDS);
			disabled_meths.push_back(GREEDY_1ST_BEST);
			disabled_meths.push_back(GREEDY_ALL_BEST_CONVDENSE);
			// even the GREEDY_ALL_BEST_GENMAT that is far better than previous ones
			disabled_meths.push_back(GREEDY_ALL_BEST_GENMAT);
#if DEBUG_OPT_MUL
			cout << "nsamples used to measure time: " << nmuls << endl;
#endif
#ifdef FAUST_TORCH
			if(this->mul_order_opt_mode >= TORCH_CPU_L2R && this->mul_order_opt_mode <= TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R && !tensor_data.size())
			{
				// init tensor data cache
				convMatGenListToTensorList(this->transform->data, tensor_data, at::kCPU, /* clone */ false, /* transpose */ ! this->is_transposed);
//				display_TensorList(tensor_data);
			}
#else
			disabled_meths.push_back(TORCH_CPU_L2R);
			disabled_meths.push_back(TORCH_CPU_GREEDY);
			disabled_meths.push_back(TORCH_CPU_DENSE_DYNPROG_SPARSE_L2R);
#endif

			for(int i=0; i < NMETS; i++)
			{
				if(std::find(std::begin(disabled_meths), std::end(disabled_meths), i) != std::end(disabled_meths))
				{
					times[i] = std::chrono::duration<double>(numeric_limits<double>::max());
					continue;
				}
				this->set_FM_mul_mode(i);
				auto start = std::chrono::system_clock::now();
				for(int j=0;j < nmuls; j++)
				{
					f();
				}
				auto end = std::chrono::system_clock::now();
				times[i] = end-start;
			}
			for(int i=0; i < NMETS-1; i++)
			{
				opt_meth = times[opt_meth]<times[i+1]?opt_meth:i+1;
			}
			if(inplace)
			{
				this->set_FM_mul_mode(opt_meth);
				t_opt = this;
			}
			else
			{
				t_opt = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0, false, false, true);
				cout << "best method measured in time on operation "<< op_name << " is: " << meth_names[opt_meth] << endl;
#if DEBUG_OPT_MUL
				cout << "all times: ";
				for(int i = 0; i < NMETS; i ++)
					cout << times[i].count() << " ";
				cout << endl;
#endif
				t_opt->set_FM_mul_mode(opt_meth);
				// leave the current Faust unchanged
				this->set_FM_mul_mode(old_meth);
			}
			delete [] times;
			t_opt->copy_transconj_state(*this);
			return t_opt;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::pruneout(const int nnz_tres, const int npasses, const bool only_forward)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			//TODO: refactor into src/algorithm/faust_pruneout.h/hpp in non-member function Faust::pruneout(TransformHelper<FPP,Cpu>& th, const int nnz_tres, const int npasses, const bool only_forward)
			int _npasses = 0;
			TransformHelper<FPP,Cpu> *pth = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0);
			MatGeneric<FPP,Cpu>* S_i, *S_j;
			MatDense<FPP,Cpu>* tmp_ds;
			MatSparse<FPP,Cpu>* tmp_sp;
			int nnz_i;
			bool factor_touched;
			while(_npasses < npasses || npasses == -1)
			{
				factor_touched = false;
				// forward pass
				for(int i = 0; (faust_unsigned_int)i < pth->size()-1; i++)
				{
					S_i = const_cast<Faust::MatGeneric<FPP,Cpu>*>(pth->get_gen_fact(i));
					S_j = const_cast<Faust::MatGeneric<FPP,Cpu>*>(pth->get_gen_fact(i+1));
					for(int offset = 0; offset<S_i->getNbCol(); offset++)
					{
						nnz_i = nnz_tres+1;
						{
							if((tmp_sp = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(S_i)))
							{
								//Matrix is read-only because it's RowMajor order
								Eigen::SparseMatrix<FPP,Eigen::ColMajor> sp_col = tmp_sp->mat.col(offset);
								nnz_i = sp_col.nonZeros();
								if(nnz_i <= nnz_tres)
								{

									//								cout << "nnz_i: " << nnz_i << " i: " << i<< endl;
									//								cout << "del col :" << offset << " fac:" << i << endl;
									tmp_sp->delete_col(offset);
									factor_touched = true;
								}
							}
							else
							{
								tmp_ds = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(S_i);
								nnz_i = tmp_ds->mat.col(offset).nonZeros();
								if(nnz_i <= nnz_tres)
								{
									//								cout << "nnz_i: " << nnz_i << " i: " << i<< endl;
									//								cout << "del col :" << offset << " fac:" << i << endl;
									tmp_ds->delete_col(offset);
									factor_touched = true;
								}
							}
							if((tmp_sp = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(S_j)))
							{
								if(nnz_i <= nnz_tres)
									tmp_sp->delete_row(offset);
							}
							else
							{
								tmp_ds = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(S_j);
								if(nnz_i <= nnz_tres)
									tmp_ds->delete_row(offset);
							}
						}
					}
				}
				// backward pass
				if(! only_forward)
					for(int i = pth->size()-1; i > 0; i--)
					{
						S_i = const_cast<Faust::MatGeneric<FPP,Cpu>*>(pth->get_gen_fact(i-1));
						S_j = const_cast<Faust::MatGeneric<FPP,Cpu>*>(pth->get_gen_fact(i));
						for(int offset = 0; offset<S_j->getNbRow(); offset++)
						{
							nnz_i = nnz_tres+1;
							{
								if(tmp_sp = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(S_j))
								{
									Eigen::SparseMatrix<FPP,Eigen::ColMajor> sp_row = tmp_sp->mat.row(offset);
									nnz_i = sp_row.nonZeros();
									if(nnz_i <= nnz_tres)
									{
										//								cout << "nnz_i: " << nnz_i << " i: " << i<< endl;
										//								cout << "del row :" << offset << " fac:" << i<< endl;
										tmp_sp->delete_row(offset);
										factor_touched = true;
									}
								}
								else
								{
									tmp_ds = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(S_j);
									nnz_i = tmp_ds->mat.row(offset).nonZeros();
									if(nnz_i <= nnz_tres)
									{
										//								cout << "nnz_i: " << nnz_i << " i: " << i<< endl;
										//								cout << "del row i:" << offset << " fac:" << i << endl;
										tmp_ds->delete_row(offset);
										factor_touched = true;
									}
								}
								if(tmp_sp = dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(S_i))
								{
									if(nnz_i <= nnz_tres)
										tmp_sp->delete_col(offset);
								}
								else
								{
									tmp_ds = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(S_i);
									if(nnz_i <= nnz_tres)
										tmp_ds->delete_col(offset);
								}

							}
						}

					}

				_npasses++;
				if(!factor_touched && npasses == -1) break;
			}
			// remove 0x0 factors
			for(int i = pth->size()-2; i > 0; i--)
			{
				auto fac_nrows = pth->get_fact_nb_rows(i);
				auto fac_ncols = pth->get_fact_nb_cols(i);
//				std::cout << "fac nrows, ncols:" << fac_nrows << " " << fac_ncols << std::endl;
				if(fac_nrows == 0 && fac_ncols == 0)
					pth->transform->erase(i);
			}
			pth->transform->update_total_nnz();
			pth->copy_transconj_state(*this);
			return pth;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::update_total_nnz()
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			this->transform->update_total_nnz();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::set_FM_mul_mode(const int mul_order_opt_mode, const bool silent /* = true */)
		{
			this->mul_order_opt_mode = mul_order_opt_mode;
			if(! silent)
			{
				std::cout << "changed mul. optimization mode to: " << this->mul_order_opt_mode;
				if(! this->mul_order_opt_mode)
					std::cout << " (opt. disabled, default mul.)";
				std::cout << std::endl;
			}

		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(const TransformHelper<FPP, Cpu>* th_right) const
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
			return new TransformHelper<FPP,Cpu>(this, th_right);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(FPP& scalar)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			const vector<MatGeneric<FPP,Cpu>*>& vec = this->transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
			//the point here is to minimize the number of copies (with direct access)
			// the constructor then will copy the factors from the vector
//			Transform<FPP,Cpu>* t = new Transform<FPP,Cpu>(vec, scalar, false, true); //optimizedCopy == false, cloning_fact == true
//			TransformHelper<FPP,Cpu>* th  = new TransformHelper<FPP,Cpu>(*t);
			TransformHelper<FPP,Cpu>* th = new TransformHelper<FPP,Cpu>(vec, scalar, false, false, true);
			th->copy_transconj_state(*this);
			th->copy_slice_state(*this);
			th->copy_fancy_idx_state(*this);
			return th;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::pop_back()
        {
            this->transform->pop_back();
        }

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::pop_front()
        {
            this->transform->pop_front();
        }

	template<typename FPP>
			void TransformHelper<FPP,Cpu>::push_back(const FPP* bdata, const int* brow_ptr, const int* bcol_inds, const int nrows, const int ncols, const int bnnz, const int bnrows, const int bncols, const bool optimizedCopy/*=false*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
			{
				this->eval_sliced_Transform();
				this->eval_fancy_idx_Transform();
				auto bsr_mat = new MatBSR<FPP, Cpu>(nrows, ncols, bnrows, bncols, bnnz, bdata, brow_ptr, bcol_inds);
				auto copying = optimizedCopy||transpose||conjugate;
				this->push_back(bsr_mat, optimizedCopy, copying, transpose, conjugate);
				if(copying) delete bsr_mat;
			}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const FPP* data, int nrows, int ncols, const bool optimizedCopy/*=famse*/, const bool transpose/*=famse*/, const bool conjugate/*=famse*/)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			auto copying = optimizedCopy||transpose||conjugate;
			auto dense_mat = new MatDense<FPP, Cpu>(data, nrows, ncols);
			this->push_back(dense_mat, optimizedCopy, copying, transpose, conjugate);
			if(copying) delete dense_mat;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy /* false by deft */, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			auto sparse_mat = new MatSparse<FPP,Cpu>(nnz, nrows, ncols, data, row_ptr, id_col);
			auto copying = optimizedCopy||transpose||conjugate;
			this->push_back(sparse_mat, optimizedCopy, copying, transpose, conjugate);
			if(copying) delete sparse_mat;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /* false by default */, const bool copying /* true to default */, const bool transpose/*=false*/, const bool conjugate/*=false*/, const bool verify_dims_agree/*=true*/)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			//warning: should not be called after initialization of factors (to respect the immutability property)
			//this function is here only for the python wrapper (TODO: see how to modify that wrapper in order to delete this function after or just use it internally -- not py/matfaust)
			this->transform->push_back(M, optimizedCopy, transpose, conjugate, copying, verify_dims_agree); // 2nd argument is for opt. (possibly converting dense <-> sparse)
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /* false by default */, const bool copying /* true to default */, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			this->push_back(M, optimizedCopy, copying, transpose, conjugate, /* verify_dims_agree */ true);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_first(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /* false by default */, const bool copying /* true to default */)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			//warning: should not be called after initialization of factors (to respect the immutability property)
			//this function is here only for python wrapper (TODO: see how to modify that wrapper in order to delete this function after or just use it internally -- not py/matfaust)
			this->transform->push_first(M, optimizedCopy, this->is_conjugate, copying); //2nd argument is for opt. (possibly converting dense <-> sparse)
		}

	template<typename FPP>
		template<typename Head, typename ... Tail>
		void TransformHelper<FPP,Cpu>::push_back_(Head& h, Tail&... t)
		{
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
//			for(auto f: h)
//				this->push_back(f, false, false);
			for(auto it=h.begin(); it < h.end(); it++)
			{
				auto f = *it;
				this->push_back(f, false, false);
			}
//			this->push_back(h, false, false);
			this->push_back_(t...);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back_()
		{
			// do nothing, here just for empty tail of above function
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNBytes() const
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
			faust_unsigned_int nbytes = 0;
			for(auto fac : this->transform->data)
			{
					nbytes += fac->getNBytes();
			}
			return nbytes;
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_total_nnz() const
		{
			return this->transform->get_total_nnz();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::resize(faust_unsigned_int sz)
		{
			return this->transform->resize(sz);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::display() const
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
			std::cout << this->to_string() << std::endl;
		}

	template<typename FPP>
		std::string TransformHelper<FPP,Cpu>::to_string() const
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
			return this->transform->to_string(this->is_transposed);
		}

	//private
	template<typename FPP>
	const MatGeneric<FPP,Cpu>* TransformHelper<FPP,Cpu>::get_gen_fact(const faust_unsigned_int id) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->data[this->is_transposed?size()-id-1:id];
	}

template<typename FPP>
	MatGeneric<FPP,Cpu>* TransformHelper<FPP,Cpu>::get_gen_fact_nonconst(const faust_unsigned_int id) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->data[this->is_transposed?size()-id-1:id];
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::update(const MatGeneric<FPP,Cpu>& M, const faust_unsigned_int fact_id)
	{
		MatGeneric<FPP,Cpu>& M_ = const_cast<MatGeneric<FPP,Cpu>&>(M);
		// I promise I won't change M
		MatSparse<FPP,Cpu> *sp_mat, *sp_fact;
		MatDense<FPP,Cpu> *ds_mat, *ds_fact;
		MatGeneric<FPP,Cpu>* fact = get_gen_fact_nonconst(fact_id);
		if(sp_mat = dynamic_cast<MatSparse<FPP,Cpu>*>(&M_))
		{
			if(! (sp_fact = dynamic_cast<MatSparse<FPP,Cpu>*>(fact)))
				throw std::runtime_error("A sparse factor can't be updated with a dense factor");
			*sp_fact = *sp_mat;
		}
		else if(ds_mat = dynamic_cast<MatDense<FPP,Cpu>*>(&M_))
		{
			if(! (ds_fact = dynamic_cast<MatDense<FPP,Cpu>*>(fact)))
				throw std::runtime_error("A dense factor can't be updated with a sparse factor");
			*ds_fact = *ds_mat;
		}
		else
		{
			throw std::runtime_error("Only MatSparse and MatDense are accepted by TransformHelper::update().");
		}
		fact->set_id(M.is_id());
		update_total_nnz();
	}

template<typename FPP>
	void TransformHelper<FPP, Cpu>::convertToSparse()
	{
		this->eval_sliced_Transform();
		this->eval_fancy_idx_Transform();
		const MatDense<FPP,Cpu> * mat_dense;
		const MatSparse<FPP,Cpu> * mat_sparse;
		const MatBSR<FPP,Cpu> * mat_bsr;
		for(int i=0;i<this->size();i++)
		{
			if(mat_dense = dynamic_cast<const MatDense<FPP,Cpu>*>(this->get_gen_fact(i)))
			{
				mat_sparse = new MatSparse<FPP,Cpu>(*mat_dense);
				this->replace(mat_sparse, i);
			}
			else if(mat_bsr = dynamic_cast<const MatBSR<FPP,Cpu>*>(this->get_gen_fact(i)))
			{
				mat_sparse =  new MatSparse<FPP, Cpu>(mat_bsr->to_sparse()); // TODO: avoid the copy
				this->replace(mat_sparse, i);
			}
		}
	}

template<typename FPP>
	void TransformHelper<FPP, Cpu>::convertToDense()
	{
		this->eval_sliced_Transform();
		this->eval_fancy_idx_Transform();
		const MatDense<FPP,Cpu> * mat_dense;
		const MatSparse<FPP,Cpu> * mat_sparse;
		const MatBSR<FPP,Cpu> * mat_bsr;
		for(int i=0;i<this->size();i++)
		{
			if(mat_sparse = dynamic_cast<const MatSparse<FPP,Cpu>*>(this->get_gen_fact(i)))
			{
				mat_dense = new MatDense<FPP,Cpu>(*mat_sparse);
				this->replace(mat_dense, i);
			}
			else if(mat_bsr = dynamic_cast<const MatBSR<FPP,Cpu>*>(this->get_gen_fact(i)))
			{
				mat_dense =  new MatDense<FPP, Cpu>(mat_bsr->to_sparse()); // TODO: avoid the copy
				this->replace(mat_dense, i);
			}
		}
	}

template<typename FPP>
	void TransformHelper<FPP, Cpu>::replace(const MatGeneric<FPP, Cpu>* M, const faust_unsigned_int fact_id)
	{
		this->transform->replace(M, fact_id);
	}

template<typename FPP>
	unsigned long long TransformHelper<FPP,Cpu>::get_fact_addr(const faust_unsigned_int id) const
	{
		return (unsigned long long) this->transform->data[this->is_transposed?size()-id-1:id];
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
			FPP* bdata,
			int* brow_ptr,
			int* bcol_inds) const
	{
		if(id < 0 || id > this->size())
			throw std::domain_error("get_fact(BSR): index out of bounds.");
		if(this->get_fact_type(id) != BSR)
			throw std::runtime_error("get_fact(BSR): matrix requested is not a MatBSR.");
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		MatBSR<FPP, Cpu> *bmat = dynamic_cast<MatBSR<FPP, Cpu>*>(this->transform->data[id]);
		bmat->copy_bdata(bdata);
		bmat->copy_browptr(brow_ptr);
		bmat->copy_bcolinds(bcol_inds);
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact_bsr_info(const faust_unsigned_int id,
			size_t& bdata_sz,
			size_t& browptr_sz,
			size_t& bcolinds_sz,
			size_t& bnnz,
			size_t& bnrows,
			size_t& bncols) const
	{
		if(id < 0 || id > this->size())
			throw std::domain_error("get_fact_bsr_info(BSR): index out of bounds.");
		if(this->get_fact_type(id) != BSR)
			throw std::runtime_error("get_fact_bsr_info(BSR): matrix requested is not a MatBSR.");
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		MatBSR<FPP, Cpu> *bmat = dynamic_cast<MatBSR<FPP, Cpu>*>(this->transform->data[id]);
		bmat->get_buf_sizes(bdata_sz, browptr_sz, bcolinds_sz);
		bnnz = bmat->getNBlocks();
		bnrows = bmat->getNbBlockRow();
		bncols = bmat->getNbBlockCol();
	}


template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
			const int** rowptr,
			const int** col_ids,
			const FPP** elts,
			faust_unsigned_int* nnz,
			faust_unsigned_int* num_rows,
			faust_unsigned_int* num_cols) const
	{
 		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols);

	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
			int* rowptr,
			int* col_ids,
			FPP* elts,
			faust_unsigned_int* nnz,
			faust_unsigned_int* num_rows,
			faust_unsigned_int* num_cols,
			const bool transpose /* = false*/) const
	{
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, this->is_transposed ^ transpose);
		if(this->is_conjugate)
			Faust::conjugate(elts, *nnz);
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
			const FPP** elts,
			faust_unsigned_int* num_rows,
			faust_unsigned_int* num_cols) const
	{
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, elts, num_rows, num_cols);
	}

template<typename FPP>
	MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
	{
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		MatDense<FPP,Cpu> dense_factor;
		MatGeneric<FPP,Cpu>* factor_generic;
		factor_generic = this->transform->get_fact(this->is_transposed?size()-id-1:id);

		switch (factor_generic->getType())
		{
			case MatType::Dense :
				{
					MatDense<FPP,Cpu>* factor_dense_ptr = dynamic_cast<MatDense<FPP,Cpu>* > (factor_generic);
					dense_factor = (*factor_dense_ptr); //copy
				}
				break;

			case MatType::Sparse :
				{
					MatSparse<FPP,Cpu>* factor_sparse_ptr = dynamic_cast<MatSparse<FPP,Cpu>* > (factor_generic);
					dense_factor = (*factor_sparse_ptr); //copy
				}
				break;

			default:
				handleError("Faust::TransformHelper", "get_fact : unknown type of the factor matrix");
		}
		delete factor_generic;
		if(this->is_transposed) dense_factor.transpose();
		if(this->is_conjugate) dense_factor.conjugate();
		return dense_factor;
	}

//TODO: delete this function when the specific bug in FaustCoreCpp (py. wrapper)
// will be fixed (normally the parent function def. should be found, but the wrapper compilation fails anyway)
template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int &id,
			FPP* elts,
			faust_unsigned_int* num_rows,
			faust_unsigned_int* num_cols,
			const bool transpose /* default to false */) const
	{
		if(id == 0 || id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, elts, num_rows, num_cols, this->is_transposed ^ transpose);
		if(this->is_conjugate)
			Faust::conjugate(elts,*num_cols*(*num_rows));
	}




template<typename FPP>
	MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product(const int mul_order_opt_mode/*=-1*/) // const
	{
		this->eval_sliced_Transform();
		this->eval_fancy_idx_Transform();
		int old_FM_mul_mod = -1;
		MatDense<FPP, Cpu> P;
		// keep current mul mode to restore it on the function end
		if(mul_order_opt_mode > -1 && mul_order_opt_mode != this->mul_order_opt_mode)
		{
			old_FM_mul_mod = this->mul_order_opt_mode;
			this->set_FM_mul_mode(mul_order_opt_mode);
		}
		if(this->mul_order_opt_mode && this->size() > 1 /* no need to do something special if this is a single factor Faust*/)
		{
			if(this->mul_order_opt_mode == DYNPROG)
			{
				std::vector<Faust::MatGeneric<FPP,Cpu>*> data = this->transform->data;
				if(this->is_transposed)
					std::reverse(data.begin(), data.end());
				P = std::move(dynprog_multiply(data, this->isTransposed2char()));
			}
			else
			{
				MatSparse<FPP,Cpu> Id(this->getNbCol(), this->getNbCol());
				Id.setEyes();
				P = this->multiply(Id);
			}
		}
		else
			P = this->transform->get_product(this->isTransposed2char(), this->is_conjugate);
		if(old_FM_mul_mod != - 1)
		{
			this->set_FM_mul_mode(old_FM_mul_mod);
		}
		return P;
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::get_product(Faust::MatDense<FPP,Cpu>& prod, const int mul_order_opt_mode/*=-1*/) //const
	{
		this->eval_sliced_Transform();
		this->eval_fancy_idx_Transform();
		if(mul_order_opt_mode != DEFAULT_L2R)
			prod = this->get_product(mul_order_opt_mode);
		else
			this->transform->get_product(prod, this->isTransposed2char(), this->is_conjugate);
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::save_mat_file(const char* filepath) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		this->transform->save_mat_file(filepath, this->is_transposed, this->is_conjugate);
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::read_from_mat_file(const char* filepath)
	{
		this->transform->read_from_mat_file(filepath);
	}

template<typename FPP>
	int TransformHelper<FPP,Cpu>::get_mat_file_type(const char* filepath)
	{
		return	Transform<FPP,Cpu>::get_mat_file_type(filepath);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
//			std::cout << "TransformHelper<FPP,Cpu>::spectralNorm" << std::endl;
		vector <MatGeneric<FPP, Cpu>*>& orig_facts = this->transform->data;
		int start_id, end_id;
		this->transform->get_nonortho_interior_prod_ids(start_id, end_id);
		//			cout << "start_id=" << start_id << "end_id=" << end_id << endl;
		if(start_id < 0)
			return 1.0;
		else if(start_id == 0)
			return this->transform->spectralNorm(nbr_iter_max, threshold, flag);
		//			cout << "optimized norm2" << endl;
		vector<MatGeneric<FPP,Cpu>*> non_ortho_start(orig_facts.begin()+start_id, orig_facts.end());
		TransformHelper<FPP, Cpu> t(non_ortho_start, 1.0, false, false);
		return t.transform->spectralNorm(nbr_iter_max, threshold, flag);
	}

template<typename FPP>
	FPP Faust::TransformHelper<FPP,Cpu>::power_iteration(const faust_unsigned_int nbr_iter_max, const Real<FPP>& threshold, int & flag) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->power_iteration(nbr_iter_max, threshold, flag);
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::transpose()
	{
		return new TransformHelper<FPP,Cpu>(this, true, false);
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::conjugate()
	{
		return new TransformHelper<FPP,Cpu>(this, false, true);
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::adjoint() const
	{
		return new TransformHelper<FPP,Cpu>(this, true, true);
	}


template<typename FPP>
	double TransformHelper<FPP,Cpu>::normL1(const bool full_array/*=true*/, const int batch_sz/*=1*/) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->normL1(this->is_transposed, full_array, batch_sz);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::normInf(const bool full_array/*=true*/, const int batch_sz/*=1*/) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->normInf(this->is_transposed, full_array, batch_sz);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::normFro(const bool full_array/*=true*/, const int batch_sz/*=1*/) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		vector <MatGeneric<FPP, Cpu>*>& orig_facts = this->transform->data;
		int start_id, end_id;
		this->transform->get_nonortho_interior_prod_ids(start_id, end_id);
		if(start_id < 0)
			return Faust::fabs(MatDense<FPP,Cpu>::eye(this->getNbCol(), this->getNbCol()).norm());
		else if(start_id == 0)
			return this->transform->normFro(full_array, batch_sz);
		vector<MatGeneric<FPP,Cpu>*> non_ortho_start(orig_facts.begin()+start_id, orig_facts.end());
		TransformHelper<FPP, Cpu> t(non_ortho_start, 1.0, false, false);
		return t.transform->normFro(full_array, batch_sz);
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		unsigned int ncols = this->getNbCol();
		unsigned int nrows = this->getNbRow();
		//2-norm parameters
		double precision =  FAUST_PRECISION;
		faust_unsigned_int nbr_iter_max=100;
		int flag;

		vector<FPP> norm_invs(ncols);
		FPP norm;
		vector<int> coords(ncols);
		TransformHelper<FPP,Cpu>* normalizedTh = nullptr;

		for(faust_unsigned_int i=0;i<ncols;i++)
		{

			//TODO: we shouldn't have to const_cast, slice() must be const
			const TransformHelper<FPP,Cpu> * col = const_cast<TransformHelper<FPP,Cpu> *>(this)->slice(0,nrows, i, i+1);
			//				cout << "TransformHelper normalize, meth=" << meth << endl;
			switch(meth)
			{
				case 1: //1-norm
					//						cout << "normalize with 1-norm" << endl;
					norm = col->normL1();
					break;
				case 2: //2-norm
					//						cout << "normalize with 2-norm" << endl;
					norm = col->spectralNorm(nbr_iter_max, precision, flag);
					break;
				case -2: // fro norm
					//						cout << "normalize with fro-norm" << endl;
					norm = col->normFro();
					break;
				case -1: //inf norm
					//						cout << "normalize with inf-norm" << endl;
					norm = col->normInf();
					break;
				default:
					handleError("Faust::TransformHelper::normalize()", "order for the norm to use is not valid");
			}
			if(norm != FPP(0))
				norm_invs[i] = (FPP)1./norm;
			else
				norm_invs[i] = 1;
			coords[i] = i;
			delete col;
		}
		MatSparse<FPP,Cpu>* norm_diag = new MatSparse<FPP,Cpu>(coords, coords,
				norm_invs, ncols, ncols);
		//			norm_diag.Display();
		vector<MatGeneric<FPP,Cpu>*> factors;
		for(int i =0; (faust_unsigned_int)i < size(); i++)
			//NOTE: about const cast: no problem, we know we won't write it
			//NOTE: don't use get_gen_fact() to avoid transpose auto-handling
			factors.push_back(const_cast<Faust::MatGeneric<FPP, Cpu>*>(this->transform->data[i]));
#ifdef NON_OPT_FAUST_NORMALIZATION
		if(this->is_transposed)
			factors.insert(factors.begin(), norm_diag);
		else
			factors.push_back(norm_diag);
#else
		MatGeneric<FPP,Cpu>* scaled_f0;
		MatSparse<FPP, Cpu>* fs;
		MatDense<FPP,Cpu>* fd;

		if(this->is_transposed)
		{
			/** this approach is abandoned because casting approach (as below) is quicker than transposing */
			// the faust is transposed
			// that's why we compute (Faust[0]^T*norm_diag)^T
			// (the Faust structure will transpose this factor again when necessary because of its this->is_transposed flag)
			//				factors[0]->transpose();
			//				factors[0]->multiplyRight(*norm_diag);
			//				factors[0]->transpose();
			/****************************************/
			fs = dynamic_cast<MatSparse<FPP,Cpu>*>(factors[0]);
			if(!fs)
			{
				fd = dynamic_cast<MatDense<FPP,Cpu>*>(factors[0]);
				scaled_f0 = new MatDense<FPP,Cpu>(*fd);
				fd = (MatDense<FPP,Cpu>*)scaled_f0; // needed befause there is no prototype of multiply(MatGeneric,...)
				norm_diag->multiply(*fd, 'N');
			}
			else
			{
				scaled_f0 = new MatSparse<FPP,Cpu>(*fs);
				fs = (MatSparse<FPP,Cpu>*) scaled_f0;
				norm_diag->multiply(*fs, 'N');
			}
			factors[0] = scaled_f0;
		}
		else
		{
			//factors[size()-1]->multiplyRight(*norm_diag);
			if((fs = dynamic_cast<MatSparse<FPP,Cpu>*>(factors[size()-1])))
				scaled_f0 = new MatSparse<FPP,Cpu>(*fs);
			else
			{
				fd = dynamic_cast<MatDense<FPP,Cpu>*>(factors[size()-1]);
				scaled_f0 = new MatDense<FPP,Cpu>(*fd);

			}
			scaled_f0->multiplyRight(*norm_diag);
			factors[size()-1] = scaled_f0;
		}

		delete norm_diag;
#endif
		normalizedTh = new TransformHelper<FPP,Cpu>(factors, FPP(1.0), false, false);
		normalizedTh->is_transposed = this->is_transposed;
		//			normalizedTh->display();
		return normalizedTh;
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>::~TransformHelper()
	{
		// transform is deleted auto. when no TransformHelper uses it (no more weak refs left)
#ifdef FAUST_VERBOSE
		cout << "Destroying Faust::TransformHelper object." << endl;
#endif
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density /* 1.f */, bool per_row /* true */)
	{
		return TransformHelper<FPP,Cpu>::randFaust(-1, -1, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randBSRFaust(unsigned int faust_nrows, unsigned int faust_ncols, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int bnrows, unsigned int bncols, float density/*=.1f*/)
	{
		if(faust_nrows != faust_ncols)
			throw std::runtime_error("randBSRFaust: currently only random square BSR Faust can be generated.");
		if(faust_nrows%bnrows || faust_ncols%bncols)
			throw std::runtime_error("randBSRFaust: the size of blocks must evenly divide the size of Faust matrices");
		unsigned int bnnz = (unsigned int) std::round(faust_nrows*faust_ncols/bnrows/bncols*density);
		if(bnnz == 0)
			throw std::runtime_error("randBSRFaust: the nonzero blocks are too large for this Faust/matrix size.");
		// pick randomly the number of factors into {min_num_factors, ..., max_num_factors}
		std::uniform_int_distribution<int> num_fac_distr(min_num_factors, max_num_factors);
		int num_factors = num_fac_distr(generator);
		// create factors
		std::vector<MatGeneric<FPP,Cpu>*> factors(num_factors);
		for(int i=0;i<num_factors;i++)
			factors[i] = MatBSR<FPP, Cpu>::randMat(faust_nrows, faust_ncols, bnrows, bncols, bnnz);
		// create the Faust
		TransformHelper<FPP,Cpu>* randFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
		return randFaust;
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density /* 1.f */, bool per_row /* true */)
	{
		unsigned int tmp;
		if(!TransformHelper<FPP,Cpu>::seed_init) {
			std::srand(std::time(NULL)); //seed init needed for MatDense rand generation
			TransformHelper<FPP,Cpu>::seed_init = true;
		}
		// pick randomly the number of factors into {min_num_factors, ..., max_num_factors}
		std::uniform_int_distribution<int> num_fac_distr(min_num_factors, max_num_factors);
		if(min_dim_size > max_dim_size)
		{
			tmp = min_dim_size;
			min_dim_size = max_dim_size;
			max_dim_size = tmp;
		}
		std::uniform_int_distribution<int> dim_distr(min_dim_size, max_dim_size);
		std::uniform_int_distribution<int> bin_distr(0,1);
		unsigned int num_factors = num_fac_distr(generator);
		// create factors randomly respecting the RandFaustType asked and dims interval
		std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) num_factors);
		unsigned int num_rows, num_cols;
		// if faust_nrows < 0 then the number of rows is random between min_dim_size and max_dim_size
		if(faust_nrows < 0)
			num_cols = dim_distr(generator);
		else // fixed number of rows
			num_cols = faust_nrows;
		float fact_density;
		for(unsigned int i=0;i<num_factors;i++) {
			num_rows = num_cols;
			// last faust factor number of columns is the faust number of columns
			if(i < num_factors-1 || faust_ncols < 0)
				// it is random for intermediary factors and also for last factor if faust_ncols < 0
				num_cols = dim_distr(generator);
			else
				// the faust has a fixed number of columns
				num_cols = faust_ncols;
#ifdef FAUST_VERBOSE
			cout << "TransformHelper<FPP,Cpu>::randFaust() per_row: " <<  per_row << endl;
#endif
			if(density == -1.)
				fact_density = per_row?5./num_cols:5./num_rows;
			else
				fact_density = density;
			switch(t){
				case DENSE:
					factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
					break;
				case SPARSE:
					factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
					break;
				case MIXED:
					if(bin_distr(generator))
						factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
					else
						factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, fact_density, per_row);
					break;
				default:
					handleError("Faust::TransformHelper", "randFaust(): Unknown RandFaustType");
					break;

			}
			if(factors[i] == NULL) return NULL;
		}
		TransformHelper<FPP,Cpu>* randFaust = new TransformHelper<FPP, Cpu>(factors,1.0,false,false);
		return randFaust;
	}


template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::hadamardFaust(unsigned int n, const bool norma)
	{
		TransformHelper<FPP,Cpu>* hadamardFaust = nullptr;
		vector<MatGeneric<FPP,Cpu>*> factors;
		bool cloning_fact = false; // big opt. allowed only because of the RefManager used in Transform class
		//this opt. avoids to duplicate the same factor
		try {
			wht_factors(n, factors, cloning_fact, norma);
			hadamardFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
		}
		catch(std::bad_alloc)
		{
			// nothing to do, out of memory return nullptr
		}
		return hadamardFaust;
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::fourierFaust(unsigned int n, const bool norma)
	{

		vector<MatGeneric<FPP,Cpu>*> factors(n+1);
		TransformHelper<FPP,Cpu>* fourierFaust = nullptr;
		try
		{
			fft_factors(n, factors);
			FPP alpha = norma?FPP(1/sqrt((double)(1 << n))):FPP(1.0);
			fourierFaust = new TransformHelper<FPP, Cpu>(factors, alpha, false, false, /* internal call */ true);
		}
		catch(std::bad_alloc e)
		{
			//nothing to do, out of memory, return nullptr
		}
		return fourierFaust;
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::eyeFaust(unsigned int n, unsigned int m)
	{
		vector<MatGeneric<FPP,Cpu>*> factors(1);
		TransformHelper<FPP,Cpu>* eyeFaust = nullptr;
		try
		{
			MatSparse<FPP,Cpu>* eye = MatSparse<FPP,Cpu>::eye(n, m);
			factors[0] = eye;
			eyeFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
		}
		catch(std::bad_alloc e)
		{
			//nothing to do, out of memory, return nullptr
		}
		return eyeFaust;
	}

template<typename FPP>
	Faust::transf_iterator<FPP> Faust::TransformHelper<FPP, Cpu>::begin() const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->begin();
	}

template<typename FPP>
	Faust::transf_iterator<FPP> Faust::TransformHelper<FPP, Cpu>::end() const
	{
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
		const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		return this->transform->end();
	}

template<typename FPP>
	void Faust::TransformHelper<FPP,Cpu>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode/*=DEFAULT_L2R*/)
	{
		if(start_id == 0 || end_id == this->size()-1)
		{
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, Cpu>*>(this)->eval_fancy_idx_Transform();
		}
		// else no end factors is concerned by the packing, keep the "virtual" slice as is
		if(start_id < 0 || start_id >= size())
			throw out_of_range("start_id is out of range.");
		if(end_id < 0 || end_id >= size())
			throw out_of_range("end_id is out of range.");
		Faust::MatDense<FPP,Cpu> * packed_fac = nullptr;
		if(end_id == start_id)
		{
			//nothing to do except converting to MatDense if start_id
			//factor is a MatSparse
			packed_fac = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(*(begin()+start_id));
			//TODO: and MatBSR, MatDiag ?
			if(packed_fac == nullptr)
			{// factor start_id is not at MatDense, convert it
				packed_fac = new MatDense<FPP,Cpu>(*dynamic_cast<Faust::MatSparse<FPP,Cpu>*>(*(begin()+start_id)));
			}
			else
			{
				return; //no change
			}
		}
		else
		{
			// we have to multiply factors from start_id to end_id into one matrix
			// simple way to do, 1) create a overhead-free TransformHelper with these factors
			// 2) call get_product() to override the start_id factors with the result on the end
			// 3) erase factors from start_id to end_id and insert packed factor too to replace them (that's Transform object responsibility).
			// 1)
			std::vector<Faust::MatGeneric<FPP,Cpu>*> topack_factors(begin()+start_id, begin()+end_id+1);
			Faust::TransformHelper<FPP,Cpu> t(topack_factors, 1.0, false, false, false);
			t.set_FM_mul_mode(mul_order_opt_mode, /*silent*/ false);
			// 2)
			packed_fac = new MatDense<FPP,Cpu>(t.get_product());
		}
		// 3)
		faust_unsigned_int i = end_id;
		while(i>=start_id)
		{
			this->transform->erase(i);
			if(i == 0) break;
			i--;
		}
		this->transform->insert(start_id, packed_fac);
	}

template <typename FPP>
	void TransformHelper<FPP,Cpu>::pack_factors(const int mul_order_opt_mode/*=DEFAULT_L2R*/)
{
	TransformHelperGen<FPP,Cpu>::pack_factors(mul_order_opt_mode/*=DEFAULT_L2R*/);
}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::swap_rows(const faust_unsigned_int id1,
			const faust_unsigned_int id2,
			const bool permutation/*=false*/,
			const bool inplace/*=false*/,
			const bool check_transpose/*=true*/)
{
	this->eval_sliced_Transform();
	this->eval_fancy_idx_Transform();
	if(check_transpose && this->is_transposed)
		return swap_cols(id1, id2, permutation, inplace, false);
	TransformHelper<FPP,Cpu>* t = nullptr;
	if(inplace)
		t = this;
	else
	{
		t = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0,
				/* optimizedCopy*/ false, /* cloning_fact */ true, /* internal_call*/ true);
		t->copy_transconj_state(*this);
		t->copy_slice_state(*this);
		t->copy_fancy_idx_state(*this);
		t->copy_mul_mode_state(*this);
	}
	// don't use get_gen_fact_nonconst in case the the TransformHelper is transposed
	auto last_fac = t->transform->data[0];
	if(permutation)
	{
		auto P = MatSparse<FPP,Cpu>::swap_matrix(last_fac->getNbRow(), id1, id2);
		t->push_first(P, /*optimizedCopy*/ false, /* copying*/ false);
	}
	else
	{
		// swap two columns id1, id2 of the last factors
		auto dlast_fac = dynamic_cast<MatDense<FPP,Cpu>*>(last_fac);
		if(dlast_fac)
		{
			dlast_fac->swap_rows(id1, id2);
		}
		else
		{
			auto slast_fac = dynamic_cast<MatSparse<FPP,Cpu>*>(last_fac);
			slast_fac->swap_rows(id1, id2);
		}
	}
	return t;
}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::swap_cols(const faust_unsigned_int id1,
			const faust_unsigned_int id2,
			const bool permutation/*=false*/,
			const bool inplace/*=false*/,
			const bool check_transpose/*=true*/)
{
	this->eval_sliced_Transform();
	this->eval_fancy_idx_Transform();
	if(check_transpose && this->is_transposed)
		return swap_rows(id1, id2, permutation, inplace, false);
	TransformHelper<FPP,Cpu>* t = nullptr;
	if(inplace)
		t = this;
	else
	{
		t = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0,
				/* optimizedCopy*/ false, /* cloning_fact */ true, /* internal_call*/ true);
		t->copy_transconj_state(*this);
		t->copy_slice_state(*this);
		t->copy_fancy_idx_state(*this);
		t->copy_mul_mode_state(*this);
	}
	// don't use get_gen_fact_nonconst in case the the TransformHelper is transposed
	auto last_fac = t->transform->data[size()-1];
	if(permutation)
	{
		auto P = MatSparse<FPP,Cpu>::swap_matrix(last_fac->getNbCol(), id1, id2);
		t->push_back(P, /*optimizedCopy*/ false, /* copying*/ false, /* transpose */ false, /* conjugate */ false);
	}
	else
	{
		// swap two columns id1, id2 of the last factors
		auto dlast_fac = dynamic_cast<MatDense<FPP,Cpu>*>(last_fac);
		if(dlast_fac)
		{
			dlast_fac->swap_cols(id1, id2);
		}
		else
		{
			auto slast_fac = dynamic_cast<MatSparse<FPP,Cpu>*>(last_fac);
			slast_fac->swap_cols(id1,id2);
		}
	}
	return t;
}


	template<typename FPP>
FPP TransformHelper<FPP,Cpu>::get_item(faust_unsigned_int i, faust_unsigned_int j)
{
	MatDense<FPP, Cpu> M;
	faust_unsigned_int out_id;
	TransformHelperGen<FPP,Cpu>::get_item(i, j, M, out_id);
	return M.getData()[out_id];
}
	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));

}

template<typename FPP>
	template<typename FPP2>
Faust::TransformHelper<Real<FPP2>, Cpu>* Faust::TransformHelper<FPP, Cpu>::real()
{
	std::vector<MatGeneric<Real<FPP2>,Cpu>*> real_data;
	MatSparse<FPP, Cpu> *curfac_sp;
	MatDense<FPP, Cpu> *curfac_ds;
	for(auto curfac: this->transform->data)
	{
		if(curfac_ds = dynamic_cast<MatDense<FPP, Cpu>*>(curfac))
		{
			auto real_fac = new MatDense<Real<FPP2>,Cpu>(curfac->getNbRow(), curfac->getNbCol());
			*real_fac = curfac_ds->template to_real<Real<FPP2>>();
			real_data.push_back(real_fac);
		}
		else if(curfac_sp = dynamic_cast<MatSparse<FPP, Cpu>*>(curfac))
		{
			auto real_fac = new MatSparse<Real<FPP2>,Cpu>(curfac->getNbRow(), curfac->getNbCol());
			*real_fac = curfac_sp->template to_real<Real<FPP2>>();
			real_data.push_back(real_fac);
		}
		else
		{
			throw std::runtime_error("real() failed because a factor is neither a MatDense or a MatSparse");
		}
	}
	auto th = new TransformHelper<Real<FPP2>, Cpu>(real_data, 1.0, false, false, true);
	return th;
}

template<typename FPP>
bool Faust::TransformHelper<FPP,Cpu>::is_zero() const
{
	return this->transform.is_zero;
}

template<typename FPP>
Faust::Vect<FPP, Cpu> Faust::TransformHelper<FPP,Cpu>::indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* x) const
{
	auto out_nrows = ids[0]==nullptr?this->getNbRow():ids_len[0];
	Vect<FPP, Cpu> out_vec(out_nrows);
	indexMultiply(ids, ids_len, x, out_vec.getData());
	return out_vec;
}

template<typename FPP>
FPP* Faust::TransformHelper<FPP,Cpu>::indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* x, FPP *out) const
{
	return indexMultiply(ids, ids_len, x, 1, out);
}

template<typename FPP>
Faust::MatDense<FPP, Cpu> Faust::TransformHelper<FPP,Cpu>::indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* X, int ncols) const
{
	MatDense<FPP, Cpu> out_mat;
	auto out_nrows = ids[0]==nullptr?this->getNbRow():ids_len[0];
	out_mat.resize(out_nrows, ncols);
	indexMultiply(ids, ids_len, X, ncols, out_mat.getData());
	return out_mat;
}

template<typename FPP>
FPP* Faust::TransformHelper<FPP,Cpu>::indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* X, int ncols, FPP* out) const
{
	std::cout << "indexMultiply" << std::endl;
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
	//TODO: refactor with sliceMultiply if it applies
	using Mat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>;
	using SpMat = Eigen::SparseMatrix<FPP, Eigen::ColMajor>;
	using Mat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>;
	using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
	MatMap X_map(const_cast<FPP*>(X), ids[1]==nullptr?this->getNbCol():ids_len[1], ncols); // the const_cast is harmless, no modif. is made
																					// getNbCol is transpose aware
	MatDense<FPP, Cpu>* ds_fac;
	MatSparse<FPP, Cpu>* sp_fac;
//	MatBSR<FPP, Cpu>* bsr_fac; // TODO
	Mat M;
	MatDense<FPP, Cpu> mat;
	auto gen_first_fac = *(this->transform->begin());
	auto gen_last_fac = *(this->transform->end()-1);
	MatGeneric<FPP, Cpu>* left_indexed_factor = nullptr;
	auto out_nrows = ids[0]==nullptr?this->getNbRow():ids_len[0];
	if(out == nullptr)
		//TODO: document that out must be deleted by the callee
		out = new FPP[out_nrows*ncols];
	MatMap out_map(out, out_nrows, ncols);
	if(this->is_transposed)
	{
		// 1. multiply the first factor indexed according to ids by X
		// 1.1 get the first factor
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(gen_first_fac))
		{
			// 1.2 index ds_fac and multiply by X
			// eigenIndexMul ignores ids_len[0]if ids[0] is nullptr, so no worry
			if(size() == 1)
				ds_fac->eigenIndexMul(ids[0], ids[1], ids_len[0], ids_len[1], X_map, out_map, this->is_transposed, this->is_conjugate);
			else
				ds_fac->eigenIndexMul(nullptr, ids[1], ids_len[0], ids_len[1], X_map, M, this->is_transposed, this->is_conjugate);
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(gen_first_fac))
		{
			// 1.2 multiply the ids columns of the factor by X
			if(size() == 1)
				sp_fac->eigenIndexMul(ids[0], ids[1], ids_len[0], ids_len[1], X_map, out_map, this->is_transposed, this->is_conjugate);
			else
				sp_fac->eigenIndexMul(nullptr, ids[1], ids_len[0], ids_len[1], X_map, M, this->is_transposed, this->is_conjugate);
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
		// 2. then create a subFaust with all factors except the first one (copy-free)
		if(size() == 1)
			return out;
		auto end_fac_it = this->end();
		if(ids[0] != nullptr)
		{
			end_fac_it -= 1;
			left_indexed_factor = gen_last_fac;
		}
		std::vector<MatGeneric<FPP,Cpu> *> other_facs(this->begin()+1, end_fac_it);
		TransformHelper<FPP, Cpu> sub_th(other_facs, (FPP) 1.0, /* opt_copy */false, /* cloning */ false, /* internal call */ true);
		// 3. finally multiply this new Faust by M and return the result
		auto sub_th_t = sub_th.transpose();
		if(ids[0] == nullptr)
		{
			sub_th_t->multiply(M.data(), M.cols(), out);
		}
		else
		{
			mat.resize(sub_th_t->getNbRow(), M.cols());
			sub_th_t->multiply(M.data(), M.cols(), mat.getData());
		}
		delete sub_th_t;
	}
	else
	{
		// 1. multiply the indexed last factor according to ids by X
		// 1.1 get the last factor
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(gen_last_fac))
		{
			// 1.2 multiply the slice of the factor by X
			// eigenIndexMul ignores ids_len[0]if ids[0] is nullptr, so no worry
			if(size() == 1)
				ds_fac->eigenIndexMul(ids[0], ids[1], ids_len[0], ids_len[1], X_map, out_map, this->is_transposed, this->is_conjugate);
			else
				ds_fac->eigenIndexMul(nullptr, ids[1], ids_len[0], ids_len[1], X_map, M, this->is_transposed, this->is_conjugate);
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(gen_last_fac))
		{
			// 1.2 multiply the indexed factor by X
			if(size() == 1)
				sp_fac->eigenIndexMul(ids[0], ids[1], ids_len[0], ids_len[1], X_map, out_map, this->is_transposed, this->is_conjugate);
			else
				sp_fac->eigenIndexMul(nullptr, ids[1], ids_len[0], ids_len[1], X_map, M, this->is_transposed, this->is_conjugate);
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
		// 2. then create a subFaust with all factors except the last one (copy-free) and the first one if indexing on the row too
		if(size() == 1)
			return out;
		auto start_fac_it = this->begin();
		if(ids[0] != nullptr)
		{
			start_fac_it += 1;
			left_indexed_factor = gen_first_fac;
		}
		std::vector<MatGeneric<FPP,Cpu> *> other_facs(start_fac_it, this->end()-1);
		TransformHelper<FPP, Cpu> sub_th(other_facs, (FPP) 1.0, /* opt_copy */false, /* cloning */ false, /* internal call */ true);
		// 3. finally multiply this new Faust by M and return the result
		if(ids[0] == nullptr)
		{
			sub_th.multiply(M.data(), M.cols(), out);
		}
		else
		{
			mat.resize(sub_th.getNbRow(), M.cols());
			sub_th.multiply(M.data(), M.cols(), mat.getData());
		}
	}
	if(ids[0] != nullptr)
	{
		MatMap X_map(const_cast<FPP*>(mat.getData()), mat.getNbRow(), mat.getNbCol());
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(left_indexed_factor))
		{
			// 1.2 multiply the slice of the factor by X
			ds_fac->eigenIndexMul(ids[0], nullptr, ids_len[0], 0, X_map, out_map, this->is_transposed, this->is_conjugate);
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(left_indexed_factor))
		{
			// 1.2 multiply the indexed factor by X
			sp_fac->eigenIndexMul(ids[0], nullptr, ids_len[0], 0, X_map, out_map, this->is_transposed, this->is_conjugate);
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
	}
#else
	throw std::runtime_error("TransformHelper::indexMultiply is not supported with eigen version < 3.9");
#endif

	return out;
}

template<typename FPP>
Faust::MatDense<FPP, Cpu> Faust::TransformHelper<FPP,Cpu>::sliceMultiply(const Slice s[2], const FPP* X, int X_ncols/*=1*/) const
{
	Faust::MatDense<FPP, Cpu> out(s[0].size(), X_ncols);
	sliceMultiply(s, X, out.getData(), X_ncols);
	return out;
}

template<typename FPP>
FPP* Faust::TransformHelper<FPP,Cpu>::sliceMultiply(const Slice s[2], const FPP* X, FPP* out/*=nullptr*/, int X_ncols/*=1*/) const
{
#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)
	// NOTE: Take care if you edit this method, it must avoid any method that calls eval_sliced_Transform or it would became useless or bugged
	//TODO: refactor this function (too long)
	//TODO: refactor with MatDense/MatSparse/eigenSliceMul, similar to eigenIndexMul
	using Mat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>;
	using MatMap = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
	MatMap X_map(const_cast<FPP*>(X), this->getNbCol(), X_ncols); // the const_cast is harmless, no modif. is made
														 // getNbCol is transpose aware
	MatDense<FPP, Cpu>* ds_fac;
	MatSparse<FPP, Cpu>* sp_fac;
//	MatBSR<FPP, Cpu>* bsr_fac; // TODO
	Mat M;
	auto s0 = s[0]; // row slice
	auto s1 = s[1]; // column slice
	int rf_row_offset, rf_nrows;
	auto gen_first_fac = *(this->transform->begin()); // don't use this->begin because it can launch eval_sliced_Transform
	auto gen_last_fac = *(this->transform->end()-1);
	MatGeneric<FPP, Cpu>* left_sliced_factor = nullptr;
	std::vector<MatGeneric<FPP,Cpu> *> other_facs;
	if(out == nullptr)
		out = new FPP[s[0].size()*X_ncols];
	MatMap out_map(out, s[0].size(), X_ncols);
	FPP* tmp_ptr;
	int tmp_nrows;

	if(size() == 1)
	{
		// might be row and column slicing on the right factor
		rf_row_offset = s0.start_id;
		rf_nrows = s0.size();
	}
	else
	{
		// only column slicing on the right factor
		rf_row_offset = 0;
		if(this->is_transposed)
			rf_nrows = gen_first_fac->getNbCol();
		else
			rf_nrows = gen_last_fac->getNbRow();
	}
	if(this->is_transposed)
	{
		// 1. multiply the first factor sliced according to s1 to the corresponding portion of x in the new vector x_
		// 1.1 get the first factor
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(gen_first_fac))
		{
			// 1.2 multiply the slice of the factor by X
			if(size() == 1)
				if(this->is_conjugate)
					out_map = ds_fac->mat.adjoint().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					out_map = ds_fac->mat.transpose().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
			else
				if(this->is_conjugate)
					M = ds_fac->mat.adjoint().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					M = ds_fac->mat.transpose().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(gen_first_fac))
		{
			// 1.2 multiply the slice of the factor by X
			if(size() == 1)
				if(this->is_conjugate)
					out_map = sp_fac->mat.adjoint().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					out_map = sp_fac->mat.transpose().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
			else
				if(this->is_conjugate)
					M = sp_fac->mat.adjoint().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					M = sp_fac->mat.transpose().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
		if(size() == 1)
			return out;
		// 2. then create a subFaust with all factors except the first one (copy-free)
		// but a slicing on rows is also possible on the left factor
		if(s0.start_id != 0 || s0.end_id != gen_last_fac->getNbCol())
		{
			other_facs.assign(this->transform->begin()+1, this->transform->end()-1);
			left_sliced_factor = gen_first_fac;
			tmp_nrows = (*(this->transform->end()-2))->getNbCol();
			// the following sub_th product is not ensured to fit out buffer, so use another one
			tmp_ptr = new FPP[tmp_nrows*X_ncols];
		}
		else
		{
			other_facs.assign(this->transform->begin()+1, this->transform->end());
			tmp_ptr = out;
		}
		TransformHelper<FPP, Cpu> sub_th(other_facs, (FPP) 1.0, /* opt_copy */false, /* cloning */ false, /* internal call */ true);
		// 3. finally multiply this new Faust by M and return the result
		TransformHelper<FPP, Cpu> *sub_th_t;
		if(this->is_conjugate)
			sub_th_t = sub_th.adjoint();
		else
			sub_th_t = sub_th.transpose();
		sub_th_t->multiply(M.data(), M.cols(), tmp_ptr);
		delete sub_th_t;
	}
	else
	{
		// 1. multiply the last factor sliced according to s1 to the corresponding portion of x in the new vector x_
		// 1.1 get the last factor
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(gen_last_fac))
		{
			// 1.2 multiply the slice of the factor by X
			if(size() == 1)
				if(this->is_conjugate)
					out_map = ds_fac->mat.conjugate().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					out_map = ds_fac->mat.block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
			else
				if(this->is_conjugate)
					M = ds_fac->mat.conjugate().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					M = ds_fac->mat.block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(gen_last_fac))
		{
			// 1.2 multiply the slice of the factor by X
			if(size() == 1)
				if(this->is_conjugate)
					out_map = sp_fac->mat.conjugate().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					out_map = sp_fac->mat.block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
			else
				if(this->is_conjugate)
					M = sp_fac->mat.conjugate().block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
				else
					M = sp_fac->mat.block(rf_row_offset, s1.start_id, rf_nrows, s1.size()) * X_map;
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
		// 2. then create a subFaust with all factors except the last one (copy-free)
		if(size() == 1)
			return out;
		// but a slicing on rows is also possible on the left factor
		if(s0.start_id != 0 || s0.end_id != gen_first_fac->getNbRow())
		{
			other_facs.assign(this->transform->begin()+1, this->transform->end()-1);
			left_sliced_factor = gen_first_fac;
			tmp_nrows = (*(this->transform->begin()+1))->getNbRow();
			// the following sub_th product is not ensured to fit out buffer, so use another one
			tmp_ptr = new FPP[tmp_nrows*X_ncols];
		}
		else
		{
			other_facs.assign(this->transform->begin(), this->transform->end()-1);
			tmp_ptr = out;
		}
		// 3. finally multiply this new Faust by M and return the result
		TransformHelper<FPP, Cpu> sub_th(other_facs, (FPP) 1.0, /* opt_copy */false, /* cloning */ false, /* internal call */ true);
		if(this->is_conjugate)
		{
			auto sub_th_c = sub_th.conjugate();
			sub_th_c->multiply(M.data(), M.cols(), tmp_ptr);
			delete sub_th_c;
		}
		else
			sub_th.multiply(M.data(), M.cols(), tmp_ptr);
	}
	if(left_sliced_factor != nullptr)
	{
		// slice left factor too
		MatMap tmp_map(tmp_ptr, tmp_nrows, X_ncols);
		if(ds_fac = dynamic_cast<MatDense<FPP, Cpu>*>(left_sliced_factor))
		{
			// multiply the slice of the left factor by the right product
			Mat tmp;
			if(this->is_transposed)
			{
				if(this->is_conjugate)
					tmp = ds_fac->mat.adjoint();
				else
					tmp = ds_fac->mat.transpose();
			}
			else
			{
				if(this->is_conjugate)
					tmp = ds_fac->mat.conjugate();
				else
					tmp = ds_fac->mat;
			}
			out_map = tmp.block(s0.start_id, 0, s0.size(), tmp.cols()) * tmp_map;
		}
		else if(sp_fac = dynamic_cast<MatSparse<FPP, Cpu>*>(left_sliced_factor))
		{
			// multiply the slice of the left factor by the right product
			Eigen::SparseMatrix<FPP, Eigen::RowMajor> tmp;
			if(this->is_transposed)
				if(this->is_conjugate)
					tmp = sp_fac->mat.adjoint();
				else
					tmp = sp_fac->mat.transpose();
			else
				if(this->is_conjugate)
					tmp = sp_fac->mat.conjugate();
				else
					tmp = sp_fac->mat;
			out_map = tmp.block(s0.start_id, 0, s0.size(), tmp.cols()) * tmp_map;
		}
		else
		{
			throw std::runtime_error("Only MatDense and MatSparse factors are handled by sliceMultiply op for the moment.");
		}
		delete [] tmp_ptr;
	}
#else
	throw std::runtime_error("TransformHelper::sliceMultiply is not supported with eigen version < 3.9");
#endif
	return out;
}

#include "faust_TransformHelper_cat.hpp"
