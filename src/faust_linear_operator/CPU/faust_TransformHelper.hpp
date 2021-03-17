/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifdef USE_GPU_MOD
		if(FaustGPU<FPP>::are_cpu_mat_all_known(this->transform->data))
			//with refman gpu mats will be reused (because already loaded)
			enable_gpu_meth_for_mul();
#endif
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : TransformHelperGen<FPP,Cpu>()
#ifdef USE_GPU_MOD
													  , gpu_faust(nullptr)
#endif
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
#ifdef USE_GPU_MOD
			if(th_left->gpu_faust && th_right->gpu_faust)
				enable_gpu_meth_for_mul();
#endif
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate) : Faust::TransformHelper<FPP,Cpu>()
	{
		this->transform = th->transform;
		this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
		this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		this->copy_slice_state(*th);
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
#ifdef USE_GPU_MOD
			if(th.gpu_faust)
				enable_gpu_meth_for_mul();
#endif
		}

#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
	template<typename FPP>
		template<typename ... GList>
		TransformHelper<FPP,Cpu>::TransformHelper(GList& ... t) : TransformHelper<FPP,Cpu>()
		{
			this->push_back_(t...);
#ifdef USE_GPU_MOD
			// it's difficult here to test if all t are TransformHelper with a gpu_faust enabled (for example t can also be a vector<MatGeneric>)
			// so verify directly that cpu matrices copied are loaded in GPU to respect the invariant: a Faust is fully loaded on GPU or is not on GPU at all
			if(FaustGPU<FPP>::are_cpu_mat_all_known(this->transform->data))
				enable_gpu_meth_for_mul();
#endif
		}
#endif
	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatSparse<FPP,Cpu> A, const bool transpose /* deft to false */, const bool conjugate)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			MatDense<FPP,Cpu> M = this->transform->multiply(A, this->isTransposed2char());
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return M;
		}


	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> &x, const bool transpose, const bool conjugate)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
#ifdef USE_GPU_MOD
			if(this->Fv_mul_mode == GPU_MOD && gpu_faust != nullptr)
					return gpu_faust->multiply(x, this->is_transposed, this->is_conjugate);
#endif
			Vect<FPP,Cpu> v = std::move(this->transform->multiply(x, this->isTransposed2char()));
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			return v;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const FPP *x, const bool transpose, const bool conjugate)
		{
			int x_size;
			// assuming that x size is valid, infer it from this size
			if(this->is_transposed && transpose || ! this->is_transposed && ! transpose)
				x_size = this->getNbCol();
			else
				x_size = this->getNbRow();
			Vect<FPP, Cpu> vx(x_size, x);
			return std::move(this->multiply(vx, transpose, conjugate));
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::multiply(const FPP *x, FPP* y, const bool transpose, const bool conjugate)
		{
			int x_size;
			// assuming that x size is valid, infer it from this size
			if(this->is_transposed && transpose || ! this->is_transposed && ! transpose)
				x_size = this->getNbCol();
			else
				x_size = this->getNbRow();
			Vect<FPP, Cpu> vx(x_size, x);
			auto y_vec = std::move(this->multiply(vx, transpose, conjugate));
			memcpy(y, y_vec.getData(), sizeof(FPP)*y_vec.size());
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> &A, const bool transpose, const bool conjugate)
		{
			this->is_transposed ^= transpose;
			this->is_conjugate ^= conjugate;
			Faust::MatDense<FPP,Cpu> M;
#ifdef FAUST_TORCH
			if(this->mul_order_opt_mode >= 7 && this->mul_order_opt_mode <= 9 && !tensor_data.size())
			{
				// init tensor data cache
				convMatGenListToTensorList(this->transform->data, tensor_data, at::kCPU, /* clone */ false, /* transpose */ ! this->is_transposed);
//				display_TensorList(tensor_data);
			}
#endif

			switch(this->mul_order_opt_mode)
			{
				case ORDER_ALL_ENDS:
				case ORDER_1ST_BEST:
				case ORDER_ALL_BEST_CONVDENSE:
				case ORDER_ALL_BEST_MIXED:
					{
						//					this->transform->data.push_back(const_cast<Faust::MatDense<FPP,Cpu>*>(&A)); // it's ok
						std::vector<Faust::MatGeneric<FPP,Cpu>*> data(this->transform->size()+1);
						if(transpose)
						{
							int i = this->transform->size();
							for(auto fac: this->transform->data)
								data[i--] = fac;
							data[0] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
						}
						else
						{
							int i = 0;
							for(auto fac: this->transform->data)
								data[i++] = fac;
							data[i] = const_cast<Faust::MatDense<FPP,Cpu>*>(&A);
						}
						Faust::multiply_order_opt(this->mul_order_opt_mode, data, M, /*alpha */ FPP(1.0), /* beta */ FPP(0.0), {this->isTransposed2char()});
						//					this->transform->data.erase(this->transform->data.end()-1);
					}
					break;
				case CPP_PROD_PAR_REDUC:
					M = Faust::multiply_par(this->transform->data, A, this->isTransposed2char());
					break;
				case OMP_PROD_PAR_REDUC:
					M = Faust::multiply_omp(this->transform->data, A, this->isTransposed2char());
					break;
#ifdef FAUST_TORCH
				case TORCH_CPU:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_BEST_ORDER:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false,  /*clone */ false,/* chain_opt */ true, /* contiguous_dense_to_torch */ false, !this->is_transposed);
					break;
				case TORCH_CPU_DENSE_ROW_TO_TORCH:
					Faust::tensor_chain_mul(tensor_data, M, &A, /* on_gpu */ false, /*clone */ false, /* chain_opt */ false, /* contiguous_dense_to_torch */ true, !this->is_transposed);
					break;
#endif
#ifdef USE_GPU_MOD
				case GPU_MOD:
					if(gpu_faust != nullptr)
						M = gpu_faust->multiply(&A, this->is_transposed, this->is_conjugate);
					break;
#endif
				default:
					M = this->transform->multiply(A, this->isTransposed2char());
					break;
			}
			this->is_conjugate ^= conjugate;
			this->is_transposed ^= transpose;
			return M;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize(const bool transp /* deft to false */)
		{
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
			// choose the quickest method for the Faust "toarray"
			auto t = this->optimize_time_full(transp, inplace, nsamples);
			// choose the quickest method for the Faust-vector mul
			t = t->optimize_time_Fv(transp, inplace, nsamples);
			return t;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_time_full(const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples)
		{
			return this->optimize_multiply([this](){this->get_product();}, transp, inplace, nsamples, "Faust-toarray");
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_time_Fv(const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples)
		{
			// generate a random vector
			Faust::Vect<FPP,Cpu> *v = nullptr;
			if(transp)
				v = Faust::Vect<FPP,Cpu>::rand(this->getNbRow());
			else
				v = Faust::Vect<FPP,Cpu>::rand(this->getNbCol());
			return this->optimize_multiply([this, &v, &transp, &inplace, &nsamples]()
					{
						this->multiply(*v, transp);
					}, transp, inplace, nsamples, "Faust-vector mul.");
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::optimize_multiply(std::function<void()> f, const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples, const char* op_name)
		{
			TransformHelper<FPP,Cpu>* t_opt = nullptr;
			int NMETS = 11;
			std::chrono::duration<double> * times = new std::chrono::duration<double>[NMETS]; //use heap because of VS14 (error C3863)
//			MatDense<FPP,Cpu>* M = MatDense<FPP,Cpu>::randMat(transp?this->getNbRow():this->getNbCol(), 2048);
			int old_meth = this->get_mul_order_opt_mode();
			int nmuls = nsamples, opt_meth=0;
			std::vector<int> disabled_meths = {CPP_PROD_PAR_REDUC, OMP_PROD_PAR_REDUC}; // skip openmp/C++ threads methods because they are unfruitful when Eigen is multithreaded
#if DEBUG_OPT_MUL
			cout << "nsamples used to measure time: " << nmuls << endl;
#endif
#ifdef FAUST_TORCH
			if(this->mul_order_opt_mode >= 7 && this->mul_order_opt_mode <= 9 && !tensor_data.size())
			{
				// init tensor data cache
				convMatGenListToTensorList(this->transform->data, tensor_data, at::kCPU, /* clone */ false, /* transpose */ ! this->is_transposed);
//				display_TensorList(tensor_data);
			}
#else
			for(int i=7;i<GPU_MOD;i++)
				disabled_meths.push_back(i);
#endif
#ifdef USE_GPU_MOD
			try
			{
				GPUModHandler::get_singleton()->check_gpu_mod_loaded();
			}
			catch(std::runtime_error& e)
			{
			// test only if enable_gpu_mod was called
#endif
				disabled_meths.push_back(GPU_MOD);
#ifdef USE_GPU_MOD
			}
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
					//				auto FM = this->multiply(*M, transp);
//					this->get_product();
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
				cout << "best method measured in time on operation "<< op_name << " is: " << opt_meth << endl;
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
			return t_opt;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::pruneout(const int nnz_tres, const int npasses, const bool only_forward)
		{
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
								else
								{
									tmp_ds = dynamic_cast<Faust::MatDense<FPP,Cpu>*>(S_j);
									if(nnz_i <= nnz_tres)
										tmp_ds->delete_row(offset);
								}
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
			pth->transform->update_total_nnz();
			pth->copy_transconj_state(*this);
			return pth;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::update_total_nnz()
		{
			this->transform->update_total_nnz();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::enable_gpu_meth_for_mul()
		{
#ifdef USE_GPU_MOD
			set_FM_mul_mode(GPU_MOD);
			set_Fv_mul_mode(GPU_MOD);
#endif
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::set_FM_mul_mode(const int mul_order_opt_mode, const bool silent /* = false */)
		{
#ifdef USE_GPU_MOD
			if(mul_order_opt_mode == GPU_MOD && gpu_faust == nullptr)
			{
//				try
//				{
					gpu_faust = new Faust::FaustGPU<FPP>(this->transform->data);
//				}
//				catch(std::runtime_error & e)
//				{
//					// creation failed, nothing todo (must be because of cuda backend is not installed/found)
//					std::cerr << "error: can't change to GPU method (backend unavailable)." << std::endl;
//					return;
//				}
			}
#endif
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
		void TransformHelper<FPP,Cpu>::set_Fv_mul_mode(const int Fv_mul_mode)
		{
#ifdef USE_GPU_MOD
			//TODO: factorize this code with set_FM_mul_mode
			if(this->Fv_mul_mode == GPU_MOD && gpu_faust == nullptr)
			{
				try
				{
					gpu_faust = new Faust::FaustGPU<FPP>(this->transform->data);
				}
				catch(std::runtime_error & e)
				{
					// creation failed, nothing todo (must be because of cuda backend is not installed/found)
					std::cerr << "error: can't change to GPU method (backend unavailable)." << std::endl;
					return;
				}
			}
#endif
			this->Fv_mul_mode = Fv_mul_mode;
			std::cout << "changed Faust-vector mul. mode to: " << this->Fv_mul_mode;
			if(! this->Fv_mul_mode)
				std::cout << " (opt. disabled, default mul.)";
			std::cout << std::endl;
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(const TransformHelper<FPP, Cpu>* th_right) const
		{
			return new TransformHelper<FPP,Cpu>(this, th_right);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(FPP& scalar)
		{
			const vector<MatGeneric<FPP,Cpu>*>& vec = this->transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
			//the point here is to minimize the number of copies (with direct access)
			// the constructor then will copy the factors from the vector
//			Transform<FPP,Cpu>* t = new Transform<FPP,Cpu>(vec, scalar, false, true); //optimizedCopy == false, cloning_fact == true
//			TransformHelper<FPP,Cpu>* th  = new TransformHelper<FPP,Cpu>(*t);
			TransformHelper<FPP,Cpu>* th = new TransformHelper<FPP,Cpu>(vec, scalar, false, false, true);
			th->copy_transconj_state(*this);
			th->copy_slice_state(*this);
			return th;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::pop_back()
        {
            this->transform->pop_back();
#ifdef USE_GPU_MOD
			if(gpu_faust)
				gpu_faust->pop_back();
#endif
        }

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::pop_front()
        {
            this->transform->pop_front();
#ifdef USE_GPU_MOD
			if(gpu_faust)
				gpu_faust->pop_front();
#endif
        }

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy /* false by deft */, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			auto sparse_mat = new MatSparse<FPP,Cpu>(nnz, nrows, ncols, data, row_ptr, id_col);
			this->push_back(sparse_mat, optimizedCopy, false, transpose, conjugate);
			if(optimizedCopy) delete sparse_mat;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /* false by default */, const bool copying /* true to default */, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			//warning: should not be called after initialization of factors (to respect the immutability property)
			//this function is here only for python wrapper (TODO: see how to modify that wrapper in order to delete this function after or just use it internally -- not py/matfaust)
			//TODO: transpose and conjugate must be passed to transform and taken into account for the factor (like it's done in GPU2 impl., i.e.: only possible if copying is true too)
			this->transform->push_back(M, optimizedCopy, this->is_conjugate, copying); //2nd argument is for opt. (possibly converting dense <-> sparse)
            if(transpose)
                get_gen_fact_nonconst(size()-1)->transpose();
#ifdef USE_GPU_MOD
			if(gpu_faust)
				gpu_faust->push_back(const_cast<MatGeneric<FPP,Cpu>*>(M));
#endif
		}
    
	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_first(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy /* false by default */, const bool copying /* true to default */)
		{
			//warning: should not be called after initialization of factors (to respect the immutability property)
			//this function is here only for python wrapper (TODO: see how to modify that wrapper in order to delete this function after or just use it internally -- not py/matfaust)
			this->transform->push_first(M, optimizedCopy, this->is_conjugate, copying); //2nd argument is for opt. (possibly converting dense <-> sparse)
#ifdef USE_GPU_MOD
			if(gpu_faust)
				gpu_faust->insert(M, 0);
#endif
		}

	template<typename FPP>
		template<typename Head, typename ... Tail>
		void TransformHelper<FPP,Cpu>::push_back_(Head& h, Tail&... t)
		{
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
			faust_unsigned_int nbytes = 0;
			for(auto fac : this->transform->data)
			{
				if(dynamic_cast<Faust::MatDense<FPP, Cpu>*>(fac))
					nbytes += fac->getNbCol() * fac->getNbRow() * sizeof(FPP);
				else if (dynamic_cast<Faust::MatSparse<FPP, Cpu>*>(fac))
					nbytes += fac->getNonZeros() * (sizeof(FPP) + sizeof(int)) + (fac->getNbRow() + 1) * sizeof(int); // by default storage index is int
				else if (dynamic_cast<Faust::MatDiag<FPP>*>(fac))
					nbytes += sizeof(FPP) * (fac->getNbCol()<fac->getNbRow()?fac->getNbCol():fac->getNbRow());
				else
					throw runtime_error("Unknown matrix type.");
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
//			cout << "is_transposed:" << this->is_transposed << std::endl;
			this->transform->Display(this->is_transposed, false);
		}

	template<typename FPP>
		std::string TransformHelper<FPP,Cpu>::to_string() const
		{
			return this->transform->to_string(this->is_transposed);
		}

	//private
	template<typename FPP>
	const MatGeneric<FPP,Cpu>* TransformHelper<FPP,Cpu>::get_gen_fact(const faust_unsigned_int id) const
	{
		return this->transform->data[this->is_transposed?size()-id-1:id];
	}

template<typename FPP>
	MatGeneric<FPP,Cpu>* TransformHelper<FPP,Cpu>::get_gen_fact_nonconst(const faust_unsigned_int id) const
	{
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
#if USE_GPU_MOD
		if(gpu_faust)
		{
			gpu_faust->update(fact, fact_id);
		}
#endif
	}

template<typename FPP>
	unsigned long long TransformHelper<FPP,Cpu>::get_fact_addr(const faust_unsigned_int id) const
	{
		return (unsigned long long) this->transform->data[this->is_transposed?size()-id-1:id];
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
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, elts, num_rows, num_cols);
	}

template<typename FPP>
	MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
	{
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
		this->transform->get_fact(this->is_transposed?this->size()-id-1:id, elts, num_rows, num_cols, this->is_transposed ^ transpose);
		if(this->is_conjugate)
			Faust::conjugate(elts,*num_cols*(*num_rows));
	}




template<typename FPP>
	MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product(const int mul_order_opt_mode/*=DEFAULT*/) // const
	{
		int old_FM_mul_mod = -1;
		MatDense<FPP, Cpu> P;
		// keep current mul mode to restore it on the function end
		if(mul_order_opt_mode != this->mul_order_opt_mode)
		{
			old_FM_mul_mod = this->mul_order_opt_mode;
			this->set_FM_mul_mode(mul_order_opt_mode);
		}
		if(this->mul_order_opt_mode)
			switch(this->mul_order_opt_mode)
			{
#ifdef USE_GPU_MOD
				case GPU_MOD:
					if(nullptr != gpu_faust)
						P = gpu_faust->get_product(this->is_transposed, this->is_conjugate);
#endif
				default:
					//TODO: avoid to add one factor for all methods if possible
					MatDense<FPP,Cpu> Id(this->getNbCol(), this->getNbCol());
					Id.setEyes();
					P = this->multiply(Id);
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
	void TransformHelper<FPP,Cpu>::get_product(Faust::MatDense<FPP,Cpu>& prod, const int mul_order_opt_mode/*=DEFAULT*/) //const
	{
		if(mul_order_opt_mode != DEFAULT)
			prod = this->get_product(mul_order_opt_mode);
		else
			this->transform->get_product(prod, this->isTransposed2char(), this->is_conjugate);
	}

template<typename FPP>
	void TransformHelper<FPP,Cpu>::save_mat_file(const char* filepath) const
	{
		this->transform->save_mat_file(filepath, this->is_transposed, this->is_conjugate);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
	{
//			std::cout << "TransformHelper<FPP,Cpu>::spectralNorm" << std::endl;
#ifdef USE_GPU_MOD
		if(gpu_faust != nullptr)
			return gpu_faust->spectral_norm(nbr_iter_max, threshold);
#endif
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
	double TransformHelper<FPP,Cpu>::normL1() const {
		return this->transform->normL1(this->is_transposed);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::normInf() const {
		return this->transform->normL1(!this->is_transposed);
	}

template<typename FPP>
	double TransformHelper<FPP,Cpu>::normFro() const {
		vector <MatGeneric<FPP, Cpu>*>& orig_facts = this->transform->data;
		int start_id, end_id;
		this->transform->get_nonortho_interior_prod_ids(start_id, end_id);
		if(start_id < 0)
			return Faust::fabs(MatDense<FPP,Cpu>::eye(this->getNbCol(), this->getNbCol()).norm());
		else if(start_id == 0)
			return this->transform->normFro();
		vector<MatGeneric<FPP,Cpu>*> non_ortho_start(orig_facts.begin()+start_id, orig_facts.end());
		TransformHelper<FPP, Cpu> t(non_ortho_start, 1.0, false, false);
		return t.transform->normFro();
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const
	{
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
				norm_invs[i] = 1./norm;
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
#ifdef USE_GPU_MOD
		if(gpu_faust != nullptr)
			delete gpu_faust;
#endif
	}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density /* 1.f */, bool per_row /* true */)
	{
		return TransformHelper<FPP,Cpu>::randFaust(-1, -1, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
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
			//			for(int i=0;i<n+1;i++)
			//				factors[i]->Display();
			float alpha = norma?1/sqrt((double)(1 << n)):1.0;
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
		return this->transform->begin();
	}

template<typename FPP>
	Faust::transf_iterator<FPP> Faust::TransformHelper<FPP, Cpu>::end() const
	{
		return this->transform->end();
	}

template<typename FPP>
	void Faust::TransformHelper<FPP,Cpu>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode/*=DEFAULT*/)
	{
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
	void TransformHelper<FPP,Cpu>::pack_factors(const int mul_order_opt_mode/*=DEFAULT*/)
{
	TransformHelperGen<FPP,Cpu>::pack_factors(mul_order_opt_mode/*=DEFAULT*/);
}

template<typename FPP>
	TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::swap_rows(const faust_unsigned_int id1,
			const faust_unsigned_int id2,
			const bool permutation/*=false*/,
			const bool inplace/*=false*/,
			const bool check_transpose/*=true*/)
{
	if(check_transpose && this->is_transposed)
		return swap_cols(id1, id2, permutation, inplace, false);
	TransformHelper<FPP,Cpu>* t = nullptr;
	if(inplace)
		t = this;
	else
	{
		t = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0,
				/* optimizedCopy*/ false, /* cloning_fact */ true, /* internal_call*/ true);
		t->copy_transconj_state(this);
		t->copy_slice_state(this);
		t->copy_mul_mode_state(this);
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
	if(check_transpose && this->is_transposed)
		return swap_rows(id1, id2, permutation, inplace, false);
	TransformHelper<FPP,Cpu>* t = nullptr;
	if(inplace)
		t = this;
	else
	{
		t = new TransformHelper<FPP,Cpu>(this->transform->data, 1.0,
				/* optimizedCopy*/ false, /* cloning_fact */ true, /* internal_call*/ true);
		t->copy_transconj_state(this);
		t->copy_slice_state(this);
		t->copy_mul_mode_state(this);
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

	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));

}

#include "faust_TransformHelper_cat.hpp"
