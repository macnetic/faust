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

#ifndef __FAUST_TRANSFORM_HELPER___
#define __FAUST_TRANSFORM_HELPER___

#include <memory>
#include "faust_TransformHelperGen.h"
#include "faust_RefManager.h"
#include "faust_exception.h"
#include "faust_Transform.h"
#include "faust_Slice.h"
#include "faust_MatSparse.h"
#include <random>
#include <limits>
#include <algorithm>
#include <vector>
#ifdef FAUST_TORCH
#include "faust_torch.h"
#endif
#ifdef USE_GPU_MOD
#include "faust_gpu_mod.h"
#endif
#include "faust_TransformHelper_cat.h"

namespace Faust
{

	template<typename FPP>
		using transf_iterator = typename Transform<FPP,Cpu>::transf_iterator;
#ifdef USE_GPU_MOD
	template<typename FPP> class FaustGPU;
#endif


	template<typename FPP>
		class TransformHelper<FPP,Cpu> : public TransformHelperGen<FPP,Cpu> {

			friend TransformHelper<FPP,Cpu>* vertcat<FPP>(const std::vector<TransformHelper<FPP,Cpu>*>& THs);
			static std::default_random_engine generator;
			static bool seed_init;

#ifdef FAUST_TORCH
			std::vector<torch::Tensor> tensor_data;
#endif
#ifdef USE_GPU_MOD
			FaustGPU<FPP> *gpu_faust;
#endif
			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
			TransformHelper();
			TransformHelper(const TransformHelper<FPP,Cpu>* th_left, const TransformHelper<FPP,Cpu>* th_right);
			TransformHelper(TransformHelper<FPP,Cpu>* th);
			TransformHelper(TransformHelper<FPP,Cpu>& th);
			void operator=(TransformHelper<FPP,Cpu>& th);
			TransformHelper(const TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
			TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]);
			TransformHelper(TransformHelper<FPP,Cpu>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			TransformHelper(Transform<FPP,Cpu> &t, const bool moving=false);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
			template<typename ...GList> TransformHelper(GList& ... t);
#endif

			void copy_mul_mode_state(const TransformHelper<FPP,Cpu>& th);
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x, const bool transpose=false, const bool conjugate=false);
			//			MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> A) const;
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> A, const bool transpose=false, const bool conjugate=false);
			void update_total_nnz();
			void set_FM_mul_mode(const int mul_order_opt_mode, const bool silent=false);
			void set_Fv_mul_mode(const int mode);
			void enable_gpu_meth_for_mul();
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> A, const bool transpose=false, const bool conjugate=false);

			TransformHelper<FPP, Cpu>* multiply(const TransformHelper<FPP, Cpu>*) const;
			TransformHelper<FPP, Cpu>* multiply(FPP& scalar);
			template<typename Head, typename ... Tail>
				void push_back_(Head& h, Tail&... t);
			//
			void push_back_();
			void push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);

			void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true, const bool transpose=false, const bool conjugate=false);
			void pop_back();
			void pop_front();
			void push_first(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true);
			faust_unsigned_int getNBytes() const;
			bool isReal() const;
			faust_unsigned_int get_total_nnz() const;
			faust_unsigned_int size() const;
			void resize(faust_unsigned_int);
			void display() const;
			string to_string() const;
			MatDense<FPP,Cpu> get_fact(faust_unsigned_int id) const;
			void get_fact(const faust_unsigned_int id,
					const int** rowptr,
					const int** col_ids,
					const FPP** elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols) const;
			void get_fact(const faust_unsigned_int id,
					int* rowptr,
					int* col_ids,
					FPP* elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose=false) const;
			void get_fact(const faust_unsigned_int id,
					const FPP** elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols) const;
			void get_fact(const faust_unsigned_int &id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			MatDense<FPP,Cpu> get_product();// const;
			void get_product(MatDense<FPP,Cpu>& prod) const;
			void save_mat_file(const char* filename) const;
			double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
			TransformHelper<FPP,Cpu>* transpose();
			TransformHelper<FPP,Cpu>* conjugate();
			TransformHelper<FPP,Cpu>* adjoint() const;
			TransformHelper<FPP,Cpu>* vertcat(const TransformHelper<FPP,Cpu>*);
			TransformHelper<FPP,Cpu>* horzcat(const TransformHelper<FPP,Cpu>*);
			bool isConjugate() const;
			double normL1() const;
			double normFro() const;
			double normInf() const;
			TransformHelper<FPP,Cpu>* normalize(const int meth = 2/* 1 for 1-norm, 2 for 2-norm, MAX for inf-norm */) const;
			/**
			 * \param only_forward: True for applying only the forward passes of removal.
			 * \param npasses: the number of passes to run, by default it goes until the optimal Faust is obtained.
			 */
			TransformHelper<FPP,Cpu>* pruneout(const int nnz_tres, const int npasses=-1, const bool only_forward=false);
			TransformHelper<FPP,Cpu>* optimize(const bool transp=false);
			TransformHelper<FPP,Cpu>* optimize_multiply(std::function<void()> f, const bool transp=false, const bool inplace=false, const int nsamples=1, const char* op_name="unamed_op");
			TransformHelper<FPP,Cpu>* optimize_time(const bool transp=false, const bool inplace=false, const int nsamples=1);
			TransformHelper<FPP,Cpu>* optimize_time_full(const bool transp=false, const bool inplace=false, const int nsamples=1);
			TransformHelper<FPP,Cpu>* optimize_time_Fv(const bool transp=false, const bool inplace=false, const int nsamples=1);
			/**
			  \brief Returns the left hand side factors of this from index 0 to id included (as a new TransformHelper obj).

*/
//			TransformHelper<FPP,Cpu>* left(const faust_unsigned_int id, const bool copy=false) const;
			/**
			  \brief Returns the right hand side factors of this from index id to the size()-1 (as a new TransformHelper obj).

*/
//			TransformHelper<FPP,Cpu>* right(const faust_unsigned_int id, const bool copy=false) const;

			transf_iterator<FPP> begin() const;
			transf_iterator<FPP> end() const;

#ifdef USE_GPU_MOD
			FaustGPU<FPP>* get_gpu_faust();
#endif

			void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id);
			void pack_factors();

			TransformHelper<FPP,Cpu>* swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);
			TransformHelper<FPP,Cpu>* swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);

			/** for testing purpose only (memory leaks enabled) */
			void disable_dtor() { this->transform->disable_dtor(); }
			void enable_dtor() { this->transform->enable_dtor(); }
			static TransformHelper<FPP,Cpu>* randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true);
			static TransformHelper<FPP,Cpu>* randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true);
			static TransformHelper<FPP,Cpu>* hadamardFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* eyeFaust(unsigned int n, unsigned int m);

			~TransformHelper();
			unsigned long long get_fact_addr(const faust_unsigned_int id) const;
			/* use carefully */
			MatGeneric<FPP,Cpu>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
			void update(const MatGeneric<FPP,Cpu>& M, const faust_unsigned_int fact_id);
			private:
			const MatGeneric<FPP,Cpu>* get_gen_fact(const faust_unsigned_int id) const;

		};


}

#include "faust_TransformHelper.hpp"
#endif
