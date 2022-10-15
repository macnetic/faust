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

#ifndef __FAUST_TRANSFORM_HELPER___
#define __FAUST_TRANSFORM_HELPER___

#define NOMINMAX // avoids VS min/max issue with std::min/max.
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
#include "faust_TransformHelper_cat.h"

namespace Faust
{

	template<typename FPP>
		using transf_iterator = typename Transform<FPP,Cpu>::transf_iterator;


	template<typename FPP>
		class TransformHelper<FPP,Cpu> : public TransformHelperGen<FPP,Cpu>
		{

			friend TransformHelper<FPP,Cpu>* vertcat<FPP>(const std::vector<TransformHelper<FPP,Cpu>*>& THs);
			static std::default_random_engine generator;
			static bool seed_init;

#ifdef FAUST_TORCH
			std::vector<torch::Tensor> tensor_data;
#endif

			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
			TransformHelper();
			TransformHelper(const TransformHelper<FPP,Cpu>* th_left, const TransformHelper<FPP,Cpu>* th_right);
			TransformHelper(TransformHelper<FPP,Cpu>* th); //TODO: it shouldn't be necessary to have a copy ctor on pointer and another one on reference
			TransformHelper(TransformHelper<FPP,Cpu>& th);
			void operator=(TransformHelper<FPP,Cpu>& th);
			TransformHelper(const TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
			TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]);
			TransformHelper(TransformHelper<FPP,Cpu>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			TransformHelper(Transform<FPP,Cpu> &t, const bool moving=false);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
			template<typename ...GList>
				TransformHelper(GList& ... t);
#endif

			void copy_mul_mode_state(const TransformHelper<FPP,Cpu>& th);
			virtual Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &x);//TODO: should be const and child redefs too
			virtual Vect<FPP,Cpu> multiply(const FPP* x);//TODO: should be const and child redefs too
			virtual void multiply(const FPP* x, FPP* y); //TODO: should be const and child redefs too
			// Multiplies the slice s of this by the vector x corresponding elements (x size is this->getNbCol())
			FPP* sliceMultiply(const Slice s[2], const FPP* X, FPP* out=nullptr, int X_ncols=1) const;
			// Note: Prefers the prototype above as it avoids copies
			MatDense<FPP, Cpu> sliceMultiply(const Slice s[2], const FPP* X, int X_ncols=1) const;
			// Multiplies the columns ids of this by the vector x corresponding elements (x size is this->getNbCol())
			Vect<FPP, Cpu> indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* x) const;
			MatDense<FPP, Cpu> indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* X, int ncols) const;
			FPP* indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* x, FPP* out) const;
			FPP* indexMultiply(faust_unsigned_int* ids[2], size_t ids_len[2], const FPP* X, int ncols, FPP* out) const;
			// \brief multiply this by A (of size: this->getNbCol()*A_ncols) into C (buffers must be properly allocated from the callee).
			virtual void multiply(const FPP* A, int A_ncols, FPP* C);
			//			MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> A) const;
			virtual MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &A);
			virtual MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A);

			virtual TransformHelper<FPP, Cpu>* multiply(const TransformHelper<FPP, Cpu>*) const;
			virtual TransformHelper<FPP, Cpu>* multiply(FPP& scalar);
			template<typename Head, typename ... Tail>
				void push_back_(Head& h, Tail&... t);
			//
			void push_back_();
			void push_back(const FPP* data, int nrows, int ncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);
			void push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);
			void push_back(const FPP* bdata, const int* brow_ptr, const int* bcol_inds, const int nrows, const int ncols, const int bnnz, const int bnrows, const int bncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);

			void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true, const bool transpose=false, const bool conjugate=false);
			void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy, const bool copying, const bool transpose, const bool conjugate, const bool verify_dims_agree);
			void pop_back();
			void pop_front();
			void push_first(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true);
			virtual faust_unsigned_int getNBytes() const;
			virtual faust_unsigned_int get_total_nnz() const;
			virtual void update_total_nnz();
			bool is_zero() const;
			faust_unsigned_int size() const;
			virtual void resize(faust_unsigned_int);
			void display() const;
			virtual string to_string() const;
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
			void get_fact_bsr_info(const faust_unsigned_int id,
					size_t& bdata_sz,
					size_t& browptr_sz,
					size_t& bcolinds_sz,
					size_t& bnnz,
					size_t& bnrows,
					size_t& bncols) const;
			void get_fact(const faust_unsigned_int id,
					FPP* bdata,
					int* brow_ptr,
					int* bcol_inds) const;
			virtual MatDense<FPP,Cpu> get_product(const int mul_order_opt_mode=-1);// const;
			virtual void get_product(MatDense<FPP,Cpu>& prod, const int mul_order_opt_mode=-1); //const;
			virtual void save_mat_file(const char* filename) const;
            /**
             * \brief Reads a TransformHelper from a .mat file.
             */
            void read_from_mat_file(const char* filepath);
            /**
             * \brief Returns the type of the Faust::Transform stored in the .mat
             * file (0: float, 1: double, 2: complex, -1: unknown/invalid).
             */
            static int get_mat_file_type(const char* filepath);
			virtual double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
			FPP power_iteration(const faust_unsigned_int nbr_iter_max, const Real<FPP>& threshold, int & flag) const;
			virtual TransformHelper<FPP,Cpu>* transpose();
			TransformHelper<FPP,Cpu>* conjugate();
			TransformHelper<FPP,Cpu>* adjoint() const;
			virtual TransformHelper<FPP,Cpu>* vertcat(const TransformHelper<FPP,Cpu>*);
			virtual TransformHelper<FPP,Cpu>* horzcat(const TransformHelper<FPP,Cpu>*);
			virtual double normL1(const bool full_array=true, const int batch_size=1) const;
			virtual double normFro(const bool full_array=true, const int batch_size=1) const;
			virtual double normInf(const bool full_array=true, const int batch_size=1) const;
			virtual TransformHelper<FPP,Cpu>* normalize(const int meth = 2/* 1 for 1-norm, 2 for 2-norm, MAX for inf-norm */) const;
			/**
			 * \param only_forward: True for applying only the forward passes of removal.
			 * \param npasses: the number of passes to run, by default it goes until the optimal Faust is obtained.
			 */
			virtual TransformHelper<FPP,Cpu>* pruneout(const int nnz_tres, const int npasses=-1, const bool only_forward=false);
			virtual TransformHelper<FPP,Cpu>* optimize_multiply(std::function<void()> f, const bool transp=false, const bool inplace=false, const int nsamples=1, const char* op_name="unamed_op");

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

			void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode=DEFAULT_L2R);
			void pack_factors(const int mul_order_opt_mode=DEFAULT_L2R);

			virtual TransformHelper<FPP,Cpu>* swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);
			virtual TransformHelper<FPP,Cpu>* swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);

			/** for testing purpose only (memory leaks enabled) */
			void disable_dtor() { this->transform->disable_dtor(); }
			void enable_dtor() { this->transform->enable_dtor(); }
			static TransformHelper<FPP, Cpu>* randBSRFaust(unsigned int faust_nrows, unsigned int faust_ncols, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int bnrows, unsigned int bncols, float density=.1f);
			static TransformHelper<FPP,Cpu>* randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true);
			static TransformHelper<FPP,Cpu>* randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true);
			static TransformHelper<FPP,Cpu>* hadamardFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* eyeFaust(unsigned int n, unsigned int m);

			virtual ~TransformHelper();
			unsigned long long get_fact_addr(const faust_unsigned_int id) const;
			/* use carefully */
			MatGeneric<FPP,Cpu>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
			const MatGeneric<FPP,Cpu>* get_gen_fact(const faust_unsigned_int id) const;
			void update(const MatGeneric<FPP,Cpu>& M, const faust_unsigned_int fact_id);
			void replace(const MatGeneric<FPP, Cpu>* M, const faust_unsigned_int fact_id);
			void convertToSparse();
			void convertToDense();
			FPP get_item(faust_unsigned_int i, faust_unsigned_int j);
			template<typename FPP2>
			TransformHelper<Real<FPP2>, Cpu>* real();
		};


}

#include "faust_TransformHelper.hpp"
#endif
