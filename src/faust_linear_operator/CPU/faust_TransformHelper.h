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
#include "faust_RefManager.h"
#include "faust_exception.h"
#include "faust_Transform.h"
#include "faust_Slice.h"
#include "faust_MatSparse.h"
#include <random>

namespace Faust {
	using namespace std;

	template<typename FPP>
		using transf_iterator = typename Transform<FPP,Cpu>::transf_iterator;
	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class MatGeneric;

	enum RandFaustType {
		DENSE,
		SPARSE,
		MIXED
	};

	enum PackDir {
		PACK_LEFT,
		PACK_RIGHT
	};

	template<typename FPP>
		class TransformHelper<FPP,Cpu> {
			static std::default_random_engine generator;
			static bool seed_init;

			bool is_transposed;
			bool is_conjugate;
			bool is_sliced;
			Slice slices[2];
			bool is_fancy_indexed;
			int mul_order_opt_mode;
			faust_unsigned_int * fancy_indices[2];
			faust_unsigned_int fancy_num_rows;
			faust_unsigned_int fancy_num_cols;
			shared_ptr<Transform<FPP,Cpu>> transform;

			void eval_sliced_Transform();
			void eval_fancy_idx_Transform();
			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
			TransformHelper();
			TransformHelper(TransformHelper<FPP,Cpu>* th_left, TransformHelper<FPP,Cpu>* th_right);
			TransformHelper(TransformHelper<FPP,Cpu>* th);
			TransformHelper(TransformHelper<FPP,Cpu>& th);
			void operator=(TransformHelper<FPP,Cpu>& th);
			TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
			TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]);
			TransformHelper(TransformHelper<FPP,Cpu>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			TransformHelper(Transform<FPP,Cpu> &t, const bool moving=false);
			template<typename ...GList> TransformHelper(GList& ... t);
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x, const bool transpose=false, const bool conjugate=false);
//			MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> A) const;
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> A, const bool transpose=false, const bool conjugate=false);
			void update_total_nnz();
			void set_mul_order_opt_mode(const int mul_order_opt_mode);
			int get_mul_order_opt_mode() const;
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> A, const bool transpose=false, const bool conjugate=false);

			TransformHelper<FPP, Cpu>* multiply(TransformHelper<FPP, Cpu>*);
			TransformHelper<FPP, Cpu>* multiply(FPP& scalar);
			template<typename Head, typename ... Tail>
				void push_back_(Head& h, Tail&... t);

			void push_back_();
			void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true);
            void pop_back();
            void pop_front();
            void push_first(const Faust::MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const bool copying=true);
			faust_unsigned_int getNbRow() const;
			faust_unsigned_int getNbCol() const;
			faust_unsigned_int getNBytes() const;
			bool isReal() const;
			faust_unsigned_int get_total_nnz() const;
			faust_unsigned_int size() const;
			void resize(faust_unsigned_int);
			void display() const;
			string to_string() const;
			faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			unsigned int get_fact_nb_rows(const faust_unsigned_int id) const;
			unsigned int get_fact_nb_cols(const faust_unsigned_int id) const;
			unsigned int get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const;
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
			void get_fact(const faust_unsigned_int id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			bool is_fact_sparse(const faust_unsigned_int id) const;
			bool is_fact_dense(const faust_unsigned_int id) const;
			TransformHelper<FPP, Cpu>* slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
					faust_unsigned_int start_col_id, faust_unsigned_int end_col_id);
			TransformHelper<FPP, Cpu>* fancy_index(faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
			MatDense<FPP,Cpu> get_product() const;
            void get_product(Faust::MatDense<FPP,Cpu>& prod) const;
			void save_mat_file(const char* filename) const;
			double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
			TransformHelper<FPP,Cpu>* transpose();
			TransformHelper<FPP,Cpu>* conjugate();
			TransformHelper<FPP,Cpu>* adjoint();
			TransformHelper<FPP,Cpu>* vertcat(const TransformHelper<FPP,Cpu>*);
			TransformHelper<FPP,Cpu>* horzcat(const TransformHelper<FPP,Cpu>*);
			bool isTransposed() const;
			bool isConjugate() const;
			const char isTransposed2char() const;
			double normL1() const;
			double normFro() const;
			double normInf() const;
			TransformHelper<FPP,Cpu>* normalize(const int meth = 2/* 1 for 1-norm, 2 for 2-norm, MAX for inf-norm */) const;

			TransformHelper<FPP,Cpu>* pruneout(const int nnz_tres, const int npasses=-1, const bool only_forward=false);
			TransformHelper<FPP,Cpu>* optimize_storage(const bool time=true);
			TransformHelper<FPP,Cpu>* optimize(const bool transp=false);
			void optimize_multiply(const bool transp=false);
			/**
			  \brief Returns the left hand side factors of this from index 0 to id included (as a new TransformHelper obj).

			  */
			TransformHelper<FPP,Cpu>* left(const faust_unsigned_int id, const bool copy=false) const;
			/**
			  \brief Returns the right hand side factors of this from index id to the size()-1 (as a new TransformHelper obj).

			  */
			TransformHelper<FPP,Cpu>* right(const faust_unsigned_int id, const bool copy=false) const;

			transf_iterator<FPP> begin() const;
			transf_iterator<FPP> end() const;

			void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id);
			void pack_factors();
			void pack_factors(const faust_unsigned_int id, const PackDir dir);
            /** for testing purpose only (memory leaks enabled) */
			void disable_dtor() { this->transform->disable_dtor(); }
            void enable_dtor() { this->transform->enable_dtor(); }
			static TransformHelper<FPP,Cpu>* randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true);
			static TransformHelper<FPP,Cpu>* hadamardFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* fourierFaust(unsigned int n, const bool norma=true);
			static TransformHelper<FPP,Cpu>* eyeFaust(unsigned int n, unsigned int m);

			~TransformHelper();
			unsigned long long get_fact_addr(const faust_unsigned_int id) const;
			/* use carefully */
			MatGeneric<FPP,Cpu>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
			private:
			void copy_slices(TransformHelper<FPP, Cpu>* th, const bool transpose = false);
			const MatGeneric<FPP,Cpu>* get_gen_fact(const faust_unsigned_int id) const;
		};
}

#include "faust_TransformHelper.hpp"
#endif
