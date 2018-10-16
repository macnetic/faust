
#ifndef __FAUST_TRANSFORM_HELPER___
#define __FAUST_TRANSFORM_HELPER___

#include <memory>
#include "faust_exception.h"
#include "faust_Transform.h"
#include "faust_Slice.h"
#include <random>

namespace Faust {
	using namespace std;

	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class TransformHelper;

	enum RandFaustType {
		DENSE,
		SPARSE,
		MIXTE
	};

	template<typename FPP>
		class TransformHelper<FPP,Cpu> {
			static std::default_random_engine generator;
			static bool seed_init;

			bool is_transposed;
			bool is_conjugate;
			bool is_sliced;
			Slice slices[2];
			Transform<FPP,Cpu>* sliced_transform;
			shared_ptr<Transform<FPP,Cpu>> transform;

			const Transform<FPP, Cpu>* eval_sliced_Transform();
			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=true, const bool cloning_fact = true);
			TransformHelper();
			TransformHelper(TransformHelper<FPP,Cpu>* th_left, TransformHelper<FPP,Cpu>* th_right);
			TransformHelper(TransformHelper<FPP,Cpu>* th);
			TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
			TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]);
			TransformHelper(Transform<FPP,Cpu> &t);

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x) const;
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x, const bool transpose);
			MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> A) const;
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> A, const bool transpose);
			TransformHelper<FPP, Cpu>* multiply(TransformHelper<FPP, Cpu>*);
			TransformHelper<FPP, Cpu>* multiply(FPP& scalar);
			void push_back(const MatGeneric<FPP,Cpu>* M);
			faust_unsigned_int getNbRow() const;
			faust_unsigned_int getNbCol() const;
			bool isReal() const;
			faust_unsigned_int get_total_nnz() const;
			faust_unsigned_int size() const;
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
			MatDense<FPP,Cpu> get_product() const;
			void save_mat_file(const char* filename) const;
			double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
			TransformHelper<FPP,Cpu>* transpose();
			TransformHelper<FPP,Cpu>* conjugate();
			TransformHelper<FPP,Cpu>* adjoint();
			bool isTransposed() const;
			bool isConjugate() const;
			const char isTransposed2char() const;
			double normL1() const;
			double normFro() const;
			static TransformHelper<FPP,Cpu>* randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f);
			static TransformHelper<FPP,Cpu>* hadamardFaust(unsigned int n);
			~TransformHelper();

			private:
			void copy_slices(TransformHelper<FPP, Cpu>* th, const bool transpose = false);
		};
}

#include "faust_TransformHelper.hpp"
#endif
