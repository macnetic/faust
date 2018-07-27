
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

			shared_ptr<Transform<FPP,Cpu>> transform;


			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=true, const bool cloning_fact = true);
			TransformHelper();
			TransformHelper(TransformHelper<FPP,Cpu>* th);
			TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
			TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2]);

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x) const;
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> x, const bool transpose);
			MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> A) const;
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> A, const bool transpose);
			void push_back(const MatGeneric<FPP,Cpu>* M);
			faust_unsigned_int getNbRow() const;
			faust_unsigned_int getNbCol() const;
			bool isReal() const;
			faust_unsigned_int get_total_nnz() const;
			faust_unsigned_int size() const;
			void display() const;
			string to_string() const;
			MatDense<FPP,Cpu> get_fact(faust_unsigned_int id) const;
			const Transform<FPP, Cpu>* eval_sliced_Transform() const;
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
			~TransformHelper();

			private:
			void copy_slices(TransformHelper<FPP, Cpu>* th, const bool transpose = false);
		};
}

#include "faust_TransformHelper.hpp"
#endif
