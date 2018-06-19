
#ifndef __FAUST_TRANSFORM_HELPER___
#define __FAUST_TRANSFORM_HELPER___

#include <memory>
#include "faust_exception.h"
#include "faust_Transform.h"

namespace Faust {
	using namespace std;

	template<typename FPP,Device DEVICE> class Transform;
	template<typename FPP,Device DEVICE> class Vect;
	template<typename FPP,Device DEVICE> class MatDense;
	template<typename FPP,Device DEVICE> class TransformHelper;

	template<typename FPP>
		class TransformHelper<FPP,Cpu> {

			bool is_transposed;
			bool is_conjugate;

			shared_ptr<Transform<FPP,Cpu>> transform;

			public:
			TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=true);
			TransformHelper();
			TransformHelper(TransformHelper<FPP,Cpu>* th);
			TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate);
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
			MatDense<FPP,Cpu> get_fact(faust_unsigned_int id) const;
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
			~TransformHelper();
		};
}
#include "faust_TransformHelper.hpp"

#endif
