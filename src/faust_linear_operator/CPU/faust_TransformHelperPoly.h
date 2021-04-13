#ifndef __FAUST_TRANSFORM_HELPER_POLY__
#define __FAUST_TRANSFORM_HELPER_POLY__
#include "faust_TransformHelper.h"

#define ERROR_ON_FAC_NUM_CHANGE() \
	throw std::runtime_error("Can't pack factors of a TransformHelperPoly (the number of factors must be kept consistent with the basis size).");

typedef unsigned int uint;

namespace Faust
{
	enum BasisLaziness
	{
		NOT_LAZY, // all the factors are instantiated at initialization time
		INSTANTIATE_ONCE_AND_FOR_ALL, // a factor instantiated is kept until the "Faust" is deleted
		INSTANTIATE_COMPUTE_AND_FREE // a factor is instantiated and used locally then is freed
	};

	template<typename FPP>
		TransformHelper<FPP, Cpu>* basisChebyshev(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0=nullptr, bool on_gpu=false, BasisLaziness lazy_instantiation=INSTANTIATE_ONCE_AND_FOR_ALL);

	template<typename FPP>
		void poly(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out, bool on_gpu=false);
	template<typename FPP>
		void poly_cpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out);
	template<typename FPP>
		void poly_gpu(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP* out);


	/**
	 * \brief This class aims to represent a Chebyshev polynomial basis as a "Faust".
	 *
	 * It overrides some operations to optimize the performance by taking into account the factors specific structure.
	 */
	template<typename FPP>
		class TransformHelperPoly : public TransformHelper<FPP,Cpu>
		{

			static RefManager ref_man;
			MatSparse<FPP, Cpu> *L;
			MatSparse<FPP, Cpu> *rR;
			std::vector<bool> is_fact_created;
			BasisLaziness laziness;
			bool T0_is_arbitrary;
			bool mul_and_combi_lin_on_gpu;

			// ctor is private on purpose (callees must call basisChebyshev() to instantiate an object
			TransformHelperPoly(uint K,
					MatSparse<FPP, Cpu> *L,
					MatSparse<FPP, Cpu> *rR=nullptr,
					MatSparse<FPP, Cpu> *T0=nullptr,
					BasisLaziness laziness=INSTANTIATE_COMPUTE_AND_FREE,
					bool on_gpu=false);

			TransformHelperPoly(uint K, const TransformHelperPoly<FPP>& src);

			public:
			faust_unsigned_int getNbRow() const;
			faust_unsigned_int getNbCol() const;
			void get_fact(const faust_unsigned_int &id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			void get_fact(const faust_unsigned_int id,
					int* rowptr,
					int* col_ids,
					FPP* elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose = false) const;
			uint get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const;
			const MatGeneric<FPP,Cpu>* get_gen_fact(const faust_unsigned_int id) const;
			faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			bool is_fact_sparse(const faust_unsigned_int id) const;
			bool is_fact_dense(const faust_unsigned_int id) const;
			void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode=DEFAULT);

			void pack_factors(const faust_unsigned_int id, const PackDir dir, const int mul_order_opt_mode=DEFAULT);
			void pack_factors(const int mul_order_opt_mode=DEFAULT);
			TransformHelper<FPP,Cpu>* left(const faust_unsigned_int id, const bool copy=false) const;
			TransformHelper<FPP,Cpu>* right(const faust_unsigned_int id, const bool copy=false) const;
			TransformHelper<FPP,Cpu>* optimize_storage(const bool time=true);
			TransformHelper<FPP,Cpu>* clone();

			void update_total_nnz();
			MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A, const bool transpose=false, const bool conjugate=false);
			TransformHelper<FPP, Cpu>* multiply(const TransformHelper<FPP, Cpu>*) const;
			TransformHelper<FPP, Cpu>* multiply(FPP& scalar);
			faust_unsigned_int getNBytes() const;
			faust_unsigned_int get_total_nnz() const;
			void resize(faust_unsigned_int);
			string to_string() const;

			MatDense<FPP,Cpu> get_product(const int mul_order_opt_mode=DEFAULT);// const;
			void get_product(MatDense<FPP,Cpu>& prod, const int mul_order_opt_mode=DEFAULT); //const;
			void save_mat_file(const char* filename) const;
			double spectralNorm(const int nbr_iter_max, double threshold, int &flag) const;
			TransformHelper<FPP,Cpu>* vertcat(const TransformHelper<FPP,Cpu>*);
			TransformHelper<FPP,Cpu>* horzcat(const TransformHelper<FPP,Cpu>*);
			double normL1() const;
			double normFro() const;
			double normInf() const;
			TransformHelper<FPP,Cpu>* normalize(const int meth = 2/* 1 for 1-norm, 2 for 2-norm, MAX for inf-norm */) const;
			TransformHelper<FPP,Cpu>* pruneout(const int nnz_tres, const int npasses=-1, const bool only_forward=false);
			TransformHelper<FPP,Cpu>* optimize_multiply(std::function<void()> f, const bool transp=false, const bool inplace=false, const int nsamples=1, const char* op_name="unamed_op");
			TransformHelper<FPP,Cpu>* optimize_time(const bool transp=false, const bool inplace=false, const int nsamples=1);
			TransformHelper<FPP,Cpu>* optimize_time_full(const bool transp=false, const bool inplace=false, const int nsamples=1);
			TransformHelper<FPP,Cpu>* optimize_time_Fv(const bool transp=false, const bool inplace=false, const int nsamples=1);
			TransformHelper<FPP,Cpu>* swap_cols(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);
			TransformHelper<FPP,Cpu>* swap_rows(const faust_unsigned_int id1, const faust_unsigned_int id2, const bool permutation=false, const bool inplace=false, const bool check_transpose=true);
            TransformHelper<FPP, Cpu>* slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
					                    faust_unsigned_int start_col_id, faust_unsigned_int end_col_id);
			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &x, const bool transpose=false, const bool conjugate=false);
			Vect<FPP,Cpu> multiply(const FPP* x, const bool transpose=false, const bool conjugate=false);
			void multiply(const FPP* x, FPP* y, const bool transpose=false, const bool conjugate=false);
			void multiply_cpu(const FPP* x, FPP* y, const bool transpose=false, const bool conjugate=false);
			void multiply_gpu(const FPP* x, FPP* y, const bool transpose=false, const bool conjugate=false);
			MatDense<FPP, Cpu> multiply(const MatDense<FPP,Cpu> &X, const bool transpose=false, const bool conjugate=false);
			void multiply(const FPP* X, int n, FPP* out, const bool transpose=false, const bool conjugate=false);
			void multiply_cpu(const FPP* X, int n, FPP* out, const bool transpose=false, const bool conjugate=false);
			void multiply_gpu(const FPP* X, int n, FPP* out, const bool transpose=false, const bool conjugate=false);
			TransformHelper<FPP, Cpu>* next(uint K);
			TransformHelper<FPP, Cpu>* next();
			Vect<FPP, Cpu> poly(MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			MatDense<FPP, Cpu> poly(int n, MatDense<FPP,Cpu> & basisX, Vect<FPP, Cpu> coeffs);
			void poly(int d, int n, const FPP* basisX, const FPP* coeffs, FPP* out);
			TransformHelper<FPP, Cpu>* polyFaust(const FPP* coeffs);
			/**
			 * Computes the multiplication of "this" by x and the linear combination of polynomials at the same time (in a memory saving purpose).
			 *
			 * \param x: a vector of size L->getNbRow()
			 * \param y: output vector of size L->getNbRow()
			 * \param coeffs: a vector of size K+1 (this->size()), coefficients of the linear combination.
			 */
			void multiplyPoly(const FPP* x, FPP* y, const FPP* coeffs);
			void multiplyPoly_gpu(const FPP* x, FPP* y, const FPP* coeffs);
			void multiplyPoly_cpu(const FPP* x, FPP* y, const FPP* coeffs);
			/**
			 * Computes the multiplication of "this" by X and the linear combination of polynomials at the same time (in a memory saving purpose).
			 *
			 * \param X: a matrix of size L->getNbRow() x n
			 * \param Y: output matrix of size L->getNbRow() x n
			 * \param coeffs: a vector of size K+1 (this->size()), coefficients of the linear combination.
			 */
			void multiplyPoly(const FPP* X, int n, FPP* Y, const FPP* coeffs);
			void multiplyPoly_cpu(const FPP* X, int n, FPP* Y, const FPP* coeffs);
			void multiplyPoly_gpu(const FPP* X, int n, FPP* Y, const FPP* coeffs);
			void basisChebyshevT0(MatSparse<FPP,Cpu>* T0=nullptr);
			void basisChebyshevT1();
			void basisChebyshevT2();
			void basisChebyshevTi(uint i);
			void basisChebyshev_facti(uint i);
			void basisChebyshev_free_facti(uint i);
			void basisChebyshev_facti2j(uint i, uint j);
			void basisChebyshev_free_facti2j(uint i, uint j);
			void basisChebyshev_fact_all();
			void basisChebyshev_free_fact_all();

			void basisChebyshev_all();
			void create_rR(const MatSparse<FPP,Cpu>* L);

			~TransformHelperPoly();
			friend TransformHelper<FPP,Cpu>* basisChebyshev<>(MatSparse<FPP,Cpu>* L, int32_t K, MatSparse<FPP, Cpu>* T0, bool on_gpu, BasisLaziness lazy_instantiation);
		};


}
#include "faust_TransformHelperPoly.hpp"
#endif

