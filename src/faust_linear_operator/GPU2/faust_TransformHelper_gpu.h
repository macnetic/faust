#ifndef __TRANSFORM_HELPER_GPU2__
#define __TRANSFORM_HELPER_GPU2__
#include "faust_constant.h"
#include "faust_Transform_gpu.h"
//#include "faust_Transform.h"
#include "faust_TransformHelperGen.h"
#include <memory>
#include "faust_TransformHelper_cat_gpu.h"
namespace Faust
{
	template<typename FPP, FDevice DEVICE> class TransformHelper;
	template<typename FPP, FDevice DEVICE> class TransformHelperGen;
	template<typename FPP>
		class TransformHelper<FPP,GPU2> : public TransformHelperGen<FPP,GPU2>
		{
			public:
				TransformHelper();
				TransformHelper(TransformHelper<FPP,GPU2>* th, Slice s[2]);
				TransformHelper(TransformHelper<FPP,GPU2>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
				TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
				TransformHelper(const TransformHelper<FPP,Cpu>& cpu_t, int32_t dev_id=-1, void* stream=nullptr);
				TransformHelper(const TransformHelper<FPP,GPU2>& th, bool transpose, bool conjugate);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
				template<typename ...GList> TransformHelper(GList& ... t);
#endif

				void push_back(const FPP* data, int nrows, int ncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);
				void push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);
				void push_back(const FPP* bdata, const int* brow_ptr, const int* bcol_inds, const int nrows, const int ncols, const int bnnz, const int bnrows, const int bncols, const bool optimizedCopy=false, const bool transpose=false, const bool conjugate=false);
				void push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy=false, const bool copying=true, const bool transpose=false, const bool conjugate=false);
				void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const int32_t dev_id=-1, const void* stream=nullptr);
				void pop_front();
				void pop_back();
				void push_first(const MatGeneric<FPP,GPU2>*, const bool optimizedCopy=false, const bool copying=true);
				faust_unsigned_int getNBytes() const;
				template<typename Head, typename ... Tail>
					void push_back_(Head& h, Tail&... t);
				void push_back_();
				void display() const;
				std::string to_string() const;
				MatDense<FPP,GPU2> get_product(int prod_mod=-1);
				void get_product(MatDense<FPP,GPU2>& M, int prod_mod=-1);
				void get_product(MatDense<FPP,Cpu>& M, int prod_mod=-1);
				virtual MatDense<FPP, GPU2> multiply(const MatDense<FPP,GPU2> &A);
				virtual MatDense<FPP, GPU2> multiply(const MatSparse<FPP,GPU2> &A) { return multiply(MatDense<FPP,GPU2>(A));}
				virtual MatDense<FPP, Cpu> multiply(const MatSparse<FPP,Cpu> &A) { return multiply(MatDense<FPP, Cpu>(A));} // TODO: avoid CPU copy to dense
				virtual MatDense<FPP,Cpu> multiply(const MatDense<FPP,Cpu> &A);
				TransformHelper<FPP,GPU2>* multiply(const FPP& a);
				TransformHelper<FPP,GPU2>* multiply(const TransformHelper<FPP,GPU2>*);
				virtual Vect<FPP,GPU2> multiply(const Faust::Vect<FPP,GPU2>& a);
				virtual Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> &x);
				virtual void multiply(const FPP* cpu_x, FPP* cpu_y);
				virtual void multiply(const FPP* cpu_x, int x_ncols, FPP* cpu_y);
				FPP* sliceMultiply(const Slice s[2], const FPP* cpu_X, FPP* cpu_out=nullptr, int X_ncols=1) const;
				FPP* indexMultiply(faust_unsigned_int* ids[2], size_t id_lens[2], const FPP* X, int ncols=1, FPP* out=nullptr) const;
				Real<FPP> normFro(const bool full_array=true, const int batch_size=1) const;
				Real<FPP> normL1(const bool full_array=true, const int batch_size=1) const;
				Real<FPP> normInf(const bool full_array=true, const int batch_size=1) const;
				faust_unsigned_int size() const;
				void update_total_nnz() const;
				Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
				FPP power_iteration(const faust_unsigned_int nbr_iter_max, const Real<FPP>& threshold, int & flag);
				const MatGeneric<FPP,GPU2>* get_gen_fact(const faust_unsigned_int id) const;
				MatGeneric<FPP,GPU2>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
				void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id, const int mul_order_opt_mode=DEFAULT_L2R);
				void update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id);
				void replace(const MatGeneric<FPP, GPU2>* M, const faust_unsigned_int id);
				void convertToSparse();
				void convertToDense();
				template<typename FPP2>
				TransformHelper<FPP2, GPU2>* cast();


				void operator=(TransformHelper<FPP,GPU2>& th);
				typename Transform<FPP,GPU2>::iterator begin() const;
				typename Transform<FPP,GPU2>::iterator end() const;
				void tocpu(TransformHelper<FPP,Cpu>& cpu_transf) const;
				TransformHelper<FPP,Cpu>* tocpu() const;
				TransformHelper<FPP,GPU2>* vertcat(const TransformHelper<FPP,GPU2>*);
				TransformHelper<FPP,GPU2>* horzcat(const TransformHelper<FPP,GPU2>*);
				TransformHelper<FPP,GPU2>* swap_cols(const faust_unsigned_int id1,
						const faust_unsigned_int id2,
						const bool permutation=false,
						const bool inplace=false,
						const bool check_transpose=true);
				TransformHelper<FPP,GPU2>* swap_rows(const faust_unsigned_int id1,
						const faust_unsigned_int id2,
						const bool permutation=false,
						const bool inplace=false,
						const bool check_transpose=true);
				faust_unsigned_int get_total_nnz() const;
				//				faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
				TransformHelper<FPP,GPU2>* normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const;
				TransformHelper<FPP,GPU2>* transpose();
				TransformHelper<FPP,GPU2>* conjugate();
				TransformHelper<FPP,GPU2>* adjoint();
				TransformHelper<FPP,GPU2>* optimize_multiply(std::function<void()> f, const bool transp=false, const bool inplace=false, const int nsamples=1, const char* op_name="unamed_op");
				TransformHelper<FPP,GPU2>* clone(int32_t dev_id=-1, void* stream=nullptr);
				void get_fact(const faust_unsigned_int id,
						int* rowptr,
						int* col_ids,
						FPP* elts,
						faust_unsigned_int* nnz,
						faust_unsigned_int* num_rows,
						faust_unsigned_int* num_cols,
						const bool transpose /* = false*/) const;
				void get_fact(const faust_unsigned_int id,
						FPP* elts,
						faust_unsigned_int* num_rows,
						faust_unsigned_int* num_cols,
						const bool transpose  = false) const;

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



				void save_mat_file(const char* filepath) const;
				TransformHelper<FPP,GPU2> read_from_mat_file(const char* filepath);
				static int get_mat_file_type(const char* filepath);
				FPP get_item(faust_unsigned_int i, faust_unsigned_int j);
				void init_fancy_idx_transform(TransformHelper<FPP,GPU2>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols);
				virtual ~TransformHelper() {};
				static TransformHelper<FPP,GPU2>* hadamardFaust(unsigned int n, const bool norma=true);
				static TransformHelper<FPP,GPU2>* fourierFaust(unsigned int n, const bool norma=true);
				static TransformHelper<FPP,GPU2>* fourierFaustOpt(unsigned int n, const bool norma=true);
				static TransformHelper<FPP,GPU2>* optButterflyFaust(const TransformHelper<FPP, GPU2>* F);
				static TransformHelper<FPP,GPU2>* eyeFaust(unsigned int n, unsigned int m);
				TransformHelper<FPP,GPU2>* pruneout(const int nnz_tres, const int npasses=-1, const bool only_forward=false);

				static TransformHelper<FPP,GPU2>* randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true, unsigned int seed=0);
				static TransformHelper<FPP,GPU2>* randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density=.1f, bool per_row=true, unsigned int seed=0);

				static TransformHelper<FPP, GPU2>* randBSRFaust(unsigned int faust_nrows, unsigned int faust_ncols, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int bnrows, unsigned int bncols, float density=.1f);

		};


}
#include "faust_TransformHelper_gpu.hpp"
#endif
