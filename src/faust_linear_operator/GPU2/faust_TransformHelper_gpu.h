#ifndef __TRANSFORM_HELPER_GPU2__
#define __TRANSFORM_HELPER_GPU2__
#include "faust_constant.h"
#include "faust_Transform_gpu.h"
//#include "faust_Transform.h"
#include "faust_TransformHelperGen.h"
#include <memory>
namespace Faust
{
	template<typename FPP, FDevice DEVICE> class TransformHelper;
	template<typename FPP, FDevice DEVICE> class TransformHelperGen;
	template<typename FPP>
		class TransformHelper<FPP,GPU2> : public TransformHelperGen<FPP,GPU2>
		{
			public:
				TransformHelper();
				TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact = true, const bool internal_call=false);
#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
				template<typename ...GList> TransformHelper(GList& ... t);
#endif
				void push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy=false, const bool copying=true);
				void push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy=false, const int32_t dev_id=-1, const void* stream=nullptr);
				void pop_front();
				void pop_back();
				void push_first(const MatGeneric<FPP,GPU2>*, const bool optimizedCopy=false, const bool copying=true);
				faust_unsigned_int getNbRow(){ return this->transform->getNbRow();}
				faust_unsigned_int getNbCol(){ return this->transform->getNbCol();}
				unsigned int get_fact_nb_rows(const faust_unsigned_int id) const;
				unsigned int get_fact_nb_cols(const faust_unsigned_int id) const;
				template<typename Head, typename ... Tail>
					void push_back_(Head& h, Tail&... t);
				void push_back_();
				void display() const;
				std::string to_string() const;
				MatDense<FPP,GPU2> get_product();
				void get_product(MatDense<FPP,GPU2>& M);
				MatDense<FPP,GPU2> multiply(const MatDense<FPP,GPU2> &A, const bool transpose=false, const bool conjugate=false);
				TransformHelper<FPP,GPU2>* multiply(const FPP& a);
				Vect<FPP,GPU2> multiply(const Faust::Vect<FPP,GPU2>& a);
				Real<FPP> normFro() const;
				Real<FPP> normL1() const;
				Real<FPP> normInf() const;
				faust_unsigned_int size() const;
				void update_total_nnz() const;
				Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
				bool is_fact_sparse(int id) const;
				bool is_fact_dense(int id) const;
				MatGeneric<FPP,GPU2>* get_gen_fact_nonconst(const faust_unsigned_int id) const;
				void pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id);
				void pack_factors();
				void update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id);
				void operator=(TransformHelper<FPP,GPU2>& th);
				typename Transform<FPP,GPU2>::iterator begin() const;
				typename Transform<FPP,GPU2>::iterator end() const;
				void tocpu(TransformHelper<FPP,Cpu>& cpu_transf);
				void save_mat_file(const char* filename) const;
				TransformHelper<FPP,GPU2>* vertcat(const TransformHelper<FPP,GPU2>*);
				TransformHelper<FPP,GPU2>* horzcat(const TransformHelper<FPP,GPU2>*);
				void set_FM_mul_mode() const;
				void set_Fv_mul_mode() const;
				faust_unsigned_int get_total_nnz() const;
				faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
				bool is_fact_sparse(const faust_unsigned_int id) const;
				TransformHelper<FPP,GPU2>* normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const;
				TransformHelper<FPP,GPU2>* transpose();
				TransformHelper<FPP,GPU2>* conjugate();
				TransformHelper<FPP,GPU2>* adjoint();
		};
}
#include "faust_TransformHelper_gpu.hpp"
#endif
