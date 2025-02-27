#ifndef __FAUST_TRANSFORM_GPU2__
#define __FAUST_TRANSFORM_GPU2__
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatSparse_gpu.h"
#include "faust_MatDense_gpu.h"
#include "faust_RefManager.h"
#include "faust_Slice.h"
#include "faust_linear_algebra_gpu.h"
#include <vector>

namespace Faust
{

	//TODO: a TransformGen class to refactor common code between Transform<FPP,Cpu> and Transform<FPP,GPU2> (similar way to what is done for TransformHelper classes)
	template<typename FPP>
		class Transform<FPP,GPU2>
		{
			std::vector<MatGeneric<FPP,GPU2>*> data; 
			bool dtor_delete_data;
			bool dtor_disabled;
			static RefManager ref_man;
			faust_unsigned_int total_nnz;
			public:
			Transform();
			Transform(const std::vector<MatGeneric<FPP,GPU2>*> &factors, const FPP lambda_ = (FPP)1.0, const bool optimizedCopy=false, const bool cloning_fact=true);
			Transform(const Transform<FPP,GPU2>& t);
			~Transform();
			void operator=(const Transform<FPP,GPU2>& t);
			void insert(int32_t id, const MatGeneric<FPP,GPU2>*, bool copying=true);
			void push_back(const MatGeneric<FPP,GPU2>*, bool copying=true, const bool transpose=false, const bool conjugate=false);
			void push_first(const MatGeneric<FPP,GPU2>*, bool copying=true);
			void erase(int32_t id);
			void pop_front();
			void pop_back();
			void clear();
			void update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id);
			void replace(const MatGeneric<FPP, GPU2>* M, const faust_unsigned_int id);
			MatGeneric<FPP,GPU2>* get_fact(int32_t id, bool cloning_fact=true) const;
			void get_fact(const faust_unsigned_int &id,
					FPP* elts,
					faust_unsigned_int* num_rows,
					faust_unsigned_int* num_cols,
					const bool transpose=false) const;
			void get_fact(const faust_unsigned_int id,
					int* d_outer_count_ptr, int* d_inner_ptr, FPP* d_elts,
					faust_unsigned_int* nnz,
					faust_unsigned_int* num_rows, faust_unsigned_int* num_cols,
					bool transpose=false) const;
			void get_facts(std::vector<MatGeneric<FPP,GPU2>*> &factors, bool cloning_facts=true) const;
			bool is_fact_sparse(int id) const;
			bool is_fact_dense(int id) const;
			bool is_fact_bsr(int id) const;
			void transpose();
			void adjoint();
			faust_unsigned_int getNbRow()const;
			faust_unsigned_int getNbCol()const;
			float getRCG() const{return ((float)(getNbRow()*getNbCol()))/((float) get_total_nnz());} //TODO: move in hpp and factorize with CPU code
			void Display(const bool transpose=false, const bool display_small_mat_elts=false) const;
			std::string to_string(const bool transpose=false, const bool display_small_mat_elts=false) const;
			faust_unsigned_int size() const;
			faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			faust_unsigned_int get_total_nnz() const;
			void update_total_nnz();
			void scalarMultiply(const FPP& alpha, long int sid=-1);
			MatDense<FPP,GPU2> get_product(const char opThis='N', const bool isConj=false) const;
			void get_product(MatDense<FPP,GPU2>& M, const char opThis='N', const bool isConj=false) const;
			MatDense<FPP,GPU2> multiply(const MatDense<FPP,GPU2> &A, const char opThis);
			void multiply(const Transform<FPP,GPU2> & A);
			MatDense<FPP, GPU2> sliceMultiply(const Slice s[2], MatDense<FPP, GPU2>& gpu_X, const char opThis) const;
			MatDense<FPP, GPU2> indexMultiply(faust_unsigned_int* ids[2], size_t id_lens[2], MatDense<FPP, GPU2>& gpu_X, const char opThis) const;
			void multiplyLeft(const Transform<FPP,GPU2> & A);
			Vect<FPP,GPU2> multiply(const Vect<FPP,GPU2>& x, const char opThis='N');
			Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag) const;
			template<typename FPP2 = double>
				FPP power_iteration(const faust_unsigned_int nbr_iter_max, const FPP2 threshold, int & flag, const bool rand_init=true) const;
			Real<FPP> normL1(const bool transpose = false, const bool full_array=true, const int batch_size=1) const;
			void tocpu(Transform<FPP, Cpu>& cpu_transf) const;
			Transform<FPP, Cpu> tocpu() const;
			void save_mat_file(const char* filename, const bool transpose, const bool conjugate) const;
			vector<MatGeneric<FPP,GPU2>*> getData() const {return data;} //TODO: delete and make TransformHelper a friend
			gm_MatArray_t asGMObj() const;

//			using transf_iterator = typename std::vector<MatGeneric<FPP,Cpu>*>::const_iterator;
//
//			transf_iterator begin() const;
//
//			transf_iterator end() const;
			public:
			class iterator : public std::iterator<std::output_iterator_tag, MatGeneric<FPP,GPU2>*>
			{
				public:
					explicit iterator(const Transform<FPP, GPU2>& container, size_t index = 0);
					MatGeneric<FPP,GPU2>* operator*() const;
					iterator & operator++();
					iterator operator+(int);
					iterator operator-(int);
					//post-increment op
					iterator operator++(int);
					bool operator!=(const iterator& it);
					bool operator<(const iterator& it);
				private:
					size_t index;
					const Transform<FPP, GPU2> & container;
			};

			typename Transform<FPP,GPU2>::iterator begin() const;
			typename Transform<FPP,GPU2>::iterator end() const;

			friend TransformHelper<FPP,GPU2>;
			friend TransformHelperGen<FPP,GPU2>;
		};

}

#include "faust_Transform_gpu.hpp"
#endif
#endif
