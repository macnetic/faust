#ifndef __FAUST_TRANSFORM_GPU2__
#define __FAUST_TRANSFORM_GPU2__
#ifdef USE_GPU_MOD
#include "faust_gpu_mod_utils.h"
#include "faust_constant.h"
#include "faust_MatGeneric_gpu.h"
#include "faust_MatSparse_gpu.h"
#include "faust_MatDense_gpu.h"
#include "faust_RefManager.h"
#include <vector>

namespace Faust
{

	//TODO: a TransformGen class to refactor common code between Transform<FPP,Cpu> and Transform<FPP,GPU2> (similar way to what is done for TransformHelper classes)
	template<typename FPP>
		class Transform<FPP,GPU2>
		{
			gm_MatArray_t gpu_mat_arr;
			std::vector<MatGeneric<FPP,GPU2>*> data; // same factors as in gpu_mat_arr
													 // just pointer copies
			bool dtor_delete_data;
			bool dtor_disabled;
			static RefManager ref_man;
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
			void transpose();
			int32_t getNbRow()const;
			int32_t getNbCol()const;
			void Display(bool transpose=false) const;
			std::string to_string(bool transpose=false) const;
			int32_t size() const;
			faust_unsigned_int get_fact_nnz(const faust_unsigned_int id) const;
			faust_unsigned_int get_total_nnz() const;
			void update_total_nnz() const;
			void scalarMultiply(const FPP& alpha);
			MatDense<FPP,GPU2> get_product(const char opThis='N', const bool isConj=false) const;
			void get_product(MatDense<FPP,GPU2>& M, const char opThis='N', const bool isConj=false) const;
			MatDense<FPP,GPU2> multiply(const MatDense<FPP,GPU2> &A, const char opThis);
			void multiply(const Transform<FPP,GPU2> & A);
			void multiplyLeft(const Transform<FPP,GPU2> & A);
			void multiply(const FPP& a, const int32_t id=-1);
			Vect<FPP,GPU2> multiply(const Vect<FPP,GPU2>& x, const char opThis='N');
			Real<FPP> spectralNorm(int32_t nb_iter_max, float threshold, int& flag);
			FPP power_iteration(int32_t nb_iter_max, float threshold, int& flag);
			Real<FPP> normL1(const bool transpose = false) const;
			void tocpu(Transform<FPP, Cpu>& cpu_transf) const;
			Transform<FPP, Cpu> tocpu() const;
			void save_mat_file(const char* filename, const bool transpose, const bool conjugate) const;
			vector<MatGeneric<FPP,GPU2>*> getData() const {return data;} //TODO: delete and make TransformHelper a friend

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
		};
}

#endif
#endif
