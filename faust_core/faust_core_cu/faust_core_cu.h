#ifndef __FAUST_CORE_CU_H__
#define __FAUST_CORE_CU_H__

#include <vector>
#include "faust_cu_spmat.h"
#include "cusparse.h"
#include "cublas_v2.h"

// class faust_cu_mat<T>;
// class faust_cu_vec;


template<typename T> class faust_cu_mat;
template<typename T> class faust_cu_spmat;
template<typename T> class faust_cu_vec;
template<typename T> class faust_core_cu;

template<typename T>
faust_cu_vec<T> operator*(const faust_core_cu<T>& f, const faust_cu_vec<T>& v);
template<typename T>
faust_cu_mat<T> operator*(const faust_core_cu<T>& f, const faust_cu_mat<T>& M);



template<typename T>
class faust_core_cu
{
	public:
		faust_core_cu();
		faust_core_cu(const std::vector<faust_cu_spmat<T> >& facts, const T lambda_ = (T)1.0);
		faust_core_cu(const faust_core_cu<T> & A);
		faust_core_cu(const std::vector<faust_cu_mat<T> >&facts);
		void get_facts(std::vector<faust_cu_spmat<T> >& sparse_facts)const{sparse_facts = data;}
		void get_facts(std::vector<faust_cu_mat<T> >& facts)const;	
		int size()const{return data.size();} 
      faust_cu_mat<T> get_product(cublasHandle_t, cusparseHandle_t)const;
		faust_cu_spmat<T> get_fact(int id) const;		
		int getNbRow() const;
		int getNbCol() const;
		void print_file(const char* filename) const;
		void init_from_file(const char* filename);
		long long int get_total_nnz()const{return totalNonZeros;}
		void clear(){data.resize(0);totalNonZeros=0;}
		void push_back(const faust_cu_spmat<T>& S);
		void mult( const faust_cu_vec<T>& cu_x, faust_cu_vec<T>& cu_y, cusparseHandle_t cusparseHandle);
		void push_first(const faust_cu_spmat<T>& S);
		void pop_back(faust_cu_spmat<T>& S);
		void pop_first(faust_cu_spmat<T>& S);
		void pop_first(faust_cu_spmat<T>& S) const;
		void Display()const;
		void transpose();
		void updateNonZeros();
		//(*this) = (*this) * A
		void multiply(const faust_core_cu<T> & A);
		//(*this) = A * (*this)
		void multiplyLeft(const faust_core_cu<T> & A);	
		void scalarMultiply(const T scalar);
		T spectralNorm(const int nbr_iter_max, T threshold, int &flag) const;
		~faust_core_cu(){}


	public:
		void operator=(const faust_core_cu<T>&  f){data=f.data;totalNonZeros=f.totalNonZeros;}
		// add all of the sparse matrices from f.data to this->data
		void operator*=(const T  scalar){scalarMultiply(scalar);};
		void operator*=(const faust_core_cu<T>&  f){multiply(f);};
		// add the sparse matrix S to this->data
		void operator*=(const faust_cu_spmat<T>&  S){push_back(S);totalNonZeros+=S.getNonZeros();}




	private:
		std::vector<faust_cu_spmat<T> > data;
		long long int totalNonZeros;
		static const char * class_name;
		

	friend faust_cu_vec<T> operator*<>(const faust_core_cu<T>& f, const faust_cu_vec<T>& v);
	friend faust_cu_mat<T> operator*<>(const faust_core_cu<T>& f, const faust_cu_mat<T>& M);
	// friend void multiply<>(const faust_core_cu<T> & A, const faust_cu_mat<T> & B, faust_cu_mat<T> & C,const T & alpha, char typeA, char typeMult);

		
};

#include "faust_core_cu.hpp"

#endif
