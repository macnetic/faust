#ifndef __FAUST_SPMAT_H__
#define __FAUST_SPMAT_H__

#include "faust_constant.h"
#include "faust_mat.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "faust_mat_generic.h"
#include "faust_vec.h"

//class faust_vec<T>;
template<typename T> class faust_mat;
template<typename T> class faust_spmat;
template<typename T> class faust_vec;

template<typename T>
void multiply(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);

/*! \class faust_spmat
   * \brief template class representing sparse matrix
   *
   *  This class implements sparse matrix multiplication
   *
   * The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
   *
   *\tparam T scalar numeric type, e.g float or double
   */


template<typename T>
class faust_spmat : public faust_mat_generic
{
	public:
		faust_spmat();
		template<typename U>
		faust_spmat(const faust_spmat<U>& M){(this)->operator=(M);}
		template<typename U>
		faust_spmat(const faust_mat<U>& M){(this)->operator=(M);}
		faust_spmat(const faust_spmat<T>& M);
		faust_spmat(const faust_mat<T>& M);
		faust_spmat(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		faust_spmat(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const T* value, const size_t* id_row, const size_t* col_ptr);
		void set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const double* value, const size_t* id_row, const size_t* col_ptr);
		void resize(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		void resize(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_){mat.resize(dim1_,dim2_);update_dim();}
		void setZeros(){mat.setZero();nnz=0;}
		void setEyes(){mat.setIdentity();update_dim();}
		void transpose();
		T norm() const {return mat.norm();}
		void operator= (const faust_spmat<T>& M);
		void operator= (const faust_mat<T>& Mdense);
		void operator*=(const T alpha);
		void operator/=(const T alpha);
		template <typename U>
		void operator=(const faust_spmat<U> &M);
		template <typename U>
		void operator=(const faust_mat<U> &M){faust_spmat<U> spM(M);(this)->operator=(spM);}

      int getNbRow()const{return dim1;}
		int getNbCol()const{return dim2;}
		faust_unsigned_int getNonZeros()const{return nnz;}
		bool isCompressedMode()const{return mat.isCompressed();}
		void makeCompression(){mat.makeCompressed();}


		T* getValuePtr(){return mat.valuePtr();}// return pointer value of length nnz
		const T* getValuePtr()const{return mat.valuePtr();}// return const pointer value of length nnz



      // if rowMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim1,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the row (n-1)
      // if colMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim2,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the col (n-1)
		int* getOuterIndexPtr(){return mat.outerIndexPtr();}//return row-index value of length equal to the number of row+1
		const int* getOuterIndexPtr()const{return mat.outerIndexPtr();}

      // if rowMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] = column index matching the element getValuePtr()[n];
      // if colMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] =   row  index matching the element getValuePtr()[n];
		int* getInnerIndexPtr(){return mat.innerIndexPtr();}//return col index of length nnz
		const int* getInnerIndexPtr()const{return mat.innerIndexPtr();}
		
      const int* getRowPtr()const{if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(class_name,"getRowPtr : matrix is not in rowMajor");}}
      const int* getColInd()const{if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(class_name,"getColInd : matrix is not in rowMajor");}}
      int* getRowPtr(){if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(class_name,"getRowPtr : matrix is not in rowMajor");}}
      int* getColInd(){if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(class_name,"getColInd : matrix is not in rowMajor");}}

		
		

		void Display() const;
		T norm(){return mat.norm();}	
		/*!
     *  \brief 
	 * write faust_spmat into text file 
	 *\tparam filename : name of the file
	 * 
	 * the first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients
	 * all the other line contains 2 integers an 1 number :  the row, the column and the value of one coefficient 
	 * in ColMajor access of the faust_mat
     */
		void print_file(const char* filename)const;
		void print_file(const char* filename, std::ios_base::openmode mode)const;
		void init_from_file(FILE* fp);
	/*!
     *  \brief 
	 * init faust_spmat into text file 
	 *\tparam filename : name of the file
	 * 
	 * the first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients
	 * all the other line contains 2 integers an 1 number : the row, the column and the value of one coefficient 
	 * in ColMajor access of the faust_mat
     */	
		void init_from_file(const char* filename);
		void init(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		~faust_spmat(){}

	private:
		void update_dim(){dim1=mat.rows();dim2=mat.cols();nnz=mat.nonZeros();}
		static const char * class_name;

		

	private:
		Eigen::SparseMatrix<T,Eigen::RowMajor> mat;	
		faust_unsigned_int nnz;


	friend void faust_mat<T>::operator=(faust_spmat<T> const& S);

	// *this = (*this) * S
	friend void faust_mat<T>::operator*=(const faust_spmat<T>& S);
	// *this = (*this) + S
	friend void faust_mat<T>::operator+=(const faust_spmat<T>& S);
	// *this = (*this) - S
	friend void faust_mat<T>::operator-=(const faust_spmat<T>& S);
	// friend void sp_solve<>(const faust_spmat<T> & A,faust_vec<T> & x, const faust_vec<T> & y);

	// *this = S * (*this) 
	friend void faust_mat<T>::multiplyLeft(const faust_spmat<T>& S);

	// *this = S * (*this) 
	friend void  faust_vec<T>::multiplyLeft(faust_spmat<T> const& A);
	
	
	friend void multiply<>(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);

	


};

#include "faust_spmat.hpp"


#endif
