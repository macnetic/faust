#ifndef __FAUST_SPMAT_H__
#define __FAUST_SPMAT_H__

#include "faust_constant.h"
#include "faust_mat.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "faust_mat_generic.h"
#include "faust_vec.h"

/*! \class faust_spmat faust_spmat.h
* \brief Class template representing sparse matrix <br>
* This class implements sparse matrix multiplication <br>
* The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
*\tparam T scalar numeric type, e.g float or double
*/

/// faust_mat class template of dense matrix
template<typename T> class faust_mat;

/// faust_spmat class template of sparse matrix
template<typename T> class faust_spmat;

/// faust_vec class template of dense vector
template<typename T> class faust_vec;

template<typename T>
void multiply(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);


template<typename T>
void spgemm(const faust_spmat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB);



template<typename T>
class faust_spmat : public faust_mat_generic<T>
{
	public:	
	faust_spmat();
	void faust_gemm(const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB)const{spgemm((*this),B,C,alpha,beta,typeA,typeB);}
	/*!
	*  \brief Constructor<br>
	* faust_spmat is a copy of an other faust_spmat
	*/	
	template<typename U>
	faust_spmat(const faust_spmat<U>& M){(this)->operator=(M);}
	
	/*!
	*  \brief Constructor<br>
	* faust_spmat is a copy of an faust_mat (dense matrix)
	*/	
	template<typename U>
	faust_spmat(const faust_mat<U>& M){(this)->operator=(M);}

	faust_spmat(const faust_spmat<T>& M);
	faust_spmat(const faust_mat<T>& M);

	/*!
	*  \brief Constructor
	*	\tparam dim1_ : number of row of the matrix
		\tparam dim2_ : number of column of the matrix
	*/	
	faust_spmat(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

	/*!
	*  \brief Constructor
	*	\tparam nnz_  : number of non-zero
		\tparam dim1_ : number of row of the matrix
		\tparam dim2_ : number of column of the matrix
	*/	
	faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

	/*!
	*  \brief Constructor
	*	\tparam rowidx : vector<int>
		\tparam colidx : vector<int>
		\tparam values : vector<T>
		\tparam dim1_ : number of row of the matrix
		\tparam dim2_ : number of column of the matrix
	*/	
	faust_spmat(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
	/*!
	*  \brief Constructor : from CRS (Compressed Row Storage) format
	*	
	*/
	faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const T* value, const size_t* id_row, const size_t* col_ptr);
	

	/*!
	*  \brief Constructor : from CCS (Compressed Column Storage) format 
	*	
	*/
	template<typename U>
	faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const U* value, const int* row_ptr, const int* id_col);

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
	/*!
     	*\brief check if the dimension and number of nonzeros of the faust_spmat are coherent */
	void check_dim_validity() const;
	
	template <typename U>
	void operator=(const faust_spmat<U> &M);
	template <typename U>
	void operator=(const faust_mat<U> &M){faust_spmat<U> spM(M);(this)->operator=(spM);}

	//int getNbRow()const{return this->dim1;}
	//int getNbCol()const{return this->dim2;}
	faust_unsigned_int getNonZeros()const{return nnz;}
	bool isCompressedMode()const{return mat.isCompressed();}
	void makeCompression(){mat.makeCompressed();}
	
	/// return pointer value of length nnz
	T* getValuePtr(){return mat.valuePtr();} 

	/// return const pointer value of length nnz
	const T* getValuePtr()const{return mat.valuePtr();}

	// if rowMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim1,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the row (n-1)
	// if colMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim2,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the col (n-1)
	///return row-index value of length equal to the number of row+1
	int* getOuterIndexPtr(){return mat.outerIndexPtr();}
	const int* getOuterIndexPtr()const{return mat.outerIndexPtr();}

	// if rowMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] = column index matching the element getValuePtr()[n];
	// if colMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] =   row  index matching the element getValuePtr()[n];
	///return col index of length nnz
	int* getInnerIndexPtr(){return mat.innerIndexPtr();}
	const int* getInnerIndexPtr()const{return mat.innerIndexPtr();}
	
	const int* getRowPtr()const{if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(class_name,"getRowPtr : matrix is not in rowMajor");}}
	const int* getColInd()const{if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(class_name,"getColInd : matrix is not in rowMajor");}}
	int* getRowPtr(){if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(class_name,"getRowPtr : matrix is not in rowMajor");}}
	int* getColInd(){if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(class_name,"getColInd : matrix is not in rowMajor");}}
	bool isRowMajor() const{return mat.IsRowMajor;}
		
	/// Display all features of faust_spmat : dim1, dim2, nnz number of nonzeros, values, etc ...
	void Display() const;
	/*!
	*\brief Display the support of faust_spmat (i.e where are the non zero entries)
	*/
	void display_support() const;
	
	T norm(){return mat.norm();}	
	
	/*!
     	*\brief Write faust_spmat into text file 
	*\tparam filename : name of the file
	* The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients <br>
	* All the other line contains 2 integers and one number :  the row, the column and the value of one coefficient in ColMajor access of the faust_mat<br>
     	*/
	void print_file(const char* filename)const;
	void print_file(const char* filename, std::ios_base::openmode mode)const;
	void init_from_file(FILE* fp);
	
	/*!
     	*\brief Initialyse faust_spmat from text file 
	*\tparam filename : name of the file
	* The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients. <br>
	* All the other line contains 2 integers and one number : the row, the column and the value of one coefficient in ColMajor access of the faust_mat
	*/	
	void init_from_file(const char* filename);
	void init(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
	
	///Destructor
	~faust_spmat(){}

	private:
	void update_dim(){this->dim1=mat.rows();this->dim2=mat.cols();nnz=mat.nonZeros();}
	static const char * class_name;

	private:
	Eigen::SparseMatrix<T,Eigen::RowMajor> mat;	
	
	/// number of non-zero
	faust_unsigned_int nnz;

	friend void faust_mat<T>::operator=(faust_spmat<T> const& S);

	/// *this = (*this) * S
	friend void faust_mat<T>::operator*=(const faust_spmat<T>& S);
	/// *this = (*this) + S
	friend void faust_mat<T>::operator+=(const faust_spmat<T>& S);
	/// *this = (*this) - S
	friend void faust_mat<T>::operator-=(const faust_spmat<T>& S);
	
	// friend void sp_solve<>(const faust_spmat<T> & A,faust_vec<T> & x, const faust_vec<T> & y);

	/// *this = S * (*this) 
	friend void faust_mat<T>::multiplyLeft(const faust_spmat<T>& S);

	/// *this = S * (*this) 
	friend void  faust_vec<T>::multiplyLeft(faust_spmat<T> const& A);
		
	friend void multiply<>(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);
	
	friend void spgemm<>(const faust_spmat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, const T & beta, char  typeA, char  typeB);

};

#include "faust_spmat.hpp"


#endif
