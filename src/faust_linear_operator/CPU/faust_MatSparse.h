/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_MATSPARSE_H__
#define __FAUST_MATSPARSE_H__

#include <algorithm>
#include "faust_constant.h"
#include "faust_MatDiag.h"
#include "faust_MatDense.h"
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <vector>
#include "faust_MatGeneric.h"
#include "faust_Vect.h"
#include <iostream>


#include "faust_Transform.h"
#include "faust_linear_algebra.h"
#ifndef NO_MATIO
#include "matio.h" // otherwise matvar_t is defined as void from MatGeneric header
#endif
#include <random>
#include <vector>
#include "faust_WHT.h"
//! \class Faust::MatSparse<FPP,Cpu> faust_MatSparse.h
//! \brief Class template representing sparse matrix <br>
//! This class implements sparse matrix multiplication <br>
//! The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
//! \param FPP scalar numeric type, e.g float or double
//!

#include "faust_MatDiag.h"
//template<typename FPP> class MatDiag;

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{
	//! MatSparse class template of sparse matrix
	template<typename FPP,FDevice DEVICE> class MatSparse;

	//! Vect class template of dense vector
	template<typename FPP,FDevice DEVICE> class Vect;

	template<typename FPP,FDevice DEVICE> class Transform;




	template<typename FPP,FDevice DEVICE> class TransformHelper;
	template<typename FPP,FDevice DEVICE> class MatGeneric;


	template<typename FPP> class MatDiag;
	template<typename FPP, FDevice DEVICE> class MatBSR;

	//! MatDense class template of dense matrix
	template<typename FPP,FDevice DEVICE> class MatDense;
	template<typename FPP, FDevice DEVICE, typename FPP2> class EigTJ;
	template<typename FPP, FDevice DEVICE, typename FPP2> class EigTJParallel;
	template<typename FPP, FDevice DEVICE, typename FPP2> class EigTJComplex;
	template<typename FPP> class TransformHelperPoly;
	template<typename FPP, FDevice DEV> class MatButterfly;
	//TODO: simplify/remove the friendship by adding/using a public setter to is_ortho
	//template<typename FPP> void wht_factors(unsigned int n, std::vector<MatGeneric<FPP,Cpu>*>&  factors, const bool, const bool);
	template<typename FPP>
		void multiply(const Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);


	template<typename FPP>
		class MatSparse<FPP,Cpu> : public MatGeneric<FPP,Cpu>
		{

			friend MatBSR<FPP, Cpu>;
			friend EigTJ<FPP,Cpu, double>;
			friend EigTJParallel<FPP,Cpu, double>;
			friend EigTJ<FPP,Cpu, float>;
			friend EigTJParallel<FPP,Cpu, float>;
			friend EigTJComplex<FPP,Cpu, double>;
			friend EigTJComplex<FPP,Cpu, float>;
			friend Transform<FPP,Cpu>; //TODO: limit to needed member functions only (multiply)
			friend TransformHelper<FPP,Cpu>; // TODO: limit to needed member functions only
			friend TransformHelperPoly<FPP>; // TODO: limit to needed member functions only
			friend MatButterfly<FPP, Cpu>; // TODO: limit to needed member functions only
			friend void wht_factors<>(unsigned int n, std::vector<MatGeneric<FPP,Cpu>*>&  factors, const bool, const bool);
			friend class MatDense<FPP,Cpu>;
			friend class MatSparse<std::complex<double>, Cpu>;
			friend class MatSparse<std::complex<float>, Cpu>;
			friend class MatSparse<double, Cpu>;
			friend class MatSparse<float, Cpu>;
			//friend void MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S);

			public:
			MatSparse();
			MatSparse(const MatDiag<FPP>&);



			void faust_gemm(const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C, const FPP & alpha, const FPP & beta, char typeA, char typeB)
				const{spgemm((*this),B,C,alpha,beta,typeA,typeB);}
			//!  \brief Constructor<br>
			//! MatSparse is a copy of an other MatSparse
			template<typename FPP1>
				MatSparse(const MatSparse<FPP1,Cpu>& M){(this)->operator=(M);}

			//!  \brief Constructor<br>
			//! MatSparse is a copy of an MatDense (dense matrix)
			template<typename FPP1>
				MatSparse(const MatDense<FPP1,Cpu>& M){(this)->operator=(M);}

			MatSparse(const MatSparse<FPP,Cpu>& M);
			MatSparse(const MatDense<FPP,Cpu>& M);

			//! \brief Constructor
			//!	\param dim1_ : number of row of the matrix
			//!	\param dim2_ : number of column of the matrix
			MatSparse(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

			//!  \brief Constructor
			//!	\param nnz_  : number of non-zero
			//!	\param dim1_ : number of row of the matrix
			//!	\param dim2_ : number of column of the matrix
			MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

			//!  \brief Constructor
			//!	\param rowidx : vector<int>
			//!	\param colidx : vector<int>
			//!	\param values : vector<FPP,Cpu>
			//!	\param dim1_ : number of row of the matrix
			//!	\param dim2_ : number of column of the matrix
			MatSparse(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

			//! \brief Constructor
			//	WARNING: using this constructor is discounraged because rowidx, colidx, values are not necessarily safe, it's the responsibility of the caller to check their allocation space according to values.size().
			//!	\param rowidx : row indices with for all k < values.size(), M[rowidx[k]][colidx[k]] = values[k];
			//!	\param colidx : column indices with for all k < values.size(), M[rowidx[k]][colidx[k]] = values[k];
			//!	\param values : nonzeros of the matrix (should be nnz).
			//!	\param dim1_ : number of row of the matrix
			//!	\param dim2_ : number of column of the matrix
			MatSparse(const unsigned int* rowidx, const unsigned int* colidx, const FPP* values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, faust_unsigned_int nnz);


			//! \brief Constructor
			// \param tripletList: STL triplets (row_id, col_id, nonzero_value) used to define the matrix tripletList.size() is the nnz of the matrix.
			//!	\param dim1 : number of row of the matrix
			//!	\param dim2 : number of column of the matrix
			MatSparse(std::vector<Eigen::Triplet<FPP>> &tripletList, const faust_unsigned_int dim1, const faust_unsigned_int dim2);

			//!  \brief Constructor
			//	WARNING: using this constructor is discounraged because rowidx, colidx are not necessarily safe, its the responsibility of the caller to check their allocation space according to values.size().
			//!	\param rowidx : row indices with for all k < values.size(), M[rowidx[k]][colidx[k]] = values[k];
			//!	\param colidx : column indices with for all k < values.size(), M[rowidx[k]][colidx[k]] = values[k];
			//!	\param values : vector<FPP,Cpu>
			//!	\param dim1_ : number of row of the matrix
			//!	\param dim2_ : number of column of the matrix
			MatSparse(const unsigned int* rowidx, const unsigned int* colidx, const std::vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);

			//! \brief Constructor : from CRS (Compressed Row Storage) format
			MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP* value, const size_t* id_row, const size_t* col_ptr);


			//!  \brief Constructor : from CCS (Compressed Column Storage) format
			template<typename FPP1>
				MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col, const bool transpose = false);

			void set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP* value, const size_t* id_row, const size_t* col_ptr);
			void resize(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
			void resize(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_){mat.resize(dim1_,dim2_);update_dim();}
			/**
			 * Resizes without erasing the previous values contrary to what resize() does.
			 */
			void conservativeResize(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_){mat.conservativeResize(dim1_, dim2_); update_dim();}
			void setZeros(){mat.setZero();nnz=0;}
			void setEyes();
			void transpose();
			void conjugate(const bool eval=true);
			void adjoint();
			typename Eigen::NumTraits<FPP>::Real norm() const {return mat.norm();}
			void operator= (const MatSparse<FPP,Cpu>& M);
			void operator= (const MatDense<FPP,Cpu>& Mdense);
			void init (const MatDense<FPP,Cpu>& Mdense)
			{(*this)=Mdense;}
			void operator*=(const FPP alpha);
			void operator/=(const FPP alpha);
			void setCoeff(const int& i, const int& j, const FPP & val);

			//! \brief check if the dimension and number of nonzeros of the MatSparse are coherent
			void check_dim_validity() const;

			template <typename FPP1>
				void operator=(const MatSparse<FPP1,Cpu> &M);
			template <typename FPP1>
				void operator=(const MatDense<FPP1,Cpu> &M){MatSparse<FPP1,Cpu> spM(M);(this)->operator=(spM);}

			//int getNbRow()const{return this->dim1;}
			//int getNbCol()const{return this->dim2;}
			faust_unsigned_int getNonZeros()const{return nnz;}
			size_t getNBytes() const;
			static size_t getNBytes(int nnz, int nrows);
			bool isCompressedMode()const{return mat.isCompressed();}
			void makeCompression(){mat.makeCompressed();}

			//! return pointer value of length nnz
			FPP* getValuePtr(){return mat.valuePtr();}

			//! return const pointer value of length nnz
			const FPP* getValuePtr()const{return mat.valuePtr();}

			// if rowMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim1,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the row (n-1)
			// if colMajor : getOuterIndexPtr()[0]=0 ; for n=1 to dim2,  getOuterIndexPtr()[n] = getOuterIndexPtr()[n-1]+ number of non-zeros elements in the col (n-1)
			//!return row-index value of length equal to the number of row+1
			int* getOuterIndexPtr(){return mat.outerIndexPtr();}
			const int* getOuterIndexPtr()const{return mat.outerIndexPtr();}

			// if rowMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] = column index matching the element getValuePtr()[n];
			// if colMajor : for n=0 to (dim1-1), getInnerIndexPtr()[n] =   row  index matching the element getValuePtr()[n];
			//! return col index of length nnz
			int* getInnerIndexPtr(){return mat.innerIndexPtr();}
			const int* getInnerIndexPtr()const{return mat.innerIndexPtr();}

			//TODO: indices should be unsigned int
			const int* getRowPtr()const{if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(m_className,"getRowPtr : matrix is not in rowMajor");}}
			const int* getColInd()const{if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(m_className,"getColInd : matrix is not in rowMajor");}}
			int* getRowPtr(){if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(m_className,"getRowPtr : matrix is not in rowMajor");}}
			int* getColInd(){if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(m_className,"getColInd : matrix is not in rowMajor");}}
			bool isRowMajor() const{return mat.IsRowMajor;}

			void copyRowPtr(size_t* out_rowptr) const;
			void copyRowPtr(int* out_rowptr) const;
			void copyColInd(size_t* out_colInd) const;
			void copyColInd(int* out_colInd) const;
			void copyValuePtr(FPP* out_values) const;
			template <typename U>
			void copyBufs(U* out_rowptr, U* out_colind, FPP* out_values) const;

			// Virtual Function (inherited method from MatGeneric)
			MatType getType() const{ return MatType::Sparse;}


			/*!  \brief return a "copy" to the given matrix
			 *  \param isOptimize (optionnal) : boolean which the style of copy <br>
			 -True, the return copy is optimized for the product <br>
			 which means dynamic type of the copy could be different from the original one <br>
			 -False, the return copy is simple, the dynamic type isn't changed <br>   						(default value is False) <br>
			//! \return  a pointer of MatGeneric
			//  \warning the dynamic type of the copy can be different from the original object
			*/
			MatGeneric<FPP,Cpu>* Clone(const bool isOptimize=false) const;


			//! Display all features of MatSparse : dim1, dim2, number of non-zeros, values (if not more than 100), etc.
			void Display() const;

			//! Returns a string with all features of MatSparse: dim1, dim2, number of non-zeros, values (if not more than 100), etc.
			std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;

			std::string to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity);

			//! \brief Display the support of MatSparse (i.e where are the non zero entries)
			void display_support() const;

			Real<FPP> norm(){return mat.norm();}

			//! \brief Write MatSparse into text file
			//! \param filename : name of the file
			//! The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients <br>
			//! All the other line contains 2 integers and one number :  the row, the column and the value of one coefficient in ColMajor access of the MatDense<br>
			void print_file(const char* filename)const;
			void print_file(const char* filename, std::ios_base::openmode mode)const;
			void init_from_file(FILE* fp);

			//! \brief Initialyse MatSparse from text file
			//!\tparam filename : name of the file
			//! The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients. <br>
			//! All the other line contains 2 integers and one number : the row, the column and the value of one coefficient in ColMajor access of the MatDense
			void init_from_file(const char* filename);
			void init(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);


			void multiply(Vect<FPP,Cpu> & vec, char opThis='N') const
			{ vec.multiplyLeft((*this),opThis);}

			Vect<FPP,Cpu> multiply(const Vect<FPP,Cpu> & vec) const
			{
				return *this * vec;
			}
			//! \brief compute MatSparse-MatDense multiplication
			//! \param M : the dense matrix
			//! \param opThis : character
			//! M = (*this) * M if opThis='N'
			// M = (*this)' * M if opThis='T'
			void multiply(MatDense<FPP,Cpu> & M, char opThis) const;
			//! \brief compute MatSparse-MatSparse multiplication
			//! \param M : the dense matrix
			//! \param opThis : character
			//! M = (*this) * M if opThis='N'
			// M = (*this)' * M if opThis='T'
			void multiply(MatSparse<FPP,Cpu> & M, char opThis) const;
			matvar_t* toMatIOVar(bool transpose, bool conjugate, const char* var_name=nullptr) const;
			//!  \brief Creates a MatSparse from at matio variable
			void from_matio_var(matvar_t* var);
			//!  \brief Creates a MatSparse from at .mat file
			void read_from_mat_file(const char *filepath, const char *var_name);
			//!  \brief Saves a MatSparse to a .mat file
			void save_to_mat_file(const char *filepath, const char *var_name);

			//! \brief multiply this by M into this.
			void multiplyRight(MatSparse<FPP,Cpu> const & M);

			void scalarMultiply(const FPP& lambda);

			//! \brief Converts the Matrix to a matio variable, especially useful for writing into a file with libmatio. The variable is encoded in full matrix format (on the contrary to toMatVarIOVar() which keeps the sparse representation).
			// \param transpose: set to true to obtain the matio variable for the transpose Matrix.
			// \param conjugate: set it to true to obtain the matio variable for the conjugate Matrix.
			// \return The matio variable matvar_t if it succeeded or NULL otherwise.
			// \see Transform::save_mat_file()
			// \see toMatIOVar()
			matvar_t* toMatIOVarDense(bool transpose, bool conjugate) const;
			Real<FPP> normL1(const bool transpose=false) const;
			Real<FPP> normL1(faust_unsigned_int&, const bool transpose=false) const;
			Real<FPP> normInf(const bool transpose=false) const;
			Real<FPP> normInf(faust_unsigned_int&, const bool transpose=false) const;

			const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return const_cast<Eigen::SparseMatrix<FPP,Eigen::RowMajor>*>(&mat)->coeffRef(i,j);}



			void submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, MatDense<FPP, Cpu> & submat) const;
			void submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, FPP* submat_data) const;
			Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;
			void get_col(faust_unsigned_int id, Vect<FPP, Cpu> &out_vec) const;
			/*
			 * \brief Returns a slice of "this" as a MatSparse. The slice start from column start_col_id and finishes to the column start_col_id+num_cols-1 of "this".
			 *
			 * Warning: using this function is discouraged as it returns a naked pointer, use preferably another prototype of get_cols except if you really need this one
			 */

			MatSparse<FPP,Cpu>* get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const;
			/*
			 * \brief Returns a MatSparse composed of num_cols columns of "this". Their indices are defined in col_ids, the first column of the returned matrix is this[:, col_ids[0]], the second this[:, col_ids[1]], ...
			 *
			 * Warning: using this function is discouraged as it returns a naked pointer, use preferably another prototype of get_cols except if you really need this one
			 */
			MatSparse<FPP,Cpu>* get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;
			/* \brief Returns a column slice of "this" into out_cols. The slice start from column start_col_id and finishes to the column start_col_id+num_cols-1 of "this".
			 *
			 * \param out_cols: the MatSparse doesn't have to be initialized (it's handled internally).
			 */
			void get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols, MatSparse<FPP, Cpu>& out_cols) const;
			/*
			 * \brief Returns a MatSparse composed of num_cols columns of "this" into out_cols. Their indices are defined in col_ids, the first column of the returned matrix is this[:, col_ids[0]], the second this[:, col_ids[1]], ...
			 *
			 * \param out_cols: the MatSparse doesn't have to be initialized (it's handled internally).
			 *
			 * Warning: using this function is discouraged as it returns a naked pointer, use preferably another prototype of get_cols except if you really need this one
			 */
			void get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols, MatSparse<FPP, Cpu>& out_cols) const;

			void delete_col(faust_unsigned_int id);
			void delete_row(faust_unsigned_int id);

			MatSparse<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			MatSparse<FPP,Cpu>* get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
			void get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows, MatSparse<FPP, Cpu>& out_rows) const;
			void get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows, MatSparse<FPP, Cpu>& out_rows) const;

			void swap_rows(faust_unsigned_int id1, faust_unsigned_int id2);
			void swap_cols(faust_unsigned_int id1, faust_unsigned_int id2);

			std::vector<int> col_nonzero_inds(faust_unsigned_int col_id) const;
			std::vector<int> row_nonzero_inds(faust_unsigned_int row_id) const;

			std::list<std::pair<int,int>> nonzeros_indices(const double& tol=0) const;

			/**
			 * \brief Concatenates vertically top and bottom matrices, resizing this if necessary.
			 */
			void vstack(MatSparse<FPP, Cpu>& top, MatSparse<FPP, Cpu>& bottom);

			/**
			 * \brief Concatenates horizontally left and right matrices, resizing this if necessary.
			 */
			void hstack(MatSparse<FPP, Cpu>& left, MatSparse<FPP, Cpu>& right);

			void real(MatSparse<Real<FPP>, Cpu> & real_mat) const;

			template<typename FPP2>
				MatSparse<Real<FPP2>, Cpu> to_real() const;


			template<typename FPP2>
				MatSparse<FPP2, Cpu> cast() const;
			/**
			 * \brief index this(row_ids, col_ids) Eigen matrix and mutliplies the result by in_mat into out_mat (two Eigen dense matrices of respective types MatType1, MatTyp2)
			 *
			 */
			template<typename MatType1, typename MatType2>
				void eigenIndexMul(const faust_unsigned_int* row_ids, const faust_unsigned_int* col_ids, size_t nrows, size_t ncols, const MatType1 &in_mat, MatType2 &out_mat, bool transpose = false, bool conjugate = false);


			bool containsNaN() const;
			void print_bufs(const std::string name="");
			void print_asarray(const std::string name="");

			MatDense<FPP, Cpu> to_dense() const;

			static MatSparse<FPP, Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, Real<FPP> density);
			//\param : per_row means the density applies for each line rather than globally for the matrix
			static MatSparse<FPP, Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, Real<FPP> density, bool per_row);
			static MatSparse<FPP, Cpu>* eye(faust_unsigned_int num_rows, faust_unsigned_int num_cols);

			// Permutation matrix set to exchange two rows or columns (depending on the multiplication side)
			static MatSparse<FPP,Cpu>* swap_matrix(faust_unsigned_int order, faust_unsigned_int id1, faust_unsigned_int id2);
			//! Destructor
			~MatSparse(){/*std::cout<<"destructor MatSparse"<<std::endl;*//*this->mat.resize(0,0);*/}

			using SpMat = Eigen::SparseMatrix<FPP,Eigen::RowMajor>;
			private:
			void update_dim(){this->dim1=mat.rows();this->dim2=mat.cols();nnz=mat.nonZeros();}
			static const char * m_className;


			Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat;

			//! number of non-zero
			faust_unsigned_int nnz;

			//  ** modif NB v1102 ** : comment friend function
			//friend void MatDense<FPP,Cpu>::operator=(MatSparse<FPP,Cpu> const& S);

			//! *this = (*this) * S
			//friend void MatDense<FPP,Cpu>::operator*=(const MatSparse<FPP,Cpu>& S);
			//! *this = (*this) + S
			//friend void MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S);

			//! *this = (*this) - S
			//friend void MatDense<FPP,Cpu>::operator-=(const MatSparse<FPP,Cpu>& S);





			//! *this = S * (*this)
			friend void  Vect<FPP,Cpu>::multiplyLeft(MatSparse<FPP,Cpu> const& S,const char TransS);
			friend double Transform<FPP,Cpu>::normL1(const bool transpose, const bool full_array/*=true*/, const int batch_sz/*=1*/) const;

			/*friend void  MatDense<FPP,Cpu>::multiplyLeft(MatSparse<FPP,Cpu> const& S,const char TransS);*/

			// MODIF AL WARNING, ERROR WITH VISUAL STUDIO 2013 compiler
			//friend void multiply<>(const Transform<FPP,Cpu> & A, const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

			//! modif NB v1102 : comment useless function

			// cf. faust_linear_algebra.h
			friend void spgemm<>(const MatSparse<FPP,Cpu> & A,const MatDense<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
			friend void spgemm<>(const MatDense<FPP,Cpu> & A,const MatSparse<FPP,Cpu> & B, MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);
			friend void gemm_gen<>(const MatGeneric<FPP, Cpu>& A, const MatGeneric<FPP, Cpu>& B, MatDense<FPP, Cpu>& out, const FPP alpha, const FPP beta, const char opA, const char opB);
		};

	template<typename FPP>
		void copy_sp_mat(MatSparse<FPP,Cpu>& src, MatSparse<FPP, Cpu>& dst);

}

#include "faust_MatSparse.hpp"


#endif
