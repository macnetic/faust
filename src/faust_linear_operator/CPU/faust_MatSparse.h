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

#include "faust_constant.h"
#include "faust_MatDiag.h"
#include "faust_MatDense.h"
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <vector>
#include "faust_MatGeneric.h"
#include "faust_Vect.h"
#include "faust_SpBlasHandle.h"


// modif AL AL
#include "faust_Transform.h"
#include "faust_TransformHelper.h"
#include "matio.h"
#include <random>
//! \class Faust::MatSparse<FPP,Cpu> faust_MatSparse.h
//! \brief Class template representing sparse matrix <br>
//! This class implements sparse matrix multiplication <br>
//! The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
//! \param FPP scalar numeric type, e.g float or double
//!

#include "faust_MatDiag.h"
//template<typename FPP> class MatDiag;

//! Faust::MatSparse class template of sparse matrix
template<typename FPP,Device DEVICE> class MatSparse;

//! Faust::Vect class template of dense vector
template<typename FPP,Device DEVICE> class Vect;

template<typename FPP,Device DEVICE> class Transform;

template<Device DEVICE> class SpBlasHandle;

template<typename FPP>
void Faust::multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

//! modif NB v1102 : comment useless function

template<typename FPP>
void Faust::spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);

template<typename FPP>
void Faust::wht_factors(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  factors, const bool);

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

	template<typename FPP,Device DEVICE> class MatGeneric;


	template<typename FPP> class MatDiag;

	//! Faust::MatDense class template of dense matrix
	template<typename FPP,Device DEVICE> class MatDense;

	template<typename FPP>
		class MatSparse<FPP,Cpu> : public Faust::MatGeneric<FPP,Cpu>
		{

			friend Faust::TransformHelper<FPP,Cpu>; // TODO: limit to needed member functions only
			friend void Faust::wht_factors<>(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  factors, const bool);
			friend class MatDense<FPP,Cpu>;
			//friend void MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S);

			public:
			MatSparse();
			MatSparse(const Faust::MatDiag<FPP>&);



			void faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C, const FPP & alpha, const FPP & beta, char typeA, char typeB)
				const{Faust::spgemm((*this),B,C,alpha,beta,typeA,typeB);}
			//!  \brief Constructor<br>
			//! Faust::MatSparse is a copy of an other Faust::MatSparse
			template<typename FPP1>
				MatSparse(const MatSparse<FPP1,Cpu>& M){(this)->operator=(M);}

			//!  \brief Constructor<br>
			//! Faust::MatSparse is a copy of an Faust::MatDense (dense matrix)
			template<typename FPP1>
				MatSparse(const Faust::MatDense<FPP1,Cpu>& M){(this)->operator=(M);}

			MatSparse(const MatSparse<FPP,Cpu>& M);
			MatSparse(const Faust::MatDense<FPP,Cpu>& M);

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
			void setZeros(){mat.setZero();nnz=0;}
			void setEyes(){mat.setIdentity();update_dim();}
			void transpose();
			void conjugate();
			FPP norm() const {return mat.norm();}
			void operator= (const MatSparse<FPP,Cpu>& M);
			void operator= (const Faust::MatDense<FPP,Cpu>& Mdense);
			void init (const Faust::MatDense<FPP,Cpu>& Mdense,Faust::SpBlasHandle<Cpu> spblas_handle)
			{(*this)=Mdense;}
			void operator*=(const FPP alpha);
			void operator/=(const FPP alpha);
			void setCoeff(const int& i, const int& j, const FPP & val);

			//! \brief check if the dimension and number of nonzeros of the Faust::MatSparse are coherent */
			void check_dim_validity() const;

			template <typename FPP1>
				void operator=(const MatSparse<FPP1,Cpu> &M);
			template <typename FPP1>
				void operator=(const Faust::MatDense<FPP1,Cpu> &M){MatSparse<FPP1,Cpu> spM(M);(this)->operator=(spM);}

			//int getNbRow()const{return this->dim1;}
			//int getNbCol()const{return this->dim2;}
			faust_unsigned_int getNonZeros()const{return nnz;}
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


			//! Display all features of Faust::MatSparse : dim1, dim2, number of non-zeros, values (if not more than 100), etc.
			void Display() const;

			//! Returns a string with all features of Faust::MatSparse: dim1, dim2, number of non-zeros, values (if not more than 100), etc.
			std::string to_string(const bool transpose=false, const bool displaying_small_mat_elts=false) const;

			//! \brief Display the support of Faust::MatSparse (i.e where are the non zero entries)
			void display_support() const;

			FPP norm(){return mat.norm();}

			//! \brief Write Faust::MatSparse into text file
			//! \param filename : name of the file
			//! The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients <br>
			//! All the other line contains 2 integers and one number :  the row, the column and the value of one coefficient in ColMajor access of the Faust::MatDense<br>
			void print_file(const char* filename)const;
			void print_file(const char* filename, std::ios_base::openmode mode)const;
			void init_from_file(FILE* fp);

			//! \brief Initialyse Faust::MatSparse from text file
			//!\tparam filename : name of the file
			//! The first line of the file contains 3 integers : the number of row and the number of column and the number of nonzeros coefficients. <br>
			//! All the other line contains 2 integers and one number : the row, the column and the value of one coefficient in ColMajor access of the Faust::MatDense
			void init_from_file(const char* filename);
			void init(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);


			void multiply(Faust::Vect<FPP,Cpu> & vec, char opThis='N') const
			{ vec.multiplyLeft((*this),opThis);}

			Faust::Vect<FPP,Cpu> multiply(const Faust::Vect<FPP,Cpu> & vec) const
			{
				return *this * vec;
			}
			//! \brief compute MatSparse-MatDense multiplication
			//! \param M : the dense matrix
			//! \param opThis : character
			//! M = (*this) * M if opThis='N'
			// M = (*this)' * M if opThis='T'
			void multiply(Faust::MatDense<FPP,Cpu> & M, char opThis) const;
			matvar_t* toMatIOVar(bool transpose, bool conjugate) const;
			//! \brief Converts the Matrix to a matio variable, especially useful for writing into a file with libmatio. The variable is encoded in full matrix format (on the contrary to toMatVarIOVar() which keeps the sparse representation).
			// \param transpose: set to true to obtain the matio variable for the transpose Matrix.
			// \param conjugate: set it to true to obtain the matio variable for the conjugate Matrix.
			// \return The matio variable matvar_t if it succeeded or NULL otherwise.
			// \see Faust::Transform::save_mat_file()
			// \see toMatIOVar()
			matvar_t* toMatIOVarDense(bool transpose, bool conjugate) const;
			FPP normL1(const bool transpose=false) const;
			FPP normL1(faust_unsigned_int&, const bool transpose=false) const;
			Faust::Vect<FPP,Cpu> get_col(faust_unsigned_int id) const;
			Faust::MatSparse<FPP,Cpu>* get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const;
			Faust::MatSparse<FPP,Cpu>* get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const;

			Faust::MatSparse<FPP,Cpu>* get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const;
			Faust::MatSparse<FPP,Cpu>* get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const;
static MatSparse<FPP, Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, double density);
			//\param : per_row means the density applies for each line rather than globally for the matrix
			static MatSparse<FPP, Cpu>* randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, double density, bool per_row);
			static MatSparse<FPP, Cpu>* eye(faust_unsigned_int num_rows, faust_unsigned_int num_cols);
			//! Destructor
			~MatSparse(){/*std::cout<<"destructor MatSparse"<<std::endl;*//*this->mat.resize(0,0);*/}

			private:
			void update_dim(){this->dim1=mat.rows();this->dim2=mat.cols();nnz=mat.nonZeros();}
			static const char * m_className;


			private:
			Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat;

			//! number of non-zero
			faust_unsigned_int nnz;

			//  ** modif NB v1102 ** : comment friend function
			//friend void Faust::MatDense<FPP,Cpu>::operator=(MatSparse<FPP,Cpu> const& S);

			//! *this = (*this) * S
			//friend void Faust::MatDense<FPP,Cpu>::operator*=(const MatSparse<FPP,Cpu>& S);
			//! *this = (*this) + S
			//friend void Faust::MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S);

			//! *this = (*this) - S
			//friend void Faust::MatDense<FPP,Cpu>::operator-=(const MatSparse<FPP,Cpu>& S);





			//! *this = S * (*this)
			friend void  Faust::Vect<FPP,Cpu>::multiplyLeft(MatSparse<FPP,Cpu> const& S,const char TransS);
			friend double Faust::Transform<FPP,Cpu>::normL1(const bool transpose) const;

			/*friend void  Faust::MatDense<FPP,Cpu>::multiplyLeft(MatSparse<FPP,Cpu> const& S,const char TransS);*/

			// MODIF AL WARNING, ERROR WITH VISUAL STUDIO 2013 compiler
			//friend void Faust::multiply<>(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

			//! modif NB v1102 : comment useless function

			friend void Faust::spgemm<>(const MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);

		};

}

#include "faust_MatSparse.hpp"


#endif
