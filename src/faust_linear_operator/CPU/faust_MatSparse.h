#ifndef __FAUST_MATSPARSE_H__
#define __FAUST_MATSPARSE_H__

#include "faust_constant.h"
#include "faust_MatDense.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "faust_MatGeneric.h"
#include "faust_Vect.h"
#include "faust_SpBlasHandle.h"


// modif AL AL
#include "faust_Transform.h"


//! \class Faust::MatSparse<FPP,Cpu> faust_MatSparse.h
//! \brief Class template representing sparse matrix <br>
//! This class implements sparse matrix multiplication <br>
//! The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
//! \param FPP scalar numeric type, e.g float or double
//!

//! Faust::MatDense class template of dense matrix
template<typename FPP,Device DEVICE> class MatDense;

//! Faust::MatSparse class template of sparse matrix
template<typename FPP,Device DEVICE> class MatSparse;

//! Faust::Vect class template of dense vector
template<typename FPP,Device DEVICE> class Vect;

template<typename FPP,Device DEVICE> class Transform;

template<Device DEVICE> class SpBlasHandle;

template<typename FPP>
void Faust::multiply(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);

template<typename FPP>
void Faust::spgemm(const Faust::MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template<typename FPP,Device DEVICE> class MatGeneric;

    template<typename FPP>
    class MatSparse<FPP,Cpu> : public Faust::MatGeneric<FPP,Cpu>
    {

        public:
        MatSparse();


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

        //! \brief Constructor : from CRS (Compressed Row Storage) format
        MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP* value, const size_t* id_row, const size_t* col_ptr);


        //!  \brief Constructor : from CCS (Compressed Column Storage) format
        template<typename FPP1>
        MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col);

        void set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const double* value, const size_t* id_row, const size_t* col_ptr);
        void resize(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
        void resize(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_){mat.resize(dim1_,dim2_);update_dim();}
        void setZeros(){mat.setZero();nnz=0;}
        void setEyes(){mat.setIdentity();update_dim();}
        void transpose();
        FPP norm() const {return mat.norm();}
        void operator= (const MatSparse<FPP,Cpu>& M);
        void operator= (const Faust::MatDense<FPP,Cpu>& Mdense);
        void init (const Faust::MatDense<FPP,Cpu>& Mdense,Faust::SpBlasHandle<Cpu> spblas_handle)
        {(*this)=Mdense;}
        void operator*=(const FPP alpha);
        void operator/=(const FPP alpha);

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

        const int* getRowPtr()const{if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(m_className,"getRowPtr : matrix is not in rowMajor");}}
        const int* getColInd()const{if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(m_className,"getColInd : matrix is not in rowMajor");}}
        int* getRowPtr(){if(mat.IsRowMajor) return mat.outerIndexPtr(); else{handleError(m_className,"getRowPtr : matrix is not in rowMajor");}}
        int* getColInd(){if(mat.IsRowMajor) return mat.innerIndexPtr(); else{handleError(m_className,"getColInd : matrix is not in rowMajor");}}
        bool isRowMajor() const{return mat.IsRowMajor;}

        //! Display all features of Faust::MatSparse : dim1, dim2, nnz number of nonzeros, values, etc ...
        void Display() const;

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

        //! Destructor
        ~MatSparse(){}

        private:
        void update_dim(){this->dim1=mat.rows();this->dim2=mat.cols();nnz=mat.nonZeros();}
        static const char * m_className;


        private:
        Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat;

        //! number of non-zero
        faust_unsigned_int nnz;

        friend void Faust::MatDense<FPP,Cpu>::operator=(MatSparse<FPP,Cpu> const& S);

        //! *this = (*this) * S
        friend void Faust::MatDense<FPP,Cpu>::operator*=(const MatSparse<FPP,Cpu>& S);
        //! *this = (*this) + S
        friend void Faust::MatDense<FPP,Cpu>::operator+=(const MatSparse<FPP,Cpu>& S);
        //! *this = (*this) - S
        friend void Faust::MatDense<FPP,Cpu>::operator-=(const MatSparse<FPP,Cpu>& S);

        // friend void sp_solve<>(const MatSparse<FPP,Cpu> & A,Faust::Vect<FPP,Cpu> & x, const Faust::Vect<FPP,Cpu> & y);

        //! *this = S * (*this)
        friend void Faust::MatDense<FPP,Cpu>::multiplyLeft(const MatSparse<FPP,Cpu>& S);

        //! *this = S * (*this)
        friend void  Faust::Vect<FPP,Cpu>::multiplyLeft(MatSparse<FPP,Cpu> const& A);

        friend void Faust::multiply<>(const Faust::Transform<FPP,Cpu> & A, const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, char typeA, char typeMult);
        friend void Faust::spgemm<>(const MatSparse<FPP,Cpu> & A,const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char  typeA, char  typeB);

    };

}

#include "faust_MatSparse.hpp"


#endif
