/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2018):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
#ifndef __FAUST_MAT_SPARSE_GPU_H__
#define __FAUST_MAT_SPARSE_GPU_H__

#include "faust_constant_gpu.h"
#include "faust_exception.h"
#include <vector>

template <typename FPP,Device DEVICE> class Vect;
template <typename FPP,Device DEVICE> class MatDense;
#include "faust_cuda.h"
#ifdef __COMPILE_TIMERS__
  #include "faust_Timer_gpu.h"
#endif


#include "faust_SpBlasHandle_gpu.h"

//! \class Faust::class MatSparse<FPP,Gpu> faust_MatSparse_gpu.h
//! \brief Class template representing sparse matrix for GPU processing <br>
//! This class implements sparse matrix multiplication <br>
//! The sparse matrix format is Compressed Column Storage (equivalent of the ColMajor storage for dense matrix).
//! \param FPP scalar numeric type, e.g float or double




//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    template <typename FPP,Device DEVICE> class MatSparse;
    template <typename FPP,Device DEVICE> class MatDense;

    template <typename FPP> class MatSparse<FPP,Gpu>
    {
        private:
            void _create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const int device_);
            void _create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int faust_unsigned_dim2_);
            void _clear();

        public:
            MatSparse();
            MatSparse(const int* csrRowPtr_, const int* csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            MatSparse(const MatSparse<FPP,Gpu>& cu_S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            MatSparse(const MatSparse<FPP,Cpu>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            MatSparse(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            MatSparse(const Faust::MatDense<FPP,Cpu>& A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_);
            void resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol);
            void copyFromHost(const int* csrRowPtr_, const int*  csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void copyFromDevice(const int* csrRowPtr_, const int* csrColInd_, const FPP* csrValues_, const faust_unsigned_int nnz_,  const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void copyToHost(int* csrRowPtr_, int* csrColInd_, FPP* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream=0)const;
            void copyToDevice(int* csrRowPtr_, int* csrColInd_, FPP* csrValues, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const;
            void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);


            void operator=(const MatSparse<FPP,Cpu>& S);
            void init(const MatSparse<FPP,Gpu>& cu_S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void init(const MatSparse<FPP,Cpu>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void init(const Faust::MatDense<FPP,Gpu>& cu_A, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
            void init(const Faust::MatDense<FPP,Cpu>& M, Faust::SpBlasHandle<Gpu> spblasHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);




          bool operator==(const MatSparse<FPP,Gpu>& cu_S)const;
          bool operator!=(const MatSparse<FPP,Gpu>& cu_S)const{return !((*this)==cu_S);}
            void operator+=(const FPP alpha);
            void operator-=(const FPP alpha);
            void operator*=(const FPP alpha);
            void operator/=(const FPP alpha);
            void setZeros(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
            void setIdentity(const faust_unsigned_int dim1_);


        public:
            void transpose(Faust::SpBlasHandle<Gpu> spblasHandle);
            void init_from_transpose(const MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> spblasHandle);
            FPP norm() const;
            void operator= (const MatSparse<FPP,Gpu>& M);


            faust_unsigned_int getNbRow()const{return dim1;}
            faust_unsigned_int getNbCol()const{return dim2;}
            faust_unsigned_int getNonZeros()const{return nnz;}
            int* getRowPtr(){if(csrRowPtr==NULL) handleError(m_className, "getRowPtr : non-allocated GPU pointer");return csrRowPtr;}
            int* getColInd(){if(csrColInd==NULL) handleError(m_className, "getColInd : non-allocated GPU pointer");return csrColInd;}
            FPP* getValues(){if(csrValues==NULL) handleError(m_className, "getValues : non-allocated GPU pointer");return csrValues;}
            const int* getRowPtr()const{if(csrRowPtr==NULL) handleError(m_className, "getRowPtr : non-allocated GPU pointer");return csrRowPtr;}
            const int* getColInd()const{if(csrColInd==NULL) handleError(m_className, "getColInd : non-allocated GPU pointer");return csrColInd;}
            const FPP* getValues()const{if(csrValues==NULL) handleError(m_className, "getRowPtr : non-allocated GPU pointer");return csrValues;}
            const cusparseMatDescr_t& getDescr()const{return *descr;}
            int getDevice()const{return device;}

          void Display()const;
        void print_file(const char* filename)const;
            void print_file(const char* filename, std::ios_base::openmode mode)const;


            ~MatSparse(){resize(0,0,0);}


        private:
            faust_unsigned_int dim1;
            faust_unsigned_int dim2;
            faust_unsigned_int nnz;
            int* csrRowPtr;
            int* csrColInd;
            FPP* csrValues;
            int device;
          cusparseMatDescr_t* descr;
          cusparseMatDescr_t descr_content;
            static const char * m_className;



    #ifdef __COMPILE_TIMERS__
      public:

          //temporary members
          static Faust::Timer_gpu t_constructor_from_device;
          static Faust::Timer_gpu t_constructor_from_host;
          static Faust::Timer_gpu t_create;
          static Faust::Timer_gpu t_clear;
          static Faust::Timer_gpu t_copy_from_host;
          static Faust::Timer_gpu t_copy_from_device;
          static Faust::Timer_gpu t_copy_to_host;
          static Faust::Timer_gpu t_copy_to_device;
          static Faust::Timer_gpu t_move_to_device;
          static Faust::Timer_gpu t_init_from_cuspmat;
          static Faust::Timer_gpu t_init_from_spmat;
          static Faust::Timer_gpu t_init_from_cumat;
          static Faust::Timer_gpu t_init_from_mat;
          static Faust::Timer_gpu t_scalar_multiply;
          static Faust::Timer_gpu t_operator_plus_equal_real;
          static Faust::Timer_gpu t_transpose;
          static Faust::Timer_gpu t_init_from_transpose;
          static Faust::Timer_gpu t_norm;
          static Faust::Timer_gpu t_display;
          static Faust::Timer_gpu t_print_file;

          static Faust::Timer_gpu t_csrmv;
          static Faust::Timer_gpu t_csrmm;

          void print_timers()const;
    #endif

    };

}

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::operator=(const Faust::MatSparse<FPP,Gpu>& cu_S)
{ init(cu_S, cu_S.device); }

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::operator=(const Faust::MatSparse<FPP,Cpu>& S)
{ init(S); }

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::operator-=(const FPP alpha)
{this->operator+=(-1.0*alpha);}

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::operator/=(const FPP alpha)
{
   if(alpha==0.0)
      handleError(m_className,"operator/= : dividing by 0");
   else
      this->operator*=(1.0/alpha);
}

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::setZeros(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{resize(0, dim1_, dim2_);}

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::setIdentity(const faust_unsigned_int dim1_)
{   handleError(m_className,"setIdentity : to set a matrix to identity, use the full-matrix format");}


template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol)
{resize(nnz_, nbRow, nbCol, device);}

template <typename FPP>
inline void Faust::MatSparse<FPP,Gpu>::_create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{_create(nnz_, dim1_, dim2_, device);}


#include "faust_MatSparse_gpu.hpp"

#endif
