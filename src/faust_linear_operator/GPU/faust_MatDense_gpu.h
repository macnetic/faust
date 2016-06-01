#ifndef MATDENSE_GPU_H
#define MATDENSE_GPU_H


#include "faust_linear_algebra_gpu.h"
#include "faust_constant_gpu.h"

#include "faust_exception.h"
#include <vector>
#include <iterator>
#include "cuda_runtime.h"
#ifdef __COMPILE_TIMERS__
  #include "faust_Timer_gpu.h"
#endif



#include "faust_Vect_gpu.h"
// #include "faust_MatDense.h"
#ifdef __COMPILE_SPMAT__
	#include "faust_MatSparse.h"
	#include "faust_MatSparse_gpu.h"
#endif


#include "faust_reduce_gpu.h"
#include "kernels.h"
#include "faust_cuda.h"
#include "faust_gpu2cpu.h"
#include "faust_gpuCore2cpuCore.h"


//template<typename FPP, Device DEVICE> class MatDense;
//template <typename FPP,Device DEVICE> class Vect;



#ifdef __COMPILE_SPMAT__
    #include "cusparse.h"
    template <typename FPP,Device DEVICE> class MatSparse;
    #include "faust_MatSparse_gpu.h"
#endif

//! \class Faust::class MatDense<FPP,Gpu> faust_MatDense_gpu.h
//! \brief Class template representing dense matrix for GPU processing <br>


//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust
{

    // modif AL AL
    template<typename FPP,Device DEVICE>
    class MatDense;

    template <typename FPP> class MatDense<FPP,Gpu>
    {
    public:
        static const char * name;

        /// Constructeurs ///
        MatDense() : dim1(0), dim2(0), isIdentity(false), isZeros(false), device(FAUST_DEFAULT_CUDA_DEVICE), data(NULL)  {}
        MatDense(const FPP  *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
        MatDense(const MatDense<FPP,Gpu>& cu_mat, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
        MatDense(const MatDense<FPP,Cpu>& M, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );

    #ifdef __COMPILE_SPMAT__
        MatDense(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu> cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
        MatDense(const Faust::MatSparse<FPP,Gpu>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
    #endif
        MatDense(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE);


      // GETTEUR SETTEUR //
      faust_unsigned_int getNbRow() const {return dim1;}
      faust_unsigned_int getNbCol() const {return dim2;}

      int getDevice()const{return device;}

      void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_);
      void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol){resize(nbRow, nbCol, device);}
      void resize(const faust_unsigned_int nbRow){resize(nbRow, nbRow, device);}
      void copyFromHost(const FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
      void copyFromDevice(const FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
      void copyToHost(FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream=0)const;
      void copyToDevice(FPP *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const;
      void init(const MatDense<FPP,Gpu>& cu_M, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
      void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
      void check_dim_validity() const{}// empty method only useful for compatibility CPU-GPU matrix

      // (*this) = la matrice nulle
      void setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const  int device_);
      void setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol){setZeros(nbRow,nbCol,device);}
      void setZeros(){setZeros(dim1,dim2,device);}

      // (*this) = identite, pas forcement carree
      void setEyes(const faust_unsigned_int nbRow, const int device_);
      void setEyes(const faust_unsigned_int nbRow){setEyes(nbRow,device);}
      void setEyes();

      //FPP& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}

      //const FPP& operator[](faust_unsigned_int i)const{return mat.data()[i];}

      //const FPP& operator()(faust_unsigned_int i)const{return mat.data()[i];}
      //const FPP& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*dim1+i];}


      //transposition
      void transpose(Faust::BlasHandle<Gpu>);
      void init_from_transpose(const MatDense<FPP,Gpu>& cu_M, Faust::BlasHandle<Gpu>);
    #ifdef __COMPILE_SPMAT__
      void init_from_transpose(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu>);
    #endif

      // multiply (*this) = (*this) * A
      //void multiplyRight(const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::BlasHandle<Gpu>);
      //void operator*=(const MatDense<FPP>& cu_M){multiplyRight(cu_M);}
      // multiply (*this) =  A * (*this)
      //void multiplyLeft(const Faust::MatSparse<FPP,Gpu>& cu_A, Faust::BlasHandle<Gpu>);

      FPP max() const;
      FPP min() const;
      void abs();
      // frobenius norm
      FPP norm() const;
      // scalarMultiply (*this) = (*this) * lambda
      void scalarMultiply(FPP const lambda);
      void operator*=(FPP lambda){scalarMultiply(lambda);}
      void operator/=(FPP lambda){scalarMultiply(1.0/lambda);}

      bool operator==(const MatDense<FPP,Gpu>& cu_M)const;
      bool operator!=(const MatDense<FPP,Gpu>& cu_M)const{return !((*this)==cu_M);}
      void operator=( const MatDense<FPP,Gpu>& cu_A);
      void operator=(const MatDense<FPP,Cpu>& A);
    #ifdef __COMPILE_SPMAT__
      void operator=(const Faust::MatSparse<FPP,Cpu>& S);
      void operator=(const Faust::MatSparse<FPP,Gpu>& cu_S);
      void init_from_cu_spmat(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu>, const FPP coeff=1.0);


      void operator+=(const Faust::MatSparse<FPP,Gpu>& cu_S);
      void add(const Faust::MatSparse<FPP,Gpu>& S, Faust::SpBlasHandle<Gpu>);
      // operator-=(const Faust::MatSparse<FPP,Gpu>&) is not defined. Use sub(const Faust::MatSparse<FPP,Gpu>&, Faust::SpBlasHandle<Gpu>)
      void operator-=(const Faust::MatSparse<FPP,Gpu>& cu_S);
      void sub(const Faust::MatSparse<FPP,Gpu>& S, Faust::SpBlasHandle<Gpu>);
    #endif

      void operator+=(const MatDense<FPP,Gpu>& cu_A);
      void add(const MatDense<FPP,Gpu>& cu_A){this->operator+=(cu_A);}

      void operator-=(const MatDense<FPP,Gpu>& cu_A);
    void sub(const MatDense<FPP,Gpu>& cu_A){this->operator-=(cu_A);}

      // (*this)(i,j)=((*this)(i,j)) * A(i,j)
      void scalarMultiply(const MatDense<FPP,Gpu>& cu_A);
    #ifdef __COMPILE_SPMAT__
      void multiplyRight(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::BlasHandle<Gpu>, Faust::SpBlasHandle<Gpu>);
      void multiplyLeft(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::SpBlasHandle<Gpu>);
    #endif

       void Display()const;
       void print_file(const char* filename)const;

      ~MatDense(){resize(0);}

      FPP* getData(){return data;}
      const FPP* getData()const {return data;}

      // return the maximum of all coefficients of this and puts in row_id and col_id its location

      // spectral norm, "norm2", equal to the largest singular value
      //FPP spectralNorm() const;
      FPP spectralNorm(const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag, Faust::BlasHandle<Gpu> cublasHandle) const;

      // trace
      FPP trace() const;


      ////////////////// friends //////////////////////
      //template <typename FPP1> friend void gemm(const MatDense<FPP1,Gpu>& cu_A, const MatDense<FPP1,Gpu>& cu_B, MatDense<FPP1,Gpu>& cu_C, const FPP1 alpha, const FPP1 beta, char  typeA, char  typeB, BlasHandle<Gpu> cublasHandle);

      bool estIdentite()const{return isIdentity;}
      bool estNulle()const{return isZeros;}

    private :
       void _create(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_);
       void _clear();



      private:
         faust_unsigned_int dim1;
         faust_unsigned_int dim2;

         FPP* data;

         int device;

    //Eigen::Matrix<FPP,0,0> mat;
         bool isIdentity;
         bool isZeros;
         static const char * class_name;



    #ifdef __COMPILE_TIMERS__
      public:

          //temporary members
          static Faust::Timer_gpu t_constructor_from_device;
          static Faust::Timer_gpu t_constructor_from_host;
          static Faust::Timer_gpu t_constructor_from_size;
          static Faust::Timer_gpu t_create;
          static Faust::Timer_gpu t_clear;
          static Faust::Timer_gpu t_copy_from_host;
          static Faust::Timer_gpu t_copy_from_device;
          static Faust::Timer_gpu t_copy_to_host;
          static Faust::Timer_gpu t_copy_to_device;
          static Faust::Timer_gpu t_init;
          static Faust::Timer_gpu t_move_to_device;
          static Faust::Timer_gpu t_set_zeros;
          static Faust::Timer_gpu t_set_eyes;
          static Faust::Timer_gpu t_transpose;
          static Faust::Timer_gpu t_init_from_transpose;
          static Faust::Timer_gpu t_max;
          static Faust::Timer_gpu t_min;
          static Faust::Timer_gpu t_abs;
          static Faust::Timer_gpu t_norm;
          static Faust::Timer_gpu t_trace;
          static Faust::Timer_gpu t_spectral_norm2;
          static Faust::Timer_gpu t_scalar_multiply;
          static Faust::Timer_gpu t_operator_equal_from_device;
          static Faust::Timer_gpu t_operator_equal_from_host;
          static Faust::Timer_gpu t_operator_plus_equal;
          static Faust::Timer_gpu t_operator_less_equal;
          static Faust::Timer_gpu t_add_cuspmat;
          static Faust::Timer_gpu t_sub;
          static Faust::Timer_gpu t_hadamard_product;
          static Faust::Timer_gpu t_display;
          static Faust::Timer_gpu t_print_file;

          static Faust::Timer_gpu t_gemm;
          static Faust::Timer_gpu t_gemv;
          static Faust::Timer_gpu t_add_ext;
          static Faust::Timer_gpu t_power_iteration;
          static Faust::Timer_gpu t_power_iteration_operator_equal;
          static Faust::Timer_gpu t_power_iteration_normalize;
          static Faust::Timer_gpu t_power_iteration_gemv;
          static Faust::Timer_gpu t_power_iteration_dot;

          void print_timers()const;
    #endif
    };

}


/*template <typename FPP>
inline void MatDense<FPP,Gpu>::multiplyRight(const MatDense<FPP>& cu_B, Faust::BlasHandle<Gpu> cublasHandle)
{gemm(*this, cu_B, *this, 1.0, 0.0, 'N', 'N', cublasHandle);}

template <typename FPP,Device DEVICE>
inline void Faust::MatSparse<FPP,Gpu>::multiplyLeft(const Faust::MatSparse<FPP,Gpu>& cu_A, Faust::BlasHandle<Gpu> cublasHandle)
{gemm(cu_A, *this, *this, 1.0, 0.0, 'N', 'N', cublasHandle);}


#ifdef __COMPILE_SPMAT__
template <typename FPP,Device DEVICE>
   inline void Faust::MatSparse<FPP,Gpu>::multiplyRight(const Faust::MatSparse<FPP,Gpu>& cu_S, Faust::BlasHandle<Gpu> cublasHandle, cusparseHandle_t cusparseHandle)
   {gemm(*this, cu_S, *this, 1.0, 0.0, 'N', 'N', cublasHandle, cusparseHandle);}

template <typename FPP,Device DEVICE>
   inline void Faust::MatSparse<FPP,Gpu>::multiplyLeft(const Faust::MatSparse<FPP,Gpu>& cu_S, cusparseHandle_t cusparseHandle)
   {gemm(cu_S, *this, *this, 1.0, 0.0, 'N', 'N', cusparseHandle);}
#endif
*/
#include "faust_MatDense_gpu.hpp"

#endif
