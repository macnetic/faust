#ifndef __FAUST_CU_MAT_H__
#define __FAUST_CU_MAT_H__


#include "faust_constant.h"
#include "faust_exception.h"
#include <vector>
#include <iterator>
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#ifdef __COMPILE_TIMERS__
  #include "faust_cu_timer.h"
#endif
#include "LinAlgebra_cu.h"


template <typename faust_real> class faust_vec;
template <typename faust_real> class faust_cu_vec;
template <typename faust_real> class faust_mat;

#ifdef __COMPILE_SPMAT__
  #include "cusparse.h"
   template <typename faust_real> class faust_spmat;
   #include "faust_cu_spmat.h"
#endif
template <typename faust_real> class faust_core;

template <typename faust_real> class faust_cu_mat
{
public:
   static const char * name;

  /// Constructeurs ///
  faust_cu_mat() : dim1(0), dim2(0), isIdentity(false), isZeros(false), device(FAUST_DEFAULT_CUDA_DEVICE), data(NULL)  {}
  faust_cu_mat(const faust_real  *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
  faust_cu_mat(const faust_cu_mat<faust_real>& cu_mat, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
  faust_cu_mat(const faust_mat<faust_real>& M, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );

#ifdef __COMPILE_SPMAT__
  faust_cu_mat(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
  faust_cu_mat(const faust_spmat<faust_real>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
#endif
  faust_cu_mat(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE);



  /// GETTEUR SETTEUR ///
  faust_unsigned_int getNbRow() const {return dim1;}
  faust_unsigned_int getNbCol() const {return dim2;}

  int getDevice()const{return device;}

  void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_);
  void resize(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol){resize(nbRow, nbCol, device);}
  void resize(const faust_unsigned_int nbRow){resize(nbRow, nbRow, device);}
  void copyFromHost(const faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
  void copyFromDevice(const faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
  void copyToHost(faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream=0)const;
  void copyToDevice(faust_real *data_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const;
  void init(const faust_cu_mat<faust_real>& cu_M, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
  void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);


  // (*this) = la matrice nulle
  void setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const  int device_);
  void setZeros(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol){setZeros(nbRow,nbCol,device);}
  void setZeros(){setZeros(dim1,dim2,device);}

  // (*this) = identite, pas forcement carree
  void setEyes(const faust_unsigned_int nbRow, const int device_);
  void setEyes(const faust_unsigned_int nbRow){setEyes(nbRow,device);}
  void setEyes();

  //faust_real& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}

  //const faust_real& operator[](faust_unsigned_int i)const{return mat.data()[i];}

  //const faust_real& operator()(faust_unsigned_int i)const{return mat.data()[i];}
  //const faust_real& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*dim1+i];}


  //transposition
  void transpose(cublasHandle_t);
  void init_from_transpose(const faust_cu_mat<faust_real>& cu_M, cublasHandle_t);
#ifdef __COMPILE_SPMAT__
  void init_from_transpose(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t);
#endif

  // multiply (*this) = (*this) * A
  void multiplyRight(const faust_cu_mat<faust_real>& cu_B, cublasHandle_t);
  //void operator*=(const faust_cu_mat<faust_real>& cu_M){multiplyRight(cu_M);}
  // multiply (*this) =  A * (*this)
  void multiplyLeft(const faust_cu_mat<faust_real>& cu_A, cublasHandle_t);

  faust_real max() const;
  faust_real min() const;
  void abs();
  // frobenius norm
  faust_real norm() const;
  // scalarMultiply (*this) = (*this) * lambda
  void scalarMultiply(faust_real const lambda);
  void operator*=(faust_real lambda){scalarMultiply(lambda);}
  void operator/=(faust_real lambda){scalarMultiply(1.0/lambda);}

  bool operator==(const faust_cu_mat<faust_real>& cu_M)const;
  bool operator!=(const faust_cu_mat<faust_real>& cu_M)const{return !((*this)==cu_M);}
  void operator=( const faust_cu_mat<faust_real>& cu_A);
  void operator=(const faust_mat<faust_real>& A);
#ifdef __COMPILE_SPMAT__
  void operator=(const faust_spmat<faust_real>& S);
  void operator=(const faust_cu_spmat<faust_real>& cu_S);
  void init_from_cu_spmat(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t, const faust_real coeff=1.0);
 
  // operator+=(const faust_cu_spmat<faust_real>&) is not defined. Use add(const faust_cu_spmat<faust_real>&, cusparseHandle_t)
  void operator+=(const faust_cu_spmat<faust_real>& cu_S);
  void add(const faust_cu_spmat<faust_real>& S, cusparseHandle_t);
  // operator-=(const faust_cu_spmat<faust_real>&) is not defined. Use sub(const faust_cu_spmat<faust_real>&, cusparseHandle_t)
  void operator-=(const faust_cu_spmat<faust_real>& cu_S);
  void sub(const faust_cu_spmat<faust_real>& S, cusparseHandle_t);
#endif

  void operator+=(const faust_cu_mat<faust_real>& cu_A);
  void add(const faust_cu_mat<faust_real>& cu_A){this->operator+=(cu_A);}

  void operator-=(const faust_cu_mat<faust_real>& cu_A);
  void sub(const faust_cu_mat<faust_real>& cu_A){this->operator-=(cu_A);}

  // (*this)(i,j)=((*this)(i,j)) * A(i,j)
  void scalarMultiply(const faust_cu_mat<faust_real>& cu_A);
#ifdef __COMPILE_SPMAT__
  void multiplyRight(const faust_cu_spmat<faust_real>& cu_S, cublasHandle_t, cusparseHandle_t);
  void multiplyLeft(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t);
#endif

   void Display()const;
   void print_file(const char* filename)const;

  ~faust_cu_mat(){resize(0);}

  faust_real* getData(){return data;}
  const faust_real* getData()const {return data;}


  // return the maximum of all coefficients of this and puts in row_id and col_id its location


  // spectral norm, "norm2", equal to the largest singular value
  //faust_real spectralNorm() const;
  faust_real spectralNorm(const faust_unsigned_int nbr_iter_max,faust_real threshold, faust_int & flag, cublasHandle_t cublasHandle) const;

  // trace
  faust_real trace() const;



  ////////////////// friends //////////////////////
  // intra classe//
  //friend void multiply(const faust_cu_mat<faust_real> & cu_A, const faust_cu_mat<faust_real> & cu_B, faust_cu_mat<faust_real> & cu_C, cublasHandle_t cublasHandle);
  template <typename T> friend void gemm(const faust_cu_mat<T>& cu_A, const faust_cu_mat<T>& cu_B, faust_cu_mat<T>& cu_C, const T alpha, const T beta, char  typeA, char  typeB, cublasHandle_t cublasHandle);
  //friend void add(const faust_cu_mat<faust_real> & cu_A, const faust_cu_mat<faust_real> & cu_B, faust_cu_mat<faust_real> & cu_C);
  //friend void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const faust_real alpha, const faust_real beta, char  typeA, char  typeB, cublasHandle_t cublasHandle);
  //friend void multiply(const faust_cu_core& A, const faust_cu_mat& cu_B, faust_cu_mat& cu_C,const faust_real alpha, char typeA, char typeMult);
  //friend void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec& x, faust_cu_vec& y,const faust_real alpha, const faust_real beta, char typeA, cublasHandle_t cublasHandle);
  //friend faust_vec solve(const faust_cu_mat<faust_real> & cu_A, const faust_cu_vec<faust_real> & v);
  ///////////friend faust_spmat::operator=(const faust_cu_mat<faust_real>& cu_S);
  bool estIdentite()const{return isIdentity;}
  bool estNulle()const{return isZeros;}

private : 
   void _create(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int device_);
   void _clear();



  private:
     faust_unsigned_int dim1;
     faust_unsigned_int dim2;

     faust_real* data;

     int device;

//Eigen::Matrix<faust_real,0,0> mat;
     bool isIdentity;
     bool isZeros;
     static const char * class_name;



#ifdef __COMPILE_TIMERS__
  public:

      //temporary members
      static faust_cu_timer t_constructor_from_device;
      static faust_cu_timer t_constructor_from_host;
      static faust_cu_timer t_constructor_from_size;
      static faust_cu_timer t_create;
      static faust_cu_timer t_clear;
      static faust_cu_timer t_copy_from_host;
      static faust_cu_timer t_copy_from_device;
      static faust_cu_timer t_copy_to_host;
      static faust_cu_timer t_copy_to_device;
      static faust_cu_timer t_init;
      static faust_cu_timer t_move_to_device;
      static faust_cu_timer t_set_zeros;
      static faust_cu_timer t_set_eyes;
      static faust_cu_timer t_transpose;
      static faust_cu_timer t_init_from_transpose;
      static faust_cu_timer t_max;
      static faust_cu_timer t_min;
      static faust_cu_timer t_abs;
      static faust_cu_timer t_norm;
      static faust_cu_timer t_trace;
      static faust_cu_timer t_spectral_norm2;
      static faust_cu_timer t_scalar_multiply;
      static faust_cu_timer t_operator_equal_from_device;
      static faust_cu_timer t_operator_equal_from_host;
      static faust_cu_timer t_operator_plus_equal;
      static faust_cu_timer t_operator_less_equal;
      static faust_cu_timer t_add_cuspmat;
      static faust_cu_timer t_sub;
      static faust_cu_timer t_hadamard_product;
      static faust_cu_timer t_display;
      static faust_cu_timer t_print_file;

      static faust_cu_timer t_gemm;
      static faust_cu_timer t_gemv;
      static faust_cu_timer t_add_ext;
      static faust_cu_timer t_power_iteration;
      static faust_cu_timer t_power_iteration_operator_equal;
      static faust_cu_timer t_power_iteration_normalize;
      static faust_cu_timer t_power_iteration_gemv;
      static faust_cu_timer t_power_iteration_dot;

      void print_timers()const;
#endif
};

template <typename faust_real>
inline void faust_cu_mat<faust_real>::multiplyRight(const faust_cu_mat<faust_real>& cu_B, cublasHandle_t cublasHandle)
{gemm(*this, cu_B, *this, 1.0, 0.0, 'N', 'N', cublasHandle);}

template <typename faust_real>
inline void faust_cu_mat<faust_real>::multiplyLeft(const faust_cu_mat<faust_real>& cu_A, cublasHandle_t cublasHandle)
{gemm(cu_A, *this, *this, 1.0, 0.0, 'N', 'N', cublasHandle);}


#ifdef __COMPILE_SPMAT__
template <typename faust_real>
   inline void faust_cu_mat<faust_real>::multiplyRight(const faust_cu_spmat<faust_real>& cu_S, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)
   {gemm(*this, cu_S, *this, 1.0, 0.0, 'N', 'N', cublasHandle, cusparseHandle);}

template <typename faust_real>
   inline void faust_cu_mat<faust_real>::multiplyLeft(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle)
   {gemm(cu_S, *this, *this, 1.0, 0.0, 'N', 'N', cusparseHandle);}
#endif

#include "faust_cu_mat.hpp"

#endif
