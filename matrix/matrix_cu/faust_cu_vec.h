#ifndef __FAUST_CU_VEC_H__
#define __FAUST_CU_VEC_H__
 
#include "faust_constant.h"
#include <vector>
#include <iterator>
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif
#include "faust_exception.h"

template <typename faust_real> class faust_mat;
template <typename faust_real> class faust_cu_mat;
template <typename faust_real> class faust_vec;
#ifdef __COMPILE_SPMAT__
  #include "cusparse.h"
  template <typename faust_real> class faust_spmat;
  template <typename faust_real> class faust_cu_spmat;
#endif

template <typename faust_real>
class faust_cu_vec
{
 public :
 faust_cu_vec() : dim(0), device(FAUST_DEFAULT_CUDA_DEVICE), data(NULL) {}
  faust_cu_vec(const faust_real  *data_, const faust_unsigned_int dim_, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
  faust_cu_vec(const faust_cu_vec<faust_real>& cu_v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
  faust_cu_vec(const faust_vec<faust_real>& v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );
  faust_cu_vec(const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0 );

 void resize(const faust_unsigned_int dim_, int device_);
 void resize(const faust_unsigned_int dim_){resize(dim_, device);}
   
  void copyFromHost(const faust_real *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0); 
  void copyFromDevice(const faust_real *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0); 
  void copyToHost(faust_real *data_, const faust_unsigned_int dim_, cudaStream_t stream=0)const; 
  void copyToDevice(faust_real *data_, const faust_unsigned_int dim_, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const; 
  void init(const faust_cu_vec<faust_real>& cu_v, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
  void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);

void operator=(const faust_cu_vec<faust_real>& cu_v);
void operator=(const faust_vec<faust_real>& v);

faust_real max() const;
faust_real min() const;
faust_real sum()const;
faust_real norm() const;
faust_real mean()const{return sum()/dim;}
faust_real dot(const faust_cu_vec<faust_real>& cu_v, cublasHandle_t cublasHandle)const;

bool operator==(const faust_cu_vec<faust_real>& cu_v)const;
bool operator!=(const faust_cu_vec<faust_real>& cu_v)const{return !((*this)==cu_v);}

void operator+=(const faust_cu_vec<faust_real>& cu_v);
void operator-=(const faust_cu_vec<faust_real>& cu_v);
void operator*=(const faust_cu_vec<faust_real>& cu_v);
void operator/=(const faust_cu_vec<faust_real>& cu_v);

void operator+=(const faust_real alpha);
void operator-=(const faust_real alpha);
void operator*=(const faust_real alpha);
void operator/=(const faust_real alpha);
void scalarMultiply(const faust_real alpha){this->operator*=(alpha);}
void normalize(){scalarMultiply(1/norm());} 

void setValues(const faust_real);
void setValues(const faust_real, const faust_unsigned_int);

void setZeros();
void setZeros(const faust_unsigned_int);

void setOnes();
void setOnes(const faust_unsigned_int);


faust_real mean_relative_error(const faust_cu_vec<faust_real>& v)const;
// multiply (*this) =  A * (*this)
inline void  multiplyLeft(const faust_cu_mat<faust_real>& cu_A, cublasHandle_t cublasHandle);
#ifdef __COMPILE_SPMAT__
void  multiplyLeft(const faust_cu_spmat<faust_real>& cu_S,cusparseHandle_t);
#endif

~faust_cu_vec(){resize(0);}

faust_unsigned_int size() const {return dim;}

int getDevice()const{return device;} 

faust_real* getData(){return data;}
const faust_real* getData()const {return data;}

   void Display()const;


//friend void gemv(const faust_cu_mat<faust_real>& A,const faust_cu_vec<faust_real>& x, faust_cu_vec<faust_real>& y,const faust_real alpha, const faust_real beta, char typeA, cublasHandle_t cublasHandle);
//friend faust_cu_vec solve(const faust_mat<faust_real>& A, const faust_cu_vec<faust_real>& v);
#ifdef __COMPILE_SPMAT__
  //friend void sp_solve(const faust_spmat<faust_real>& A,faust_cu_vec<faust_real>& x, const faust_cu_vec<faust_real>& y);
#endif

private:
   void _create(const faust_unsigned_int dim_, int device_);
   void _clear();


private:
  faust_unsigned_int dim;
  static const char * class_name;
  faust_real* data;
  int device;   

   #ifdef __COMPILE_TIMERS__
      public: 
         faust_timer t_local_multiplyLeft;
   #endif
  
};

template<typename faust_real>
inline void faust_cu_vec<faust_real>::multiplyLeft(const faust_cu_mat<faust_real>& cu_A, cublasHandle_t cublasHandle)
{
    faust_cu_vec cu_y_copy;
    gemv(cu_A, *this, cu_y_copy, 1.0, 0.0, 'N', cublasHandle);
    this->operator=(cu_y_copy);
}

#include "faust_cu_vec.hpp"

#endif
