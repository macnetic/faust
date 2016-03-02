#ifndef __FAUST_CU_SPMAT_H__
#define __FAUST_CU_SPMAT_H__

#include "faust_constant.h"
#include "faust_exception.h"
#include <vector>

template <typename faust_real> class faust_cu_vec;
template <typename faust_real> class faust_cu_mat;
#include "faust_cuda.h"
#ifdef __COMPILE_TIMERS__
  #include "faust_cu_timer.h"
#endif

#include "cusparse.h"

template <typename faust_real> class faust_spmat;
template <typename faust_real> class faust_mat;

template <typename faust_real> class faust_cu_spmat
{
	private:
		void _create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const int device_);	
		void _create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int faust_unsigned_dim2_);	
		void _clear();

	public:
		faust_cu_spmat();
		faust_cu_spmat(const int* csrRowPtr_, const int* csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, bool dataFromGPU, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		faust_cu_spmat(const faust_cu_spmat<faust_real>& cu_S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		faust_cu_spmat(const faust_spmat<faust_real>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0); 
		faust_cu_spmat(const faust_cu_mat<faust_real>& cu_A, cusparseHandle_t cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		faust_cu_spmat(const faust_mat<faust_real>& A, cusparseHandle_t cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, const int device_);
		void resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol);
		void copyFromHost(const int* csrRowPtr_, const int*  csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void copyFromDevice(const int* csrRowPtr_, const int* csrColInd_, const faust_real* csrValues_, const faust_unsigned_int nnz_,  const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, int srcDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void copyToHost(int* csrRowPtr_, int* csrColInd_, faust_real* csrValues_, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, cudaStream_t stream=0)const;
		void copyToDevice(int* csrRowPtr_, int* csrColInd_, faust_real* csrValues, const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0)const;
		void moveToDevice(int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);


		void operator=(const faust_spmat<faust_real>& S);
		void init(const faust_cu_spmat<faust_real>& cu_S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void init(const faust_spmat<faust_real>& S, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void init(const faust_cu_mat<faust_real>& cu_A, cusparseHandle_t cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);
		void init(const faust_mat<faust_real>& M, cusparseHandle_t cusparseHandle, int dstDevice=FAUST_DEFAULT_CUDA_DEVICE, cudaStream_t stream=0);


		//void gemm(const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const char opA, const char opB, const  faust_real alpha, const faust_real beta, cusparseHandle_t cusparseHandle);
		//void gmv(const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const char opA, const  faust_real alpha, const faust_real beta, cusparseHandle_t cusparseHandle);

      bool operator==(const faust_cu_spmat<faust_real>& cu_S)const;
      bool operator!=(const faust_cu_spmat<faust_real>& cu_S)const{return !((*this)==cu_S);}
		void operator+=(const faust_real alpha);
		void operator-=(const faust_real alpha);
		void operator*=(const faust_real alpha);
		void operator/=(const faust_real alpha);
		void setZeros(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_);
		void setIdentity(const faust_unsigned_int dim1_);


	public:
		void transpose(cusparseHandle_t cusparseHandle);
		void init_from_transpose(const faust_cu_spmat<faust_real>& cu_S, cusparseHandle_t cusparseHandle);
		faust_real norm() const;
		void operator= (const faust_cu_spmat<faust_real>& M);
		//void operator= (const faust_cu_mat<faust_real>& Mdense);
		
		faust_unsigned_int getNbRow()const{return dim1;}
		faust_unsigned_int getNbCol()const{return dim2;}
		faust_unsigned_int getNonZeros()const{return nnz;}
		int* getRowPtr(){if(csrRowPtr==NULL) handleError(class_name, "getRowPtr : non-allocated GPU pointer");return csrRowPtr;}
		int* getColInd(){if(csrColInd==NULL) handleError(class_name, "getColInd : non-allocated GPU pointer");return csrColInd;}
		faust_real* getValues(){if(csrValues==NULL) handleError(class_name, "getValues : non-allocated GPU pointer");return csrValues;}
		const int* getRowPtr()const{if(csrRowPtr==NULL) handleError(class_name, "getRowPtr : non-allocated GPU pointer");return csrRowPtr;}
		const int* getColInd()const{if(csrColInd==NULL) handleError(class_name, "getColInd : non-allocated GPU pointer");return csrColInd;}
		const faust_real* getValues()const{if(csrValues==NULL) handleError(class_name, "getRowPtr : non-allocated GPU pointer");return csrValues;}
		const cusparseMatDescr_t& getDescr()const{return *descr;}
		int getDevice()const{return device;}
		
      void Display()const;
   	void print_file(const char* filename)const;
		void print_file(const char* filename, std::ios_base::openmode mode)const;


		~faust_cu_spmat(){resize(0,0,0);}


	private:
		faust_unsigned_int dim1;
		faust_unsigned_int dim2;
		faust_unsigned_int nnz;
		int* csrRowPtr;
		int* csrColInd;
		faust_real* csrValues;			
		int device;
      cusparseMatDescr_t* descr;
      cusparseMatDescr_t descr_content;
		static const char * class_name;

	//friend void faust_cu_mat::operator=(const faust_cu_spmat<faust_real>& S);

	// *this = (*this) + S
	//friend void faust_cu_mat::operator+=(const faust_cu_spmat<faust_real>& cu_S);
	// *this = (*this) - S
	//friend void faust_cu_mat::operator-=(const faust_cu_spmat<faust_real>& cu_S);
	//friend void sp_solve(const faust_cu_spmat<faust_real>& cu_A,faust_cu_vec<faust_real>& cu_x, const faust_cu_vec<faust_real>& cu_y);

	// *this = S * (*this) 
	//friend void faust_cu_mat::multiplyLeft(const faust_cu_spmat<faust_real>& S, cusparseHandle_t);

	// *this = S * (*this) 
	//friend void  faust_cu_vec::multiplyLeft(const faust_cu_spmat<faust_real>& A, cusparseHandle_t);
	//friend void multiply(const faust_core<faust_real>& A, const faust_cu_mat<faust_real>& B, faust_cu_mat<faust_real>& C,const faust_real & alpha, char typeA, char typeMult);

#ifdef __COMPILE_TIMERS__
  public:

      //temporary members
      static faust_cu_timer t_constructor_from_device;
      static faust_cu_timer t_constructor_from_host;
      static faust_cu_timer t_create;
      static faust_cu_timer t_clear;
      static faust_cu_timer t_copy_from_host;
      static faust_cu_timer t_copy_from_device;
      static faust_cu_timer t_copy_to_host;
      static faust_cu_timer t_copy_to_device;
      static faust_cu_timer t_move_to_device;
      static faust_cu_timer t_init_from_cuspmat;
      static faust_cu_timer t_init_from_spmat;
      static faust_cu_timer t_init_from_cumat;
      static faust_cu_timer t_init_from_mat;
      static faust_cu_timer t_scalar_multiply;
      static faust_cu_timer t_operator_plus_equal_real;
      static faust_cu_timer t_transpose;
      static faust_cu_timer t_init_from_transpose;
      static faust_cu_timer t_norm;
      static faust_cu_timer t_display;
      static faust_cu_timer t_print_file;

      static faust_cu_timer t_csrmv;
      static faust_cu_timer t_csrmm;

      void print_timers()const;
#endif
	
};

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::operator=(const faust_cu_spmat<faust_real>& cu_S)
{ init(cu_S, cu_S.device); }

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::operator=(const faust_spmat<faust_real>& S)
{ init(S); }

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::operator-=(const faust_real alpha)
{this->operator+=(-1.0*alpha);}

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::operator/=(const faust_real alpha)
{
   if(alpha==0.0)
      handleError(class_name,"operator/= : dividing by 0");
   else
      this->operator*=(1.0/alpha);
}

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::setZeros(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{resize(0, dim1_, dim2_);}

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::setIdentity(const faust_unsigned_int dim1_)
{   handleError(class_name,"setIdentity : to set a matrix to identity, use the full-matrix format");}


template <typename faust_real>
inline void faust_cu_spmat<faust_real>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int nbRow, const faust_unsigned_int nbCol)
{resize(nnz_, nbRow, nbCol, device);}

template <typename faust_real>
inline void faust_cu_spmat<faust_real>::_create(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{_create(nnz_, dim1_, dim2_, device);}	


#include "faust_cu_spmat.hpp"

#endif
