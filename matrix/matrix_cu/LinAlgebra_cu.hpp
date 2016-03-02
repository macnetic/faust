#ifndef __LINALGEBRA_CU_HPP__
#define __LINALGEBRA_CU_HPP__

#include <iostream>
#include "faust_cu_mat.h"
#include "faust_cu_vec.h"
#include "faust_cuda.h"
#include "kernels.h"

#ifdef __COMPILE_SPMAT__
   #include "faust_cu_spmat.h"
#endif
#include "faust_exception.h"


#ifdef __COMPILE_TIMERS__
	#include "faust_timer.h"
#endif

using namespace std;



template <typename faust_real>
void setOp(const faust_cu_mat<faust_real>& cu_A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2)
{
    if(opA == 'N')
    {
        opA1=cu_A.getNbRow();
        opA2=cu_A.getNbCol();
    }
    else if(opA == 'T')
    {
        opA1=cu_A.getNbCol();
        opA2=cu_A.getNbRow();
    }
    else
        handleError("LinAlgebra_cu","setOp : invalid character");
}


template <typename faust_real>
void add(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C)
{   
#ifdef __COMPILE_TIMERS__
cu_A.t_add_ext.start();
#endif
	if ((cu_A.getNbCol() != cu_B.getNbCol()) || (cu_A.getNbRow() != cu_B.getNbRow()))
		handleError("LinAlgebra_cu"," add : matrix dimension not equal");

	cu_C = cu_A;
	cu_C += cu_B;	

#ifdef __COMPILE_TIMERS__
cu_A.t_add_ext.stop();
#endif
}


// compute the biggest eigenvalue of A, A must be semi-definite positive 
template<typename T>	
T power_iteration(const  faust_cu_mat<T> & cu_A, const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag, cublasHandle_t cublasHandle)
{	

	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration.start();
	#endif


   const int nb_col = cu_A.getNbCol();
   int i = 0;
   flag = 0;
	 
   if (nbr_iter_max <= 0)
      handleError("LinAlgebra "," power_iteration :  nbr_iter_max <= 0");
   if (nb_col != cu_A.getNbRow())
      handleError("LinAlgebra "," power_iteration : faust_core<T> 1 must be a squared matrix"); 	
	 
   faust_cu_vec<T> cu_xk(nb_col);
   cu_xk.setOnes();
   faust_cu_vec<T> cu_xk_norm(nb_col);
   T lambda_old=1.0;
   T lambda = 0.0;
   while(fabs(lambda_old-lambda)>threshold && i<nbr_iter_max)
   {
      i++;
      lambda_old = lambda;
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_operator_equal.start();
	#endif
      cu_xk_norm = cu_xk;
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_operator_equal.stop();
      cu_A.t_power_iteration_normalize.start();
	#endif
      cu_xk_norm.normalize();
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_normalize.stop();
      cu_A.t_power_iteration_gemv.start();
	#endif
      gemv(cu_A, cu_xk_norm, cu_xk, cublasHandle);
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_gemv.stop();
      cu_A.t_power_iteration_dot.start();
	#endif
      lambda = cu_xk_norm.dot(cu_xk, cublasHandle);
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_dot.stop();
	#endif
      //std::cout << "i = " << i << " ; lambda=" << lambda << std::endl;
   }
   flag = (i<nbr_iter_max)?i:-1;

	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration.stop();
	#endif
   return lambda;
}


template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y,const faust_real alpha, const faust_real beta, char typeA, cublasHandle_t cublasHandle)
{
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.start();
	#endif
	faust_unsigned_int opA1,opA2;
	int currentGPU;
	faust_cudaGetDevice(&currentGPU);
	
   setOp(cu_A, typeA, opA1, opA2);
   const cublasOperation_t opA = (typeA=='T')?CUBLAS_OP_T:CUBLAS_OP_N;
	
	if   (opA2 != cu_x.size() )
		handleError("LinAlgebra_cu", "gemv : dimension conflict  between matrix op(A) and input vector x");
	
	if ( (beta!=0)  &&  (cu_y.size() != opA1))
		handleError("LinAlgebra_cu", "gemv : dimension conflict  between matrix op(A) and output vector y");
	

	if (opA1==0)
	{
	   cu_y.resize(0);
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.stop();
	#endif
		return;
	}

	if(opA2==0 || cu_A.estNulle() || alpha==0.0)
	{
		if (beta == 0.0)
			cu_y.setZeros(opA1);
		else if (beta != 1.0)
			cu_y *= beta;
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.stop();
	#endif
		return;
	}

	if(cu_A.estIdentite())
	{
		if(beta == 0.0)
		{
			cu_y = cu_x;
			if(alpha!=1.0) 
				cu_y *= alpha;
		}
		else
		{
			faust_cu_vec<faust_real> cu_tmp(cu_x);
			if(alpha!=1.0) 
				cu_tmp *= alpha;

         if(beta != 1.0)
			   cu_y *= beta;
			cu_y += cu_tmp;
		}
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.stop();
	#endif
		return;
	}
	

	int best_device = cu_A.getDevice();
   faust_cudaSetDevice(best_device);
	
	const faust_cu_mat<faust_real>* cu_A_ptr;
	if(cu_A.getDevice() == best_device)
		cu_A_ptr = &cu_A;
	else
		cu_A_ptr = new faust_cu_mat<faust_real>(cu_A, best_device);

	
	const faust_cu_vec<faust_real>* cu_x_ptr = &cu_x;
	if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
		cu_x_ptr = new faust_cu_vec<faust_real>(cu_x, best_device);

	faust_cu_vec<faust_real>* cu_y_ptr;
	if(cu_y.getDevice() == best_device)
		cu_y_ptr = &cu_y;
	else
		cu_y_ptr = new faust_cu_vec<faust_real>(cu_y, best_device);

   if(beta==0.0)
	   cu_y_ptr->resize(opA1);

	if(cu_A_ptr->getData()==NULL || cu_x_ptr->getData()==NULL || cu_y_ptr->getData()==NULL)
		handleError("LinAlgebra_cu", "gemv : A, x or y is not allocated");


		faust_cu_gemv(cublasHandle, opA,
                           opA1, opA2,
                           &alpha,
                           cu_A_ptr->getData(), cu_A_ptr->getNbRow(),
                           cu_x_ptr->getData(), 1,
                           &beta,
                           cu_y_ptr->getData(), 1);
	
	if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
		delete cu_x_ptr;
	if(cu_A.getDevice() != best_device)
		delete cu_A_ptr;
	if(cu_y.getDevice() != best_device)
	{
		cu_y.init(*cu_y_ptr, cu_y.getDevice());
		delete cu_y_ptr;
	}
		
	faust_cudaSetDevice(currentGPU);
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.stop();
	#endif
}




	
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const faust_real alpha, const faust_real beta, char  typeA, char  typeB, cublasHandle_t cublasHandle)
{

#ifdef __COMPILE_TIMERS__
cu_A.t_gemm.start();
#endif
	faust_unsigned_int opA1,opB1,opA2,opB2;

   setOp(cu_A, typeA, opA1, opA2);
   setOp(cu_B, typeB, opB1, opB2);
	cublasOperation_t opA = (typeA=='T')?CUBLAS_OP_T:CUBLAS_OP_N;
	cublasOperation_t opB = (typeB=='T')?CUBLAS_OP_T:CUBLAS_OP_N;


	if (opA2 != opB1)
		handleError("LinAlgebra_cu", "gemm : dimension conflict  between matrix op(A) and matrix op(B)");

	if ( (beta!=0) && (cu_C.getNbRow()!=opA1 || cu_C.getNbCol()!=opB2) )
		handleError("LinAlgebra_cu", "gemm : invalid dimension for output matrix C");



	if (opA1*opB2 == 0)
	{
		cu_C.resize(opA1,opB2);
		#ifdef __COMPILE_TIMERS__
		cu_A.t_gemm.stop();
		#endif
		return;
	}

	if(opA2==0 || cu_A.estNulle() || cu_B.estNulle() || alpha==0.0)
	{
		if (beta == 0.0)
			cu_C.setZeros(opA1,opB2);
		else if (beta != 1.0)
			cu_C *= beta;
		#ifdef __COMPILE_TIMERS__
		cu_A.t_gemm.stop();
		#endif
		return;
	}

	if(cu_A.estIdentite())
	{
		if(beta == 0.0)
		{
			cu_C = cu_B;
			if(typeB == 'T') 
				cu_C.transpose(cublasHandle);
			if(alpha!=1.0) 
				cu_C *= alpha;
		}
		else
		{
			faust_cu_mat<faust_real> cu_tmp(cu_B);
			if(typeB == 'T') 
				cu_tmp.transpose(cublasHandle);
			if(alpha!=1.0) 
				cu_tmp *= alpha;

         if(beta != 1.0)
			   cu_C *= beta;
			cu_C += cu_tmp;
		}
		#ifdef __COMPILE_TIMERS__
			cu_A.t_gemm.stop();
		#endif
		return;
	}

	if(cu_B.estIdentite())
	{
		if(beta == 0.0)
		{
			cu_C = cu_A;
			if(typeA == 'T') 
				cu_C.transpose(cublasHandle);
			if(alpha!=1.0) 
				cu_C *= alpha;
		}
		else
		{
			faust_cu_mat<faust_real> cu_tmp(cu_A);
			if(typeA == 'T') 
				cu_tmp.transpose(cublasHandle);
			if(alpha!=1.0) 
				cu_tmp *= alpha;

         if (beta != 1.0)
			   cu_C *= beta;
			cu_C += cu_tmp;
		}
		#ifdef __COMPILE_TIMERS__
			cu_A.t_gemm.stop();
		#endif
		return;
	}
	

	int currentGPU;
	faust_cudaGetDevice(&currentGPU);
	const int best_device = cu_A.getDevice();
   faust_cudaSetDevice(best_device);

	const faust_cu_mat<faust_real>* cu_A_ptr = &cu_A;
     if((cu_A.getDevice() != best_device) || (cu_A.getDevice() == best_device && cu_A == cu_C))
		cu_A_ptr = new faust_cu_mat<faust_real>(cu_A, best_device);

	const faust_cu_mat<faust_real>* cu_B_ptr = &cu_B;
	if((cu_B.getDevice() != best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C && cu_B!=cu_A))
		cu_B_ptr = new faust_cu_mat<faust_real>(cu_B, best_device);


	faust_cu_mat<faust_real>* cu_C_ptr;
	if(cu_C.getDevice() == best_device)
		cu_C_ptr = &cu_C;
	else
		cu_C_ptr = new faust_cu_mat<faust_real>(cu_C, best_device);


   if(beta == 0.0)
      cu_C_ptr->resize(opA1,opB2);
   else if(cu_C_ptr->estNulle())
      cu_C_ptr->setZeros(opA1,opB2);
   else if (cu_C_ptr->estIdentite())   
   {
      if(opA1 != opB2)
		   handleError("LinAlgebra_cu", "gemm : cu_C has been set to identite but is not a square matrix");
      cu_C_ptr->resize(opA1,opB2);
      kernel_memset(cu_C_ptr->getData(), (faust_real)0.0, cu_C_ptr->getNbRow()*cu_C_ptr->getNbCol());
      kernel_add_diag_const(cu_C_ptr->getData(), (faust_real)1.0, cu_C_ptr->getNbRow());
   }

	if(cu_A_ptr->getData()==NULL || cu_B_ptr->getData()==NULL || cu_C_ptr->getData()==NULL)
		handleError("LinAlgebra_cu", "gemm : A, B or C is not allocated");


   faust_cu_gemm(cublasHandle,
       opA, opB, opA1, opB2, opA2,
       &alpha, cu_A.getData(), cu_A.getNbRow(),
       cu_B_ptr->getData(), cu_B_ptr->getNbRow(),
       &beta, cu_C_ptr->getData(), cu_C_ptr->getNbRow());

	if((cu_A.getDevice() != best_device) || (cu_A.getDevice()==best_device && cu_A==cu_C))
		delete cu_A_ptr;
	if((cu_B.getDevice() != best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C && cu_B!=cu_A))
		delete cu_B_ptr;
	if(cu_C.getDevice() != best_device)
	{
		cu_C.init(*cu_C_ptr, cu_C.getDevice());
		delete cu_C_ptr;
	}
		
   faust_cudaSetDevice(currentGPU);

#ifdef __COMPILE_TIMERS__
cu_A.t_gemm.stop();
#endif
}




#ifdef __COMPILE_SPMAT__
///// FUNCTIONS with faust_cu_spmat matrices /////

template <typename faust_real>
void setOp(const faust_cu_spmat<faust_real>& cu_S, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2)
{
    if(opA == 'N')
    {
        opA1=cu_S.getNbRow();
        opA2=cu_S.getNbCol();
    }
    else if(opA == 'T')
    {
        opA1=cu_S.getNbCol();
        opA2=cu_S.getNbRow();
    }
    else
        handleError("setOp","gemm(const faust_cu_mat&, faust_cu_mat&, ...) : invalid character");
}

template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const  faust_real alpha, const faust_real beta, const char opA, cusparseHandle_t cusparseHandle)
{
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmv.start();
#endif
   faust_unsigned_int opA1,opA2;

   setOp(cu_A,opA,opA1,opA2);
   const cusparseOperation_t cusparseOpA = (opA=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;

   if(opA2 != cu_x.size())
        handleError("LinAlgebra_cu","gemv(const faust_cu_spmat&, const faust_cu_vec&, faust_cu_vec&, ...) : incorrect dimension between A and x");

	if ( (beta!=0)  &&  (cu_y.size() != opA1))
		handleError("LinAlgebra_cu", "gemv : dimension conflict  between matrix op(A) and output vector y");


   if (opA1==0)
   {
      cu_y.resize(0);
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmv.stop();
#endif
      return;
         
   }
   if (opA2==0 || cu_A.getNonZeros()==0 || alpha==0.0)
   {
      if(beta==0.0)
			cu_y.setZeros(opA1);
		else if (beta !=1.0)
			cu_y *= beta;
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmv.stop();
#endif
		return;
   }
         
        //handleError("LinAlgebra_cu","gemv(const faust_cu_vec&, faust_cu_vec&, ...) : vector x has not been initialized");
   else
   {


       // device of cu_A is used to compute
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      const int best_device = cu_A.getDevice();
      faust_cudaSetDevice(best_device);

      const faust_cu_vec<faust_real>* cu_x_ptr = &cu_x;
      if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
         cu_x_ptr = new faust_cu_vec<faust_real>(cu_x, best_device);

      const faust_cu_spmat<faust_real>* cu_A_ptr = &cu_A;
      if(cu_A.getDevice() != best_device)
         cu_A_ptr = new faust_cu_spmat<faust_real>(cu_A, best_device);


      faust_cu_vec<faust_real>* cu_y_ptr = &cu_y;
      if (cu_y.getDevice() != best_device)
         if(beta!=0.0)
            cu_y_ptr = new faust_cu_vec<faust_real>(cu_y, best_device);
         else
            cu_y_ptr = new faust_cu_vec<faust_real>(opA1, best_device);
      else if(beta==0.0)
         cu_y_ptr->resize(opA1);


      faust_cu_csrmv(cusparseHandle,cusparseOpA,
                cu_A_ptr->getNbRow(), cu_A_ptr->getNbCol(),
                cu_A_ptr->getNonZeros(),
                &alpha, cu_A_ptr->getDescr(), cu_A_ptr->getValues(),
                cu_A_ptr->getRowPtr(), cu_A_ptr->getColInd(),
                cu_x_ptr->getData(),
                &beta, cu_y_ptr->getData());


      if(cu_y.getDevice() != best_device)
      {
         cu_y = faust_cu_vec<faust_real>(*cu_y_ptr, cu_y.getDevice());
         delete cu_y_ptr ;
      }

      if(cu_A.getDevice() != best_device)
         delete cu_A_ptr;
      if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
         delete cu_x_ptr;

	   if (currentGPU != FAUST_GPU_ALREADY_SET)
		   faust_cudaSetDevice(currentGPU);
   }
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmv.stop();
#endif

}

template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const  faust_real alpha, const faust_real beta, const char opA, const char opB, cusparseHandle_t cusparseHandle)
{

#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.start();
#endif
    faust_unsigned_int opA1,opA2,opB1,opB2;


    if(opA=='T' && opB=='T')
        handleError("LinAlgebra_cu","gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : cannot compute transpose(A)*transpose(B). Try doing transpose(B*A) after converting A into a full matrix");

   setOp(cu_A,opA,opA1,opA2);
   setOp(cu_B,opB,opB1,opB2);
   const cusparseOperation_t cusparseOpA = (opA=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;
   const cusparseOperation_t cusparseOpB = (opB=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;




   if(opA2 != opB1)
        handleError("LinAlgebra_cu","gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : incorrect dimension between A and B");

	if ( (beta!=0) && (cu_C.getNbRow()!=opA1 || cu_C.getNbCol()!=opB2) )
		handleError("LinAlgebra_cu", "gemm : invalid dimension for output matrix C");



   if (opA1*opB2 == 0)
   {
	   cu_C.resize(opA1,opB2);
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.stop();
#endif
      return;
   }

   if(opA2==0 || cu_A.getNonZeros()==0 || cu_B.estNulle() || alpha==0.0)
   {
      if(beta == 0.0)
         cu_C.setZeros(opA1,opB2,cu_C.getDevice());
      else if(beta != 1.0)
         cu_C *= beta;
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.stop();
#endif
      return;
   }

   if(cu_B.estIdentite())
   {
      if (beta==0.0)
      {
         if (opA=='T') 
            cu_C.init_from_transpose(cu_A, cusparseHandle);
         else
            cu_C = cu_A;
         if (alpha!=1.0)
            cu_C *= alpha;
      }
      else
      {  
         faust_cu_spmat<faust_real> cu_S(cu_A, cu_C.getDevice());
         if (opA=='T') 
            cu_S.transpose(cusparseHandle);
         if (alpha!=1.0)
            cu_S *= alpha;

         cu_C *= beta;
         cu_C += cu_S;
      }   
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.stop();
#endif
      return;
   }
   else if (cu_B.getData()==NULL)
      handleError("LinAlgebra_cu","gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : matrix B has not been initialized");
        

   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   const int best_device = cu_B.getDevice();
   faust_cudaSetDevice(best_device);

   const faust_cu_mat<faust_real>* cu_B_ptr = &cu_B;
	if((cu_B.getDevice()!=best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C))
      cu_B_ptr = new faust_cu_mat<faust_real>(cu_B, best_device);

   const faust_cu_spmat<faust_real>* cu_A_ptr = &cu_A;
   if(cu_A.getDevice() != best_device)
      cu_A_ptr = new faust_cu_spmat<faust_real>(cu_A, best_device);

   faust_cu_mat<faust_real>* cu_C_ptr = &cu_C;
   if  (beta!=0)
   {
      if(cu_C.getDevice() != best_device)
      cu_C_ptr = new faust_cu_mat<faust_real>(cu_C, best_device);
   }
   else if(cu_C.getDevice() != best_device)
      cu_C_ptr = new faust_cu_mat<faust_real>(opA1, opB2, best_device);
   else
      cu_C_ptr->resize(opA1, opB2, best_device);



   faust_cu_csrmm(cusparseHandle,cusparseOpA,cusparseOpB,
         cu_A_ptr->getNbRow(), cu_C_ptr->getNbCol(), cu_A_ptr->getNbCol(),
         cu_A_ptr->getNonZeros(),
         &alpha, cu_A_ptr->getDescr(), cu_A_ptr->getValues(),
         cu_A_ptr->getRowPtr(), cu_A_ptr->getColInd(),
         cu_B_ptr->getData(), cu_B_ptr->getNbRow(),
         &beta, cu_C_ptr->getData(), cu_C_ptr->getNbRow());


   if(cu_C_ptr->getDevice()!=cu_C.getDevice())
   {
      cu_C = faust_cu_mat<faust_real>(*cu_C_ptr, cu_C.getDevice());
      delete cu_C_ptr ;
   }

   if(cu_A.getDevice() != best_device)
      delete cu_A_ptr;

	if((cu_B.getDevice()!=best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C))
      delete cu_B_ptr;

   faust_cudaSetDevice(currentGPU);
    
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.stop();
#endif
}
#endif



////// FUNCTIONS with faust_cu_mat matrices //////

// y = op(A) * x
template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, char typeA, cublasHandle_t cublasHandle)
{gemv(cu_A, cu_x, cu_y, faust_real(1.0), faust_real(0.0), typeA, cublasHandle);}
// y = A * x
template <typename faust_real>
void gemv(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cublasHandle_t cublasHandle)
{gemv(cu_A, cu_x, cu_y, (faust_real)1.0, (faust_real)0.0, 'N', cublasHandle);}
// y = A * x
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cublasHandle_t cublasHandle)
{gemv(cu_A, cu_x, cu_y, faust_real(1.0), faust_real(0.0), 'N', cublasHandle);}


// C = op(A) * op(B)
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, char typeA, char typeB, cublasHandle_t cublasHandle)
{gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), typeA, typeB, cublasHandle);}
// C = A * B
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle)
{gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), 'N', 'N', cublasHandle);}
// C = A * B
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle)
{gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), 'N', 'N', cublasHandle);}

//////////////////////////////////////////////////



#ifdef __COMPILE_SPMAT__


///// FUNCTIONS with faust_cu_spmat matrices /////

// y = op(A) * x
template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, const char opA, cusparseHandle_t cusparseHandle)
{gemv(cu_A, cu_x, cu_y, faust_real(1.0), faust_real(0.0), opA, cusparseHandle);}
// y = A * x
template <typename faust_real>
void gemv(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cusparseHandle_t cusparseHandle)
{gemv(cu_A, cu_x, cu_y, faust_real(1.0), faust_real(0.0), 'N', cusparseHandle);}
// y = A * x
template <typename faust_real>
void multiply(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_vec<faust_real>& cu_x, faust_cu_vec<faust_real>& cu_y, cusparseHandle_t cusparseHandle)
{gemv(cu_A, cu_x, cu_y, faust_real(1.0), faust_real(0.0), 'N', cusparseHandle);}

// C = op(A) * op(B) (with A sparse)
template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const char opA, const char opB, cusparseHandle_t cusparseHandle)
{gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), opA, opB, cusparseHandle);}
// C = A * B (with A sparse)
template <typename faust_real>
void gemm(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cusparseHandle_t cusparseHandle)
{gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), 'N', 'N', cusparseHandle);}
template <typename faust_real>
void multiply(const faust_cu_spmat<faust_real>& cu_A, const faust_cu_mat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cusparseHandle_t cusparseHandle)
{ gemm(cu_A, cu_B, cu_C, faust_real(1.0), faust_real(0.0), 'N', 'N', cusparseHandle);}

// C = alpha*op(A)*op(B) + beta*C (with B sparse)
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const faust_real alpha, const faust_real beta, const char opA, const char opB, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)
{char opAt=(opA=='N')?'T':'N';char opBt=(opB=='N')?'T':'N';gemm(cu_B, cu_A, cu_C, alpha, beta, opBt, opAt, cusparseHandle);cu_C.transpose(cublasHandle);}
// C = op(A) * op(B) (with B sparse)
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, const char opA, const char opB, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)
{char opAt=(opA=='N')?'T':'N';char opBt=(opB=='N')?'T':'N';gemm(cu_B, cu_A, cu_C, faust_real(1.0), faust_real(0.0), opBt, opAt, cusparseHandle);cu_C.transpose(cublasHandle);}
// C = A * B (with B sparse)
template <typename faust_real>
void gemm(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)
{gemm(cu_B, cu_A, cu_C, faust_real(1.0), faust_real(0.0), 'T', 'T', cusparseHandle);cu_C.transpose(cublasHandle);}
// C = A * B (with B sparse)
template <typename faust_real>
void multiply(const faust_cu_mat<faust_real>& cu_A, const faust_cu_spmat<faust_real>& cu_B, faust_cu_mat<faust_real>& cu_C, cublasHandle_t cublasHandle, cusparseHandle_t cusparseHandle)
{ gemm(cu_B, cu_A, cu_C, faust_real(1.0), faust_real(0.0), 'T', 'T', cusparseHandle);cu_C.transpose(cublasHandle);}

//////////////////////////////////////////////////
#endif // __COMPILE_SPMAT__


#endif	
