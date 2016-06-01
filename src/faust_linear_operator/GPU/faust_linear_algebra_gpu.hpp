#ifndef __LINALGEBRA_GPU_HPP__
#define __LINALGEBRA_GPU_HPP__

#define __COMPILE_SPMAT__

//modif AL
//#include "faust_constant.h"

using namespace std;


template <typename FPP>
void Faust::setOp(const Faust::MatDense<FPP,Gpu>& cu_A, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2)
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
        handleError("LinAlgebra_gpu","Faust::setOp : invalid character");
}


template <typename FPP>
void Faust::add(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C)
{
#ifdef __COMPILE_TIMERS__
cu_A.t_add_ext.start();
#endif
	if ((cu_A.getNbCol() != cu_B.getNbCol()) || (cu_A.getNbRow() != cu_B.getNbRow()))
		handleError("LinAlgebra_gpu"," Faust::add : matrix dimension not equal");

	cu_C = cu_A;
	cu_C += cu_B;

#ifdef __COMPILE_TIMERS__
cu_A.t_add_ext.stop();
#endif
}


// compute the biggest eigenvalue of A, A must be semi-definite positive
template<typename FPP>
FPP Faust::power_iteration(const  Faust::MatDense<FPP,Gpu> & cu_A, const faust_unsigned_int nbr_iter_max,FPP threshold, faust_int & flag, Faust::BlasHandle<Gpu> blasHandle)
{

	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration.start();
	#endif


   const int nb_col = cu_A.getNbCol();
   int i = 0;
   flag = 0;

   if (nbr_iter_max <= 0)
      handleError("linear_algebra "," Faust::power_iteration :  nbr_iter_max <= 0");
   if (nb_col != cu_A.getNbRow())
      handleError("linear_algebra "," Faust::power_iteration : Faust::Transform<T> 1 must be a squared matrix");

   Faust::Vect<FPP,Gpu> cu_xk(nb_col);
   cu_xk.setOnes();
   Faust::Vect<FPP,Gpu> cu_xk_norm(nb_col);
   FPP lambda_old=1.0;
   FPP lambda = 0.0;
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
      Faust::gemv(cu_A, cu_xk_norm, cu_xk, blasHandle);
	#ifdef __COMPILE_TIMERS__
		cu_A.t_power_iteration_gemv.stop();
      cu_A.t_power_iteration_dot.start();
	#endif
      lambda = cu_xk_norm.dot(cu_xk, blasHandle);
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


template <typename FPP>
void Faust::gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y,const FPP alpha, const FPP beta, char typeA, Faust::BlasHandle<Gpu> blasHandle)
{
	#ifdef __COMPILE_TIMERS__
		cu_A.t_gemv.start();
	#endif
	faust_unsigned_int opA1,opA2;
	int currentGPU;
	faust_cudaGetDevice(&currentGPU);

   Faust::setOp(cu_A, typeA, opA1, opA2);
   const cublasOperation_t opA = (typeA=='T')?CUBLAS_OP_T:CUBLAS_OP_N;

	if   (opA2 != cu_x.size() )
		handleError("LinAlgebra_gpu", "Faust::gemv : dimension conflict  between matrix op(A) and input vector x");

	if ( (beta!=0)  &&  (cu_y.size() != opA1))
		handleError("LinAlgebra_gpu", "Faust::gemv : dimension conflict  between matrix op(A) and output vector y");


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
			Faust::Vect<FPP,Gpu> cu_tmp(cu_x);
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

	const Faust::MatDense<FPP,Gpu>* cu_A_ptr;
	if(cu_A.getDevice() == best_device)
		cu_A_ptr = &cu_A;
	else
		cu_A_ptr = new Faust::MatDense<FPP,Gpu>(cu_A, best_device);


	const Faust::Vect<FPP,Gpu>* cu_x_ptr = &cu_x;
	if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
		cu_x_ptr = new Faust::Vect<FPP,Gpu>(cu_x, best_device);

	Faust::Vect<FPP,Gpu>* cu_y_ptr;
	if(cu_y.getDevice() == best_device)
		cu_y_ptr = &cu_y;
	else
		cu_y_ptr = new Faust::Vect<FPP,Gpu>(cu_y, best_device);

   if(beta==0.0)
	   cu_y_ptr->resize(opA1);

	if(cu_A_ptr->getData()==NULL || cu_x_ptr->getData()==NULL || cu_y_ptr->getData()==NULL)
		handleError("LinAlgebra_gpu", "Faust::gemv : A, x or y is not allocated");


		faust_cu_gemv(blasHandle.GetHandle(), opA,
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





template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const FPP alpha, const FPP beta, char  typeA, char  typeB, Faust::BlasHandle<Gpu> blasHandle)
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
		handleError("LinAlgebra_gpu", "Faust::gemm : dimension conflict  between matrix op(A) and matrix op(B)");

	if ( (beta!=0) && (cu_C.getNbRow()!=opA1 || cu_C.getNbCol()!=opB2) )
		handleError("LinAlgebra_gpu", "Faust::gemm : invalid dimension for output matrix C");



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
				cu_C.transpose(blasHandle);
			if(alpha!=1.0)
				cu_C *= alpha;
		}
		else
		{
			Faust::MatDense<FPP,Gpu> cu_tmp(cu_B);
			if(typeB == 'T')
				cu_tmp.transpose(blasHandle);
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
				cu_C.transpose(blasHandle);
			if(alpha!=1.0)
				cu_C *= alpha;
		}
		else
		{
			Faust::MatDense<FPP,Gpu> cu_tmp(cu_A);
			if(typeA == 'T')
				cu_tmp.transpose(blasHandle);
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

	const Faust::MatDense<FPP,Gpu>* cu_A_ptr = &cu_A;
     if((cu_A.getDevice() != best_device) || (cu_A.getDevice() == best_device && cu_A == cu_C))
		cu_A_ptr = new Faust::MatDense<FPP,Gpu>(cu_A, best_device);

	const Faust::MatDense<FPP,Gpu>* cu_B_ptr = &cu_B;
	if((cu_B.getDevice() != best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C && cu_B!=cu_A))
		cu_B_ptr = new Faust::MatDense<FPP,Gpu>(cu_B, best_device);


	Faust::MatDense<FPP,Gpu>* cu_C_ptr;
	if(cu_C.getDevice() == best_device)
		cu_C_ptr = &cu_C;
	else
		cu_C_ptr = new Faust::MatDense<FPP,Gpu>(cu_C, best_device);


   if(beta == 0.0)
      cu_C_ptr->resize(opA1,opB2);
   else if(cu_C_ptr->estNulle())
      cu_C_ptr->setZeros(opA1,opB2);
   else if (cu_C_ptr->estIdentite())
   {
      if(opA1 != opB2)
		   handleError("LinAlgebra_gpu", "Faust::gemm : cu_C has been set to identite but is not a square matrix");
      cu_C_ptr->resize(opA1,opB2);
      kernel_memset(cu_C_ptr->getData(), (FPP)0.0, cu_C_ptr->getNbRow()*cu_C_ptr->getNbCol());
      kernel_add_diag_const(cu_C_ptr->getData(), (FPP)1.0, cu_C_ptr->getNbRow());
   }

	if(cu_A_ptr->getData()==NULL || cu_B_ptr->getData()==NULL || cu_C_ptr->getData()==NULL)
		handleError("LinAlgebra_gpu", "Faust::gemm : A, B or C is not allocated");

   faust_cu_gemm(blasHandle.GetHandle(),
       opA, opB, opA1, opB2, opA2,
       &alpha, cu_A_ptr->getData(), cu_A_ptr->getNbRow(),
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

template <typename FPP>
void Faust::setOp(const Faust::MatSparse<FPP,Gpu>& cu_S, const char opA, faust_unsigned_int& opA1, faust_unsigned_int& opA2)
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
        handleError("Faust::setOp","Faust::gemm(const faust_cu_mat&, faust_cu_mat&, ...) : invalid character");
}

template <typename FPP>
void gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const  FPP alpha, const FPP beta, const char opA, Faust::SpBlasHandle<Gpu> spblasHandle)
{
#ifdef __COMPILE_TIMERS__
cu_A.t_csrmv.start();
#endif
   faust_unsigned_int opA1,opA2;

   Faust::setOp(cu_A,opA,opA1,opA2);
   const cusparseOperation_t cusparseOpA = (opA=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;

   if(opA2 != cu_x.size())
        handleError("LinAlgebra_gpu","gemv(const faust_cu_spmat&, const faust_cu_vec&, faust_cu_vec&, ...) : incorrect dimension between A and x");

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

        //handleError("LinAlgebra_cu","Faust::gemv(const faust_cu_vec&, faust_cu_vec&, ...) : vector x has not been initialized");
   else
   {


       // device of cu_A is used to compute
      int currentGPU;
      faust_cudaGetDevice(&currentGPU);
      const int best_device = cu_A.getDevice();
      faust_cudaSetDevice(best_device);

      const Faust::Vect<FPP,Gpu>* cu_x_ptr = &cu_x;
      if((cu_x.getDevice()!=best_device) || (cu_x.getDevice()==best_device && cu_x==cu_y))
         cu_x_ptr = new Faust::Vect<FPP,Gpu>(cu_x, best_device);

      const Faust::MatSparse<FPP,Gpu>* cu_A_ptr = &cu_A;
      if(cu_A.getDevice() != best_device)
         cu_A_ptr = new Faust::MatSparse<FPP,Gpu>(cu_A, best_device);


      Faust::Vect<FPP,Gpu>* cu_y_ptr = &cu_y;
      if (cu_y.getDevice() != best_device)
         if(beta!=0.0)
            cu_y_ptr = new Faust::Vect<FPP,Gpu>(cu_y, best_device);
         else
            cu_y_ptr = new Faust::Vect<FPP,Gpu>(opA1, best_device);
      else if(beta==0.0)
         cu_y_ptr->resize(opA1);


      faust_cu_csrmv(spblasHandle.GetHandle(),cusparseOpA,
                cu_A_ptr->getNbRow(), cu_A_ptr->getNbCol(),
                cu_A_ptr->getNonZeros(),
                &alpha, cu_A_ptr->getDescr(), cu_A_ptr->getValues(),
                cu_A_ptr->getRowPtr(), cu_A_ptr->getColInd(),
                cu_x_ptr->getData(),
                &beta, cu_y_ptr->getData());


      if(cu_y.getDevice() != best_device)
      {
         cu_y = Faust::Vect<FPP,Gpu>(*cu_y_ptr, cu_y.getDevice());
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

template <typename FPP>
void Faust::gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const  FPP alpha, const FPP beta, const char opA, const char opB, Faust::SpBlasHandle<Gpu> spblasHandle)
{

#ifdef __COMPILE_TIMERS__
cu_A.t_csrmm.start();
#endif
    faust_unsigned_int opA1,opA2,opB1,opB2;


    if(opA=='T' && opB=='T')
        handleError("LinAlgebra_cu","Faust::gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : cannot compute transpose(A)*transpose(B). Try doing transpose(B*A) after converting A into a full matrix");

   setOp(cu_A,opA,opA1,opA2);
   setOp(cu_B,opB,opB1,opB2);
   const cusparseOperation_t cusparseOpA = (opA=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;
   const cusparseOperation_t cusparseOpB = (opB=='T')?CUSPARSE_OPERATION_TRANSPOSE:CUSPARSE_OPERATION_NON_TRANSPOSE;




   if(opA2 != opB1)
        handleError("LinAlgebra_cu","Faust::gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : incorrect dimension between A and B");

	if ( (beta!=0) && (cu_C.getNbRow()!=opA1 || cu_C.getNbCol()!=opB2) )
		handleError("LinAlgebra_cu", "Faust::gemm : invalid dimension for output matrix C");



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
            cu_C.init_from_transpose(cu_A, spblasHandle);
         else
            cu_C = cu_A;
         if (alpha!=1.0)
            cu_C *= alpha;
      }
      else
      {
         Faust::MatSparse<FPP,Gpu> cu_S(cu_A, cu_C.getDevice());
         if (opA=='T')
            cu_S.transpose(spblasHandle);
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
      handleError("LinAlgebra_cu","Faust::gemm(const faust_cu_spmat&, const faust_cu_mat&, faust_cu_mat&, ...) : matrix B has not been initialized");


   int currentGPU;
   faust_cudaGetDevice(&currentGPU);
   const int best_device = cu_B.getDevice();
   faust_cudaSetDevice(best_device);

   const Faust::MatDense<FPP,Gpu>* cu_B_ptr = &cu_B;
	if((cu_B.getDevice()!=best_device) || (cu_B.getDevice()==best_device && cu_B==cu_C))
      cu_B_ptr = new Faust::MatDense<FPP,Gpu>(cu_B, best_device);

   const Faust::MatSparse<FPP,Gpu>* cu_A_ptr = &cu_A;
   if(cu_A.getDevice() != best_device)
      cu_A_ptr = new Faust::MatSparse<FPP,Gpu>(cu_A, best_device);

   Faust::MatDense<FPP,Gpu>* cu_C_ptr = &cu_C;
   if  (beta!=0)
   {
      if(cu_C.getDevice() != best_device)
      cu_C_ptr = new Faust::MatDense<FPP,Gpu>(cu_C, best_device);
   }
   else if(cu_C.getDevice() != best_device)
      cu_C_ptr = new Faust::MatDense<FPP,Gpu>(opA1, opB2, best_device);
   else
      cu_C_ptr->resize(opA1, opB2, best_device);



   faust_cu_csrmm(spblasHandle.GetHandle(),cusparseOpA,cusparseOpB,
         cu_A_ptr->getNbRow(), cu_C_ptr->getNbCol(), cu_A_ptr->getNbCol(),
         cu_A_ptr->getNonZeros(),
         &alpha, cu_A_ptr->getDescr(), cu_A_ptr->getValues(),
         cu_A_ptr->getRowPtr(), cu_A_ptr->getColInd(),
         cu_B_ptr->getData(), cu_B_ptr->getNbRow(),
         &beta, cu_C_ptr->getData(), cu_C_ptr->getNbRow());


   if(cu_C_ptr->getDevice()!=cu_C.getDevice())
   {
      cu_C = Faust::MatDense<FPP,Gpu>(*cu_C_ptr, cu_C.getDevice());
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

#endif // __COMPILE_SPMAT__



////// FUNCTIONS with faust_cu_mat matrices //////

// y = op(A) * x
template <typename FPP>
void Faust::gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, char typeA, Faust::BlasHandle<Gpu> blasHandle)
{gemv(cu_A, cu_x, cu_y, FPP(1.0), FPP(0.0), typeA, blasHandle);}
// y = A * x
template <typename FPP>
void Faust::gemv(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle)
{Faust::gemv(cu_A, cu_x, cu_y, (FPP)1.0, (FPP)0.0, 'N', blasHandle);}
// y = A * x
template <typename FPP>
void Faust::multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::BlasHandle<Gpu> blasHandle)
{Faust::gemv(cu_A, cu_x, cu_y, FPP(1.0), FPP(0.0), 'N', blasHandle);}


// C = op(A) * op(B)
template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, char typeA, char typeB, Faust::BlasHandle<Gpu> blasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), typeA, typeB, blasHandle);}
// C = A * B
template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', blasHandle);}
// C = A * B
template <typename FPP>
void Faust::multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', blasHandle);}

//////////////////////////////////////////////////

#ifdef __COMPILE_SPMAT__
///// FUNCTIONS with faust_cu_spmat matrices /////
// y = op(A) * x
template <typename FPP>
void Faust::gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, const char opA, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemv(cu_A, cu_x, cu_y, FPP(1.0), FPP(0.0), opA, spblasHandle);}
// y = A * x
template <typename FPP>
void Faust::gemv(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemv(cu_A, cu_x, cu_y, FPP(1.0), FPP(0.0), 'N', spblasHandle);}
// y = A * x
template <typename FPP>
void Faust::multiply(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::Vect<FPP,Gpu>& cu_x, Faust::Vect<FPP,Gpu>& cu_y, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemv(cu_A, cu_x, cu_y, FPP(1.0), FPP(0.0), 'N', spblasHandle);}

// C = op(A) * op(B) (with A sparse)
template <typename FPP>
void Faust::gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const char opA, const char opB, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), opA, opB, spblasHandle);}
// C = A * B (with A sparse)
template <typename FPP>
void Faust::gemm(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', spblasHandle);}
template <typename FPP>
void Faust::multiply(const Faust::MatSparse<FPP,Gpu>& cu_A, const Faust::MatDense<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::SpBlasHandle<Gpu> spblasHandle)
{ Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', spblasHandle);}

// C = alpha*op(A)*op(B) + beta*C (with B sparse)
template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const FPP alpha, const FPP beta, const char opA, const char opB, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle)
{
	char opAt=(opA=='N')?'T':'N';
	char opBt=(opB=='N')?'T':'N';
	if(opAt=='T' && opBt=='T')
	{
		Faust::MatDense<FPP,Gpu> cu_B_full(cu_B, spblasHandle);
		Faust::gemm(cu_A, cu_B_full, cu_C, alpha, beta, 'N', 'N', blasHandle);
	}
	else
	{
		Faust::gemm(cu_B, cu_A, cu_C, alpha, beta, opBt, opAt, spblasHandle);
		cu_C.transpose(blasHandle);
	}
}
// C = op(A) * op(B) (with B sparse)
template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, const char opA, const char opB, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), opA, opB, blasHandle, spblasHandle);}
// C = A * B (with B sparse)
template <typename FPP>
void Faust::gemm(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle)
{Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', blasHandle, spblasHandle);}
// C = A * B (with B sparse)
template <typename FPP>
void Faust::multiply(const Faust::MatDense<FPP,Gpu>& cu_A, const Faust::MatSparse<FPP,Gpu>& cu_B, Faust::MatDense<FPP,Gpu>& cu_C, Faust::BlasHandle<Gpu> blasHandle, Faust::SpBlasHandle<Gpu> spblasHandle)
{ Faust::gemm(cu_A, cu_B, cu_C, FPP(1.0), FPP(0.0), 'N', 'N', blasHandle, spblasHandle);}

//////////////////////////////////////////////////
#endif // __COMPILE_SPMAT__


#endif
