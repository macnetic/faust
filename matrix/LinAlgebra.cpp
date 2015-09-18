#include <iostream>
#include "LinAlgebra.h"
#include "faust_mat.h"
#include "faust_vec.h"
#include <Eigen/QR>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#include "faust_spmat.h"
#include "faust_exception.h"

	//////////FONCTION faust_mat - faust_mat ////////////////////

#ifdef __COMPILE_TIMERS__
	#include "faust_timer.h"
#endif

#ifdef __GEMM_WITH_OPENBLAS__
	#include "cblas.h"
#endif


const char * interface_name="LinAlgebra";

 void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C)
{   
#ifdef __COMPILE_TIMERS__
A.t_multiply.start();
#endif

	if (A.getNbCol() != B.getNbRow())
	{
		//handleError("Linalgebra : multiply :  nbCol of A = %d while nbRow of B = %d",A.getNbCol(),B.getNbRow());
		handleError(interface_name,"multiply : invalid matrix dimensions");	
	}
	 
	if(A.isZeros || B.isZeros)
	{
		C.resize(A.dim1,B.dim2);
		faust_real *const ptr_data_dst = C.getData();
		memset(ptr_data_dst, 0, sizeof(faust_real) * C.dim1*C.dim2);
		C.isZeros = true;
		C.isIdentity = false;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}

	if(B.isIdentity)
	{
		C=A;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}

	if(A.isIdentity)
	{
		C=B;
		#ifdef __COMPILE_TIMERS__
			A.t_multiply.stop();
		#endif
		return;
	}


	if (((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))))
	{
		handleError(interface_name," multiply : C is the same object as A or B");		
	}else
	{
		C.resize(A.getNbRow(),B.getNbCol());
		C.mat.noalias() = A.mat * B.mat;
	}
	C.isZeros = false;
	C.isIdentity = false;

#ifdef __COMPILE_TIMERS__
A.t_multiply.stop();
#endif
}


faust_vec solve(const faust_mat & A, const faust_vec & v)
{
	if (A.getNbRow() != v.getDim())
	{
		handleError(interface_name," : solve : number of row of the  dense matrix A is different from dimension of the vector v");	
	}
	faust_vec sol(A.getNbCol());
	sol.Display();
	v.Display();
	A.Display();
	sol.vec=A.mat.colPivHouseholderQr().solve(v.vec);	
	return sol;
}

void solve(const faust_spmat & A,faust_vec & x, const faust_vec & y)
{
	if (A.getNbRow() != y.getDim())
	{
		handleError(interface_name,"solve : number of row of the sparse matrix A is different from dimension of the vector y");	
	}
	x.resize(A.getNbCol());
	x.resize(A.getNbCol());
	Eigen::SparseQR<Eigen::SparseMatrix<faust_real>, Eigen::COLAMDOrdering<int> >   solver(A.mat);
	x.vec=solver.solve(y.vec);
}



void gemv(const faust_mat & A,const faust_vec & x,faust_vec & y,const faust_real & alpha, const faust_real & beta, char typeA)
{
	faust_unsigned_int nbRowOpA,nbColOpA;
	const faust_vec* px;
	if  ((&x) == (&y))
		px = new faust_vec(x);
	else
		px = &x;
		
	
	if (typeA == 'T')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{	
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}
	
	if   (nbColOpA != px->getDim() )
	{
		//handleError("Linalgebra : gemv : nbCol of op(A) = %d while dim of x = %d",nbColOpA,px->getDim());
		handleError(interface_name, "gemv : dimension conflict  between matrix op(A) and input vector x");
	}
	
	if ( (beta!=0)  &&  (y.getDim() != nbRowOpA))
	{
		handleError(interface_name, "gemv : dimension conflict  between matrix op(A) and output vector y");
	}
	
	y.resize(nbRowOpA);
	
	
	
	#ifdef __GEMM_WITH_OPENBLAS__
		CBLAS_TRANSPOSE transA,transB;
		if (typeA=='T')
			transA = CblasTrans;
		else 
			transA = CblasNoTrans;	
	#endif
	
	#ifndef __GEMM_WITH_OPENBLAS__
	if (beta == 0.0)
	{	
		if (typeA == 'N')
		{
			y.vec.noalias() = alpha * A.mat * px->vec;			
		}else
		{
		
			y.vec.noalias() = alpha * A.mat.transpose() * px->vec;
		}
	}else
	{	
		if (typeA == 'N')
		{
			y.vec = alpha * A.mat * px->vec + beta * y.vec;			
		}else
		{
			y.vec = alpha * A.mat.transpose() * px->vec + beta * y.vec;
		}
	}
	#else
		#ifdef FAUST_SINGLE	
				cblas_sgemv(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);
		#else
			cblas_dgemv(CblasColMajor,transA,A.getNbRow(),A.getNbCol(),alpha,A.getData(),A.getNbRow(),px->getData(),1,beta,y.getData(),1);
		#endif
	#endif
							
	if  ((&x) == (&y))
		delete px; 
	px=NULL;
	
}
	
	
	
void gemm(const faust_mat & A,const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA, char  typeB)
{
#ifdef __COMPILE_TIMERS__
A.t_gemm.start();
#endif
	faust_unsigned_int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{
		handleError(interface_name, "gemv : gemm : C is the same object as A or B");		
	}

	if (typeA == 'T')
	{
		nbRowOpA = A.getNbCol();
		nbColOpA = A.getNbRow();
	}else
	{
		nbRowOpA = A.getNbRow();
		nbColOpA = A.getNbCol();
	}


	if (typeB == 'T')
	{
		nbRowOpB = B.getNbCol();
		nbColOpB = B.getNbRow();
	}else
	{
		nbRowOpB = B.getNbRow();
		nbColOpB = B.getNbCol();
	}


	if (nbColOpA != nbRowOpB)
	{
		handleError(interface_name, "gemm : dimension conflict  between matrix op(A) and matrix op(B)");
		
	}







	if ( (beta!=0)  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{

		//handleError("Linalgebra : gemm : nbRow of op(A) = %d while nbRow of op(C) = %d\n or nbCol of op(B) = %d  while nbCol of C = %d",nbRowOpA,C.getNbRow(),nbColOpB,C.getNbCol());
		handleError(interface_name, "gemm : invalid dimension for output matrix C");
	
	}
		
        C.resize(nbRowOpA,nbColOpB);


	#ifdef __GEMM_WITH_OPENBLAS__
		CBLAS_TRANSPOSE transA,transB;
		if (typeA=='T')
			transA = CblasTrans;
		else 
			transA = CblasNoTrans;
		if (typeB=='T')
			transB = CblasTrans;
		else 
			transB = CblasNoTrans;
	#endif



	if (beta == 0.0)
	{

		if(A.isZeros || B.isZeros)
		{
			
			faust_real *const ptr_data_dst = C.getData();
			memset(ptr_data_dst, 0, sizeof(faust_real) * C.dim1*C.dim2);
			C.isZeros = true;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
	
		if(A.isIdentity)
		{
			C=B;
			if(typeB == 'T') 
				C.transpose();
			if(alpha!=1.0) 
				C*= alpha;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
		if(B.isIdentity)
		{
			C=A;
			if(typeA == 'T') 
				C.transpose();
			if(alpha!=1.0) 
				C*= alpha;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}


		#ifndef __GEMM_WITH_OPENBLAS__
			if (typeA == 'N')
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat * B.mat;			
				else
					C.mat.noalias() = alpha * A.mat * B.mat.transpose();
			}else
			{
				if (typeB == 'N')
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat;			
				else
					C.mat.noalias() = alpha * A.mat.transpose() * B.mat.transpose();
			}
		#else
			#ifdef FAUST_SINGLE
				cblas_sgemm(CblasColMajor, transA, transB, C.dim1, C.dim2, nbColOpA, alpha, A.getData(), A.dim1, B.getData(), B.dim1, 0.0f, C.getData(), C.dim1);
			#else
				cblas_dgemm(CblasColMajor, transA, transB, C.dim1, C.dim2, nbColOpA, alpha, A.getData(), A.dim1, B.getData(), B.dim1, 0.0, C.getData(), C.dim1);
			#endif
		#endif


	}else
	{
		if(A.isZeros || B.isZeros)
		{
			C *= beta;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
	
		if(A.isIdentity)
		{
			C *= beta;
			if(typeB == 'N' && alpha == 1.0)
			{
				C += B;
				#ifdef __COMPILE_TIMERS__
					A.t_gemm.stop();
				#endif
				return;
			}
			faust_mat B_tmp(B);
			if(typeB == 'T') 
				B_tmp.transpose();
			if(alpha != 1.0) 
				B_tmp *= alpha;
			C += B_tmp;
			C.isZeros = false;
			C.isIdentity = false;
			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}
		if(B.isIdentity)
		{
			C *= beta;
			if(typeA == 'N' && alpha == 1.0)
			{
				C += A;
				C.isZeros = false;
				C.isIdentity = false;
				#ifdef __COMPILE_TIMERS__
					A.t_gemm.stop();
				#endif
				return;
			}
			faust_mat A_tmp(A);
			if(typeA == 'T') 
				A_tmp.transpose();
			if(alpha != 1.0) 
				A_tmp *= alpha;
			C += A_tmp;
			C.isZeros = false;
			C.isIdentity = false;

			#ifdef __COMPILE_TIMERS__
				A.t_gemm.stop();
			#endif
			return;
		}


		#ifndef __GEMM_WITH_OPENBLAS__
			if (typeA == 'N')
			{
				if (typeB == 'N')
						C.mat = alpha * A.mat * B.mat + beta * C.mat;
				else
					C.mat = alpha * A.mat * B.mat.transpose() + beta * C.mat;
			}else
			{
				if (typeB == 'N')
					C.mat = alpha * A.mat.transpose() * B.mat + beta * C.mat ;			
				else
					C.mat = alpha * A.mat.transpose() * B.mat.transpose() + beta * C.mat;
			}
		#else
                
			#ifdef FAUST_SINGLE
				cblas_sgemm(CblasColMajor, transA, transB, C.dim1, C.dim2, nbColOpA, alpha, A.getData(), A.dim1, B.getData(), B.dim1, beta, C.getData(), C.dim1);
			#else
				cblas_dgemm(CblasColMajor, transA, transB, C.dim1, C.dim2, nbColOpA, alpha, A.getData(), A.dim1, B.getData(), B.dim1, beta, C.getData(), C.dim1);
			#endif
		#endif

	}
	C.isZeros = false;
	C.isIdentity = false;
#ifdef __COMPILE_TIMERS__
A.t_gemm.stop();
#endif
}




void add(const faust_mat & A, const faust_mat & B, faust_mat & C)
{   
#ifdef __COMPILE_TIMERS__
A.t_add_ext.start();
#endif
	if ((A.getNbCol() != B.getNbCol()) || (A.getNbRow() != B.getNbRow()) || (A.getNbRow() != C.getNbRow()) || (A.getNbCol() != C.getNbCol()))
	{
		handleError(interface_name," add : matrix dimension not equal");
		
	}else
	{

			C.mat = A.mat + B.mat;

	}
	C.isZeros = false;
	C.isIdentity = false;
#ifdef __COMPILE_TIMERS__
A.t_add_ext.stop();
#endif
}


	
	
// compute the biggest eigenvalue of A, A must be semi-definite positive 	
faust_real power_iteration(const  faust_mat & A, const faust_unsigned_int nbr_iter_max,faust_real threshold, faust_int & flag)
{	
	#ifdef __COMPILE_TIMERS__
		A.t_power_iteration.start();
	#endif 	
	 faust_unsigned_int nb_col = A.getNbCol();
	 faust_unsigned_int nb_row = A.getNbRow();
	 faust_unsigned_int i = 0;
	 faust_unsigned_int k;
	 bool do_continue=true;
	 faust_real abs_eigen_value;
	 
	 bool stop_crit;
	 flag = 0;

	 
	 if (nbr_iter_max <= 0)
	 {
		handleError(interface_name,"power_iteration : nbr_iter_max <= 0");
	 }
	 if (nb_col != nb_row)
	 {
		handleError(interface_name,"power_iteration : A must be square");
        
	 }
	 

	 
	 faust_vec xk(nb_col);
	 faust_vec xk_pp(nb_col);
	 xk.setOnes();
	 xk.normalize();
	 
	 while(do_continue)
	 {	
		i++;
		gemv(A,xk,xk_pp,1.,0,'N');
		abs_eigen_value = xk_pp.norm();
		//std::cout<<"current_norm : "<< abs_eigen_value<<std::endl;
		xk_pp.scalarMultiply(1/abs_eigen_value);
		
		stop_crit =true;
		k=0;
		while ((stop_crit) && (k<nb_col))
		{
			if (fabs(xk_pp(k) - xk(k))>threshold)
			{
				stop_crit = false;
			}
			k++;
		}
		
		if (stop_crit)
		{
			do_continue = false;
			flag = i;
			//std::cout<<"convergence : "<< i <<std::endl;

		}
		
		if (i >= nbr_iter_max)
		{	//std::cout<<"divergence"<<std::endl;
			do_continue = false;
			flag = -1;
		}
		
		
		xk = xk_pp;
	 }
	 //std::cout<<" flag :"<<flag<<std::endl;
	 #ifdef __COMPILE_TIMERS__
		A.t_power_iteration.stop();
	#endif
	/*std::cout<<"flag inside power_it : "<<flag<<std::endl;
	std::cout<<"threshold inside power_it : "<<threshold<<std::endl;
	std::cout<<"max_it inside power_it : "<<nbr_iter_max<<std::endl;*/		
	 return abs_eigen_value;
	 
}
































	






// non-member operators definitions

faust_vec operator*(const faust_mat& M, const faust_vec& v)
{
	faust_vec vec(v);
	vec.multiplyLeft(M);
	return vec;
}

	

