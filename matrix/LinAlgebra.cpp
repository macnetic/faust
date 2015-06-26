#include <iostream>
#include "LinAlgebra.h"
//////////FONCTION faust_mat - faust_mat ////////////////////

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif



 void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C)
{   
#ifdef __COMPILE_TIMERS__
A.t_multiply.start();
#endif

	if (A.getNbCol() != B.getNbRow())
	{
		std::cerr << "ERREUR multiply : nbCol of A = " << A.getNbCol(); 
       		std::cerr <<" while nbRow of B = " << B.getNbRow() << std::endl;
        	exit( EXIT_FAILURE);		
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
		std::cerr << "ERREUR multiply : C pointe vers A ou B" << std::endl; 
		exit( EXIT_FAILURE);	
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


/*
void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta)
{
	if (A.getNbCol() != B.getNbRow())
	{
		std::cerr << "ERREUR gemm : nbCol of A = " << A.getNbCol(); 
       	std::cerr <<" while nbRow of B = " << B.getNbRow() << std::endl;
        exit( EXIT_FAILURE);
	}
	if ( (C.getNbRow() != A.getNbRow())	|| (C.getNbCol() != B.getNbCol()) )
	{
		std::cerr << "ERREUR gemm : nbRow of A = "<< A.getNbRow() << " while nbRow of C = " << C.getNbRow() << std::endl;
		std::cerr << "or nbCol of B = "<< B.getNbCol() << " while nbCol of C = " << C.getNbCol() << std::endl;
		exit( EXIT_FAILURE);
	}
	 
	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{
		std::cerr << "ERREUR multiply : C pointe vers A ou B" << std::endl; 
		exit( EXIT_FAILURE);	
	}
		
	if ((alpha == 0) && (beta == 0))
	{
		C.setZeros();
	}else
				{
					
					if (beta == 0)
					{
						C.mat.noalias() = alpha * A.mat * B.mat;
					}else
					{
						if (alpha == 0)
						{
							C.mat = beta * C.mat;
						}else
						{
							C.mat = alpha * A.mat * B.mat + beta * C.mat;
						}
					}
				}	
			}
		}

	}
}
*/






/*
void gemm(faust_mat & A, faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA,char  typeB)
{	
	if ( (typeA != 'N') && (typeA != 'T') )
	{
		std::cerr << "ERREUR gemm : typeA different de 'N' et 'T' " << std::endl; 
        exit( EXIT_FAILURE);	
	}
	
	if  ( (typeB != 'N') && (typeB != 'T') )
	{
		std::cerr << "ERREUR gemm : typeB different de 'N' et 'T' " << std::endl; 
        exit( EXIT_FAILURE);
	}

	if (typeA == 'T')
	{
		A.transpose();		
	}
	if (typeB == 'T')
	{
		B.transpose();		
	}
	
	
	if (A.getNbCol() != B.getNbRow())
	{
		std::cerr << "ERREUR gemm : nbCol of op(A) = " << A.getNbCol(); 
       	std::cerr <<" while nbRow of op(B) = " << B.getNbRow() << std::endl;
        exit( EXIT_FAILURE);
	}
	
	
	if ( (beta!=0)  && ( (C.getNbRow() != A.getNbRow())	|| (C.getNbCol() != B.getNbCol()) ) )
	{
		std::cerr << "ERREUR gemm : nbRow of op(A) = "<< A.getNbRow() << " while nbRow of C = " << C.getNbRow() << std::endl;
		std::cerr << "or nbCol of op(B) = "<< B.getNbCol() << " while nbCol of C = " << C.getNbCol() << std::endl;
		exit( EXIT_FAILURE);
	
	}else
	{
		C.resize(A.getNbRow(),B.getNbCol());
	}
	
	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{
		std::cerr << "ERREUR gemm : C pointe vers A ou B" << std::endl; 
		exit( EXIT_FAILURE);	
	}
		
	if ((alpha == 0) && (beta == 0))
	{
		C.setZeros();
	}else
	{
					
		if (beta == 0)
		{
			C.mat.noalias() = alpha * A.mat * B.mat;
		}else
		{
			if (alpha == 0)
			{
				C.mat = beta * C.mat;
			}else
			{
				C.mat = alpha * A.mat * B.mat + beta * C.mat;
			}
		}
	}	

	
	if (typeA == 'T')
	{
		A.transpose();		
	}
	if (typeB == 'T')
	{
		B.transpose();		
	}
	
}*/


void gemm(const faust_mat & A,const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA, char  typeB)
{
#ifdef __COMPILE_TIMERS__
A.t_gemm.start();
#endif
	int nbRowOpA,nbRowOpB,nbColOpA,nbColOpB;

	if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
	{
		std::cerr << "ERREUR newgemm : C pointe vers A ou B" << std::endl; 
		exit( EXIT_FAILURE);	
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
		std::cerr << "ERREUR gemm : nbCol of op(A) = " << nbColOpA; 
		std::cerr <<" while nbRow of op(B) = " << nbColOpB << std::endl;
		exit( EXIT_FAILURE);
	}







	if ( (beta!=0)  && ( (C.getNbRow() != nbRowOpA)	|| (C.getNbCol() != nbColOpB) ) )
	{
		std::cerr << "ERREUR gemm : nbRow of op(A) = "<< A.getNbRow() << " while nbRow of C = " << C.getNbRow() << std::endl;
		std::cerr << "or nbCol of op(B) = "<< B.getNbCol() << " while nbCol of C = " << C.getNbCol() << std::endl;
	exit( EXIT_FAILURE);
	}
		
        C.resize(nbRowOpA,nbColOpB);

	if (beta == 0)
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
		std::cerr << "ERREUR add : matrix dimension not equal" << std::endl; 
        exit( EXIT_FAILURE);		
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



