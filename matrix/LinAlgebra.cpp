#include <iostream>
#include "LinAlgebra.h"
//////////FONCTION faust_mat - faust_mat ////////////////////


 void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C)
{   
	if (A.getNbCol() != B.getNbRow())
	{
		std::cerr << "ERREUR multiply : nbCol of A = " << A.getNbCol(); 
       	std::cerr <<" while nbRow of B = " << B.getNbRow() << std::endl;
        exit( EXIT_FAILURE);		
	}else
	{ 
		if (((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))))
		{
			std::cerr << "ERREUR multiply : C pointe vers A ou B" << std::endl; 
			exit( EXIT_FAILURE);	
		}else
		{
			C.resize(A.getNbRow(),B.getNbCol());
			C.mat.noalias() = A.mat * B.mat;
		}

	}
}


void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta)
{
	if (A.getNbCol() != B.getNbRow())
	{
		std::cerr << "ERREUR gemm : nbCol of A = " << A.getNbCol(); 
       	std::cerr <<" while nbRow of B = " << B.getNbRow() << std::endl;
        exit( EXIT_FAILURE);
	}else
	{
		if ( (C.getNbRow() != A.getNbRow())	|| (C.getNbCol() != B.getNbCol()) )
		{
			std::cerr << "ERREUR gemm : nbRow of A = "<< A.getNbRow() << " while nbRow of C = " << C.getNbRow() << std::endl;
			std::cerr << "or nbCol of B = "<< B.getNbCol() << " while nbCol of C = " << C.getNbCol() << std::endl;
			exit( EXIT_FAILURE);
		}else
		{ 
			if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
			{
				std::cerr << "ERREUR multiply : C pointe vers A ou B" << std::endl; 
				exit( EXIT_FAILURE);	
			}else
			{	
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
	}else
	{
		if ( (C.getNbRow() != A.getNbRow())	|| (C.getNbCol() != B.getNbCol()) )
		{
			std::cerr << "ERREUR gemm : nbRow of op(A) = "<< A.getNbRow() << " while nbRow of C = " << C.getNbRow() << std::endl;
			std::cerr << "or nbCol of op(B) = "<< B.getNbCol() << " while nbCol of C = " << C.getNbCol() << std::endl;
			exit( EXIT_FAILURE);
		}else
		{ 
			if ( ((&(C.mat)) == (&(A.mat))) || ((&(C.mat)) == (&(B.mat))) )
			{
				std::cerr << "ERREUR multiply : C pointe vers A ou B" << std::endl; 
				exit( EXIT_FAILURE);	
			}else
			{	
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
	
	if (typeA == 'T')
	{
		A.transpose();		
	}
	if (typeB == 'T')
	{
		B.transpose();		
	}
	
	
	
	
	
}

	
	
	
void add(const faust_mat & A, const faust_mat & B, faust_mat & C)
{   
	if ((A.getNbCol() != B.getNbCol()) || (A.getNbRow() != B.getNbRow()) || (A.getNbRow() != C.getNbRow()) || (A.getNbCol() != C.getNbCol()))
	{
		std::cerr << "ERREUR add : matrix dimension not equal" << std::endl; 
        exit( EXIT_FAILURE);		
	}else
	{

			C.mat = A.mat + B.mat;

	}
}



