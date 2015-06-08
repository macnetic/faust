#include "faust_mat.h"
//#include <cmath>
#include <iostream>
using namespace std;



/// CONSTRUCTEUR ///
  faust_mat::faust_mat(const Eigen::Matrix<faust_real,Eigen::Dynamic,Eigen::Dynamic> & mat_) : mat(mat_), dim1(mat_.rows()), dim2(mat_.cols())
  {
	  
  }
  
  faust_mat::faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol )
  {
	  int i,j;
	  if ((nbRow < 0) || (nbCol < 0))
	  {
		cerr << "ERREUR une dimension ne peut etre negative" << endl;; 
        exit( EXIT_FAILURE);  
	  }else
	  {
		dim1 = nbRow;
		dim2 = nbCol;
		Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> m(nbRow,nbCol);
		
		for (int j=0;j<nbCol;j++)
		{
			for (int i=0;i<nbRow;i++)
			{
				m(i,j) = mat_[i+(j*nbRow)];
			}
		}
		
		mat = m;
	  }
  }
  faust_mat::faust_mat(const int nbRow, const int nbCol) : mat(nbRow,nbCol),dim1(nbRow),dim2(nbCol)
  {}
  
  
  
  

  
 
/// GETTEUR SETTEUR ///


 
 void faust_mat::resize(const int nbRow,const int nbCol)
{	
		if ( (nbRow <0) || (nbCol <0) )
		{
			cerr << "ERREUR resize : les nouvelles dimensions doivent etre strictement positive" << endl;
			exit( EXIT_FAILURE);
		}
		else
		{
			if ((dim1 != nbRow) || (dim2 != nbCol))
			{	
				dim1 = nbRow;
				dim2 = nbCol;
				mat.resize(nbRow,nbCol);
			}
		}
	
}







 
 

 

 
 
 /// EGALITE ///

 
 
 bool faust_mat::isEqual(const faust_mat & B) const 
 {
	if ((getNbCol() != B.getNbCol()) || (getNbRow() != B.getNbRow()))
	{
		cerr << "ERREUR isEqual : dimension of the matrix are not the same" << endl; 
        exit( EXIT_FAILURE);	
	
	}else
	{
		if (isZeros())
		{
			return B.isZeros();
		}else
		{
			if (B.isZeros())
			{
				return isZeros();
			}else
			{
			return mat.isApprox(B.mat,FAUST_PRECISION);
			}
		
		}
	}
		
 }



 
 /// OPERATION BASIQUE ///
 
 // frobenius norm

 
 
 
 void faust_mat::transpose()
 {
	/*Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat_copy = mat; 
	int dim1_copy = dim1;
	int dim2_copy = dim2;
	resize(dim2_copy,dim1_copy);
	 
	 mat = mat_copy.transpose();*/
	mat.transpose();
	int dim1_copy = dim1;
	int dim2_copy = dim2;
	resize(dim2_copy,dim1_copy);
 }
 
 void faust_mat::multiplyRight(faust_mat const& A)
 {  if (dim2 != A.dim1)
	{
		std::cerr << "ERREUR multiply : nbCol of this = " << getNbCol(); 
       	std::cerr <<" while nbRow of A = " << A.getNbRow() << std::endl;
        exit( EXIT_FAILURE);	
		 
	}else
	{
		/*int dim1_copy = dim1;
		Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat_copy = mat; 
		resize(dim1_copy,A.dim2);
		if (&(A.mat) != &(mat))
		{	
			mat.noalias() = mat_copy * A.mat;
		}else
		{
			mat = mat_copy * A.mat;
		}*/
		
		mat = mat * A.mat;
		int dim1_copy = dim1;
		resize(dim1_copy,A.dim2);
	}
 }
 
 void faust_mat::scalarMultiply(faust_real const lambda)
 {
	 mat = lambda * mat;
 }
 
 
 
 void faust_mat::add(faust_mat const& A)
 {  
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		std::cerr << "ERREUR add : matrix dimension not equal" << std::endl; 
        exit( EXIT_FAILURE);
		 
	}else
	{ 
		mat = mat + A.mat;

	}
 }



 
  // Affichage
  void faust_mat::Display() const
  {     std::cout << "nb_row=" << getNbRow() << endl;
        std::cout << "nb_col=" << getNbCol()   <<endl;  
		std::cout << mat <<endl; 
  }

  
    /// SURCHARGE OPERATEUR ///
  // affectation
  void faust_mat::operator=(faust_mat const& A)
  {
	  mat = A.mat;
	  dim1 = A.dim1;
	  dim2 = A.dim2;
  }
  
  



bool operator==(faust_mat const& A, faust_mat const& B)
{
	return A.isEqual(B);
}



 
