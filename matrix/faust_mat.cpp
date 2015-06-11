#include "faust_mat.h"
//#include <cmath>
#include <iostream>
using namespace std;



/// CONSTRUCTEUR ///
  faust_mat::faust_mat(const Eigen::Matrix<faust_real,Eigen::Dynamic,Eigen::Dynamic> & mat_) : 
     mat(mat_), dim1(mat_.rows()), dim2(mat_.cols()){}
  
  faust_mat::faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol )
  {
	  int i,j;
	  if ((nbRow < 0) || (nbCol < 0))
	  {
		cerr << "ERREUR une dimension ne peut etre negative" << endl;; 
        	exit( EXIT_FAILURE);  
	  }
	  
	dim1 = nbRow;
	dim2 = nbCol;
	Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> m(nbRow,nbCol);
		
	for (int j=0;j<nbCol;j++)
		for (int i=0;i<nbRow;i++)
			m(i,j) = mat_[i+(j*nbRow)];
	mat = m;
	  
  }
 

/// GETTEUR SETTEUR ///


faust_real faust_mat::getCoeff(const int i,const int j) const
{
	if  ( (i<0) || (i>=dim1) )
	{
		cerr << "ERROR getCoeff : index out of range : 0>i  or i>=nb_row" << endl;
		exit( EXIT_FAILURE);		
	}
	
	if  ( (j<0) || (j>=dim2) )
	{
		cerr << "ERROR getCoeff : index out of range : 0>j  or j>=nb_col" << endl;
		exit( EXIT_FAILURE);		
	}
	
	return mat(i,j);
}



void faust_mat::getCoeffs(std::vector<faust_real> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const
{
	if ( (id_row.size()!= id_col.size()) || (id_col.size() != valueS.size()) )
	{
		cerr << "ERROR getCoeffs : id_row, id_col, valueS don't have the same length" << endl;
		exit( EXIT_FAILURE);
	}
	
	int n=id_col.size();
	int current_id_row,current_id_col;

	for (int i=0;i<n;i++)
	{
		current_id_row =  id_row[i];
		current_id_col =  id_col[i];
		valueS[i]=getCoeff(current_id_row,current_id_col);
	}
}

void faust_mat::setCoeff(const faust_real & value,const int id_row, const int id_col)
{
	if ( (id_row < 0) || (id_row >= dim1) )
	{
		cerr << "ERROR setCoeff : index out of range : 0>id_row  or id_row>=nb_row" << endl;
		exit( EXIT_FAILURE);	
	}
			
	if ( (id_col < 0) || (id_col >= dim2) )
	{
		cerr << "ERROR setCoeff : index out of range : (0>id_col)  or (id_col >= nb_col)" << endl;
		exit( EXIT_FAILURE);	
	}
			
	mat(id_row,id_col) = value;
}


void faust_mat::setCoeffs(const faust_real value,const std::vector<int> & id_row,const std::vector<int>  & id_col)
{
	if (id_row.size()!= id_col.size())
	{
		cerr << "ERREUR setCoeffs : id_row and id_col don't have the same length" << endl;
		exit( EXIT_FAILURE);
	}
	
	int n = id_row.size();
	for (int i=0;i<n;i++)
		setCoeff(value,id_row[i],id_col[i]);
}


void faust_mat::setCoeffs(const std::vector<faust_real> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col)
{
	if ( (id_row.size()!= id_col.size()) || (id_row.size()!= valueS.size()) )
	{
		cerr << "ERREUR setCoeffs : id_row,id_col,valueS don't have the same length" << endl;
		exit( EXIT_FAILURE);
	}
	
	int n = id_row.size();
	for (int i=0;i<n;i++)
		setCoeff(valueS[i],id_row[i],id_col[i]);
	
}
	

 
 void faust_mat::resize(const int nbRow,const int nbCol)
{	
		if ( (nbRow <0) || (nbCol <0) )
		{
			cerr << "ERREUR resize : les nouvelles dimensions doivent etre strictement positive" << endl;
			exit( EXIT_FAILURE);
		}
		else if ((dim1 != nbRow) || (dim2 != nbCol))
		{
			dim1 = nbRow;
			dim2 = nbCol;
			mat.resize(nbRow,nbCol);
		}
}
 
 
 /// EGALITE ///

 bool faust_mat::isEqual(const faust_mat & B) const 
 {
	if ((getNbCol() != B.getNbCol()) || (getNbRow() != B.getNbRow()))
	{
		cerr << "ERREUR isEqual : dimension of the matrix are not the same" << endl; 
        	exit( EXIT_FAILURE);	
	}

	if (isZeros())
		return B.isZeros();
	else if (B.isZeros())
		return isZeros();
	else
		return mat.isApprox(B.mat,FAUST_PRECISION);
 }


 
 /// OPERATION BASIQUE ///
 
//arithmetique
faust_real faust_mat::max(std::vector<int> & id_row,std::vector<int> & id_col) const
{
	faust_real maxi = max();
	int i,j,k;
	if ( (id_row.size() != 0) || (id_col.size() != 0) )
	{
		cerr << "ERREUR max : sizes of id_row and id_col must be equal to zero" << endl;
		exit( EXIT_FAILURE);	
	}
	
	for (j=0;j<getNbCol();j++)
		for (i=0;i<getNbRow();i++)
			if (mat(i,j) == maxi)
			{
				id_row.push_back(i);
				id_col.push_back(j);
			}
	return maxi;
}
 
 
 
 void faust_mat::transpose()
 {
	/*Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat_copy = mat; 
	int dim1_copy = dim1;
	int dim2_copy = dim2;
	resize(dim2_copy,dim1_copy);
	 
	 mat = mat_copy.transpose();*/
	mat = mat.transpose().eval();
	int dim1_copy = dim1;
	dim1 = dim2;
        dim2 = dim1_copy; 
        
	//resize(dim2_copy,dim1_copy);
 }
 
 void faust_mat::multiplyRight(faust_mat const& A)
 {  
	if (dim2 != A.dim1)
	{
		std::cerr << "ERREUR multiply : nbCol of this = " << getNbCol(); 
       		std::cerr <<" while nbRow of A = " << A.getNbRow() << std::endl;
        	exit( EXIT_FAILURE);	
	}
	
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
	resize(dim1, A.dim2);
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
	}
	mat = mat + A.mat;
 }

 
 void faust_mat::sub(faust_mat const& A)
 {  
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		std::cerr << "ERREUR sub : matrix dimension not equal" << std::endl; 
        	exit( EXIT_FAILURE);
	}
	mat = mat - A.mat;
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

 
