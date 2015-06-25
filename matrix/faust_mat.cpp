#include "faust_mat.h"
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;



/// CONSTRUCTEUR ///
  faust_mat::faust_mat(const Eigen::Matrix<faust_real,Eigen::Dynamic,Eigen::Dynamic> & mat_) : 
     mat(mat_), dim1(mat_.rows()), dim2(mat_.cols()){}
  
  faust_mat::faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol )
  {

#ifdef __COMPILE_TIMERS__
t_constr.start();
#endif

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
	  
#ifdef __COMPILE_TIMERS__
t_constr.stop();
#endif

  }
 

/// GETTEUR SETTEUR ///


faust_real faust_mat::getCoeff(const int i,const int j) const
{

#ifdef __COMPILE_TIMERS__
t_get_coeff.start();
#endif

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
	

#ifdef __COMPILE_TIMERS__
t_get_coeff.stop();
#endif

	return mat(i,j);
}



void faust_mat::getCoeffs(std::vector<faust_real> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const
{

#ifdef __COMPILE_TIMERS__
t_get_coeffs.start();
#endif

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

#ifdef __COMPILE_TIMERS__
t_get_coeffs.stop();
#endif

}

void faust_mat::setCoeff(const faust_real & value,const int id_row, const int id_col)
{

#ifdef __COMPILE_TIMERS__
t_set_coeff.start();
#endif

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

#ifdef __COMPILE_TIMERS__
t_set_coeff.stop();
#endif

}


void faust_mat::setCoeffs(const faust_real value,const std::vector<int> & id_row,const std::vector<int>  & id_col)
{
#ifdef __COMPILE_TIMERS__
t_set_coeffs.start();
#endif
	if (id_row.size()!= id_col.size())
	{
		cerr << "ERREUR setCoeffs : id_row and id_col don't have the same length" << endl;
		exit( EXIT_FAILURE);
	}
	
	int n = id_row.size();
	for (int i=0;i<n;i++)
		setCoeff(value,id_row[i],id_col[i]);
#ifdef __COMPILE_TIMERS__
t_set_coeffs.stop();
#endif
}


void faust_mat::setCoeffs(const std::vector<faust_real> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col)
{
#ifdef __COMPILE_TIMERS__
t_set_coeffs2.start();
#endif
	if ( (id_row.size()!= id_col.size()) || (id_row.size()!= valueS.size()) )
	{
		cerr << "ERREUR setCoeffs : id_row,id_col,valueS don't have the same length" << endl;
		exit( EXIT_FAILURE);
	}
	
	int n = id_row.size();
	for (int i=0;i<n;i++)
		setCoeff(valueS[i],id_row[i],id_col[i]);
	
#ifdef __COMPILE_TIMERS__
t_set_coeffs2.stop();
#endif
}
	

 
 void faust_mat::resize(const int nbRow,const int nbCol)
{	
#ifdef __COMPILE_TIMERS__
t_resize.start();
#endif
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
#ifdef __COMPILE_TIMERS__
t_resize.stop();
#endif
}

 void faust_mat::check_dim_validity()
 {

#ifdef __COMPILE_TIMERS__
t_check_dim.start();
#endif

	bool verifSize = (getNbCol() == mat.cols()) &&  (getNbRow() == mat.rows());
	
	if (!verifSize)
	{
		std::cerr << "Error in faust_mat::check_dim_validity : Size incompatibility in the faust_mat" << std::endl;
		exit(EXIT_FAILURE);
	}
#ifdef __COMPILE_TIMERS__
t_check_dim.stop();
#endif

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

#ifdef __COMPILE_TIMERS__
t_max.start();
#endif

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

#ifdef __COMPILE_TIMERS__
t_max.stop();
#endif

	return maxi;
}

faust_real faust_mat::min(std::vector<int> & id_row,std::vector<int> & id_col) const
{
	faust_real mini = min();
	int i,j,k;
	if ( (id_row.size() != 0) || (id_col.size() != 0) )
	{
		cerr << "ERREUR max : sizes of id_row and id_col must be equal to zero" << endl;
		exit( EXIT_FAILURE);	
	}
	
	for (j=0;j<getNbCol();j++)
		for (i=0;i<getNbRow();i++)
			if (mat(i,j) == mini)
			{
				id_row.push_back(i);
				id_col.push_back(j);
			}
	return mini;
}




 
 
 
 void faust_mat::transpose()
 {

#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif

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

#ifdef __COMPILE_TIMERS__
t_transpose.stop();
#endif

 }
 
 void faust_mat::multiplyRight(faust_mat const& A)
 {  

#ifdef __COMPILE_TIMERS__
t_mult_right.start();
#endif

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

#ifdef __COMPILE_TIMERS__
t_mult_right.stop();
#endif

 }
 
 void faust_mat::multiplyLeft(faust_mat const& A)
 { 

#ifdef __COMPILE_TIMERS__
 t_mult_left.start();
#endif


	if (dim1 != A.dim2)
	{
		std::cerr << "ERREUR multiply : nbRow of this = " << getNbRow(); 
       		std::cerr <<" while nbCol of A = " << A.getNbCol() << std::endl;
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
		
	mat = A.mat * mat;
	resize(A.dim1, dim2);

#ifdef __COMPILE_TIMERS__
 t_mult_left.stop();
#endif

 }
 
 void faust_mat::scalarMultiply(faust_real const lambda)
 {
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.start();
#endif
	 mat = lambda * mat;
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
 }
 
 
 
 void faust_mat::add(faust_mat const& A)
 {  
#ifdef __COMPILE_TIMERS__
t_add.start();
#endif
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		std::cerr << "ERREUR add : matrix dimension not equal" << std::endl; 
        	exit( EXIT_FAILURE);
	}
	mat = mat + A.mat;
#ifdef __COMPILE_TIMERS__
t_add.stop();
#endif
 }

 
 void faust_mat::sub(faust_mat const& A)
 {  
#ifdef __COMPILE_TIMERS__
t_sub.start();
#endif
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		std::cerr << "ERREUR sub : matrix dimension not equal" << std::endl; 
        	exit( EXIT_FAILURE);
	}
	mat = mat - A.mat;
#ifdef __COMPILE_TIMERS__
t_sub.stop();
#endif
 }


 
  // Affichage
  void faust_mat::Display() const
  {     std::cout << "nb_row=" << getNbRow() << endl;
        std::cout << "nb_col=" << getNbCol()   <<endl;  
	std::cout << mat <<endl; 
  }
  
  void faust_mat::print_file(const char* filename)const
  {
	ofstream fichier;
	fichier.open(filename);
	for (int i=0 ; i<getNbRow() ; i++)
	{
		for (int j=0 ; j<getNbCol() ; j++)
			fichier << setprecision(20) <<mat(i,j) << " ";
		fichier << endl;
	}
	fichier.close();
  }

  
    /// SURCHARGE OPERATEUR ///
  // affectation
  void faust_mat::operator=(faust_mat const& A)
  {
	  mat = A.mat;
	  dim1 = A.dim1;
	  dim2 = A.dim2;
  }
  

void faust_mat::init_from_file(const char* filename)
{
#ifdef __COMPILE_TIMERS__
t_print_file.start();
#endif
  // la premiere ligne contient 2 entiers : dim1 et dim2
  // chacune des autres lignes contient une valeur par ligne
  // suivant la premiere dimension puis la deuxieme 

  ifstream* vec_stream;
  vec_stream = new ifstream(filename);
  istream_iterator<faust_real> start(*vec_stream), eos;
  vector<faust_real> vec(start, eos); 

  if((vec[0]*vec[1]+2) != vec.size())
  {
	  cerr << "Error in faust_mat::init_from_file : impossible to read matrix from file " << filename << endl;
	  exit(EXIT_FAILURE);
  }
  resize(vec[0],vec[1]);
  memcpy(getData(), &vec[2], sizeof(faust_real) * dim1 * dim2); 
#ifdef __COMPILE_TIMERS__
t_print_file.stop();
#endif
}


bool operator==(faust_mat const& A, faust_mat const& B)
{
	return A.isEqual(B);
}


#ifdef __COMPILE_TIMERS__
faust_timer faust_mat::t_constr;
faust_timer faust_mat::t_get_coeff;
faust_timer faust_mat::t_get_coeffs;
faust_timer faust_mat::t_set_coeff;
faust_timer faust_mat::t_set_coeffs;
faust_timer faust_mat::t_set_coeffs2;
faust_timer faust_mat::t_resize;
faust_timer faust_mat::t_check_dim;
faust_timer faust_mat::t_max;
faust_timer faust_mat::t_transpose;
faust_timer faust_mat::t_mult_right;
faust_timer faust_mat::t_mult_left;
faust_timer faust_mat::t_scalar_multiply;
faust_timer faust_mat::t_add;
faust_timer faust_mat::t_sub;
faust_timer faust_mat::t_print_file;
faust_timer faust_mat::t_multiply;
faust_timer faust_mat::t_gemm;
faust_timer faust_mat::t_add_ext;

void faust_mat::print_timers()const
{
   cout << "timers in faust_mat :" << endl;
   cout << "t_constr          = " << t_constr.get_time()          << " s for "<< t_constr.get_nb_call()          << " calls" << endl;
   cout << "t_get_coeff       = " << t_get_coeff.get_time()       << " s for "<< t_get_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_get_coeffs      = " << t_get_coeffs.get_time()      << " s for "<< t_get_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeff       = " << t_set_coeff.get_time()       << " s for "<< t_set_coeff.get_nb_call()       << " calls" << endl;
   cout << "t_set_coeffs      = " << t_set_coeffs.get_time()      << " s for "<< t_set_coeffs.get_nb_call()      << " calls" << endl;
   cout << "t_set_coeffs2     = " << t_set_coeffs2.get_time()     << " s for "<< t_set_coeffs2.get_nb_call()     << " calls" << endl;
   cout << "t_resize          = " << t_resize.get_time()          << " s for "<< t_resize.get_nb_call()          << " calls" << endl;
   cout << "t_check_dim       = " << t_check_dim.get_time()       << " s for "<< t_check_dim.get_nb_call()       << " calls" << endl;
   cout << "t_max             = " << t_max.get_time()             << " s for "<< t_max.get_nb_call()             << " calls" << endl;
   cout << "t_transpose       = " << t_transpose.get_time()       << " s for "<< t_transpose.get_nb_call()       << " calls" << endl;
   cout << "t_mult_right      = " << t_mult_right.get_time()      << " s for "<< t_mult_right.get_nb_call()      << " calls" << endl;
   cout << "t_mult_left       = " << t_mult_left.get_time()       << " s for "<< t_mult_left.get_nb_call()       << " calls" << endl;
  cout << "t_scalar_multiply = " << t_scalar_multiply.get_time() << " s for "<< t_scalar_multiply.get_nb_call() << " calls" << endl;
   cout << "t_add             = " << t_add.get_time()             << " s for "<< t_add.get_nb_call()             << " calls" << endl;
   cout << "t_sub             = " << t_sub.get_time()             << " s for "<< t_sub.get_nb_call()             << " calls" << endl;
cout << "t_print_file      = " << t_print_file.get_time()      << " s for "<< t_print_file.get_nb_call()      << " calls" << endl<<endl;
                                                      
   cout << "timers in faust_mat / LinearAlgebra :" << endl;
   cout << "t_multiply        = " << t_multiply.get_time()        << " s for "<< t_multiply.get_nb_call()        << " calls" << endl;
   cout << "t_gemm            = " << t_gemm.get_time()            << " s for "<< t_gemm.get_nb_call()            << " calls" << endl;
   cout << "t_add_ext         = " << t_add_ext.get_time()         << " s for "<< t_add_ext.get_nb_call()         << " calls" << endl<<endl<<endl;
}
#endif
