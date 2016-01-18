
//#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Sparse>
#include "faust_exception.h"

#ifdef __GEMM_WITH_OPENBLAS__
	#include "cblas_algebra.h"
#endif




#ifdef __GEMM_WITH_MKL__
	#include "mkl_spblas.h"
#endif

template<typename U> class faust_mat;
using namespace std;
template<typename T>
const char * faust_mat<T>::class_name = "faust_mat<T>::";


/// CONSTRUCTEUR ///
template<typename T>
faust_mat<T>::faust_mat(const Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> & mat_) : 
     // dim1(mat_.rows()), dim2(mat_.cols()),mat(mat_),isIdentity(false),isZeros(false){}
     faust_mat_generic(mat_.rows(),mat_.cols()),mat(mat_),isIdentity(false),isZeros(false){}

	 

	 
	 
	 
	 
template<typename T>	 
faust_mat<T>::faust_mat(const T  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol ) : faust_mat_generic(nbRow,nbCol), mat(nbRow,nbCol), isIdentity(false),isZeros(false)
  {

#ifdef __COMPILE_TIMERS__
t_constr.start();
#endif
          /*
	  int i,j;
	  if ((nbRow < 0) || (nbCol < 0))
	  {
		cerr << "ERREUR une dimension ne peut etre negative" << endl;; 
        	exit( EXIT_FAILURE);  
	  }
	  
	dim1 = nbRow;
	dim2 = nbCol;
	Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> m(nbRow,nbCol);
		
	for (int j=0;j<nbCol;j++)
		for (int i=0;i<nbRow;i++)
			m(i,j) = mat_[i+(j*nbRow)];
	mat = m;*/

	memcpy(getData(), data_, nbRow*nbCol*sizeof(T));	

	  
#ifdef __COMPILE_TIMERS__
t_constr.stop();
#endif

  }

/// GETTEUR SETTEUR ///


	

template<typename T> 
void faust_mat<T>::resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol)
{	
#ifdef __COMPILE_TIMERS__
t_resize.start();
#endif
		if ( (nbRow <0) || (nbCol <0) )
		{
			//cerr << "ERREUR resize : les nouvelles dimensions doivent etre strictement positive" << endl;
			//exit( EXIT_FAILURE);
			handleError(class_name, "resize : new dimensions must be positive");
		}
		else if ((dim1 != nbRow) || (dim2 != nbCol))
		{
			faust_mat_generic::resize(nbRow,nbCol);
			mat.resize(nbRow,nbCol);
		}

        isZeros = false;
        isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_resize.stop();
#endif
}


template<typename T>
void faust_mat<T>::check_dim_validity()
 {

#ifdef __COMPILE_TIMERS__
t_check_dim.start();
#endif

	bool verifSize = (getNbCol() == mat.cols()) &&  (getNbRow() == mat.rows());
	
	if (!verifSize)
	{
		handleError(class_name, "check_dim_validity : Size incompatibility in the faust_mat");
	}
#ifdef __COMPILE_TIMERS__
t_check_dim.stop();
#endif

 }

template<typename T>
void faust_mat<T>::setZeros()
{
	memset(getData(), 0, sizeof(T) * dim1*dim2);
	isZeros = true;
	isIdentity = false;
}

template<typename T>
void faust_mat<T>::setEyes()
{
	setZeros();
	T* ptr_data = getData();
	for (int i=0 ; i<std::min(dim1,dim2); i++)
		ptr_data[i*dim1+i] = 1.0;
	if (dim1 == dim2)
		isIdentity = true;
	isZeros = false;
	
}

 /// EGALITE ///
template<typename T>
bool faust_mat<T>::isEqual(const faust_mat<T> & B) const 
 {
	if ((getNbCol() != B.getNbCol()) || (getNbRow() != B.getNbRow()))	
		handleError(class_name, "isEqual : dimension of the 2 matrix are not the same\n");	
	

	if (isZeros)
		return B.isZeros;
	else if (B.isZeros)
		return isZeros;
	else
		return mat.isApprox(B.mat,FAUST_PRECISION);
 }

 
template<typename T> 
bool faust_mat<T>::isEqual(const faust_mat<T> & B, T threshold) const
{
	if ((getNbCol() != B.getNbCol()) || (getNbRow() != B.getNbRow()))
	{
		handleError(class_name, "isEqual : dimension of the 2 matrix are not the same\n");	
	}
	bool egalite =true;
	for (int i=0;i<getNbRow();i++)
	{
		for (int j=0;j<getNbCol();j++)
		{
			if (std::abs(mat(i,j)==0))
			{
			
				if (	(std::abs(mat(i,j)-B.mat(i,j))) > threshold )
				{
					egalite = false;
					std::cout<<" i : "<<i<<" j : "<<j<<std::endl;
				}
			}else
			{
				if (	(std::abs(mat(i,j)-B.mat(i,j))/std::abs(mat(i,j))) > threshold )
				{
					egalite = false;
					std::cout<<" i : "<<i<<" j : "<<j<<std::endl;
				}
			}
		}
	}
	
	return egalite;
	
}
 
 /// OPERATION BASIQUE ///
 




 
 
template<typename T> 
void faust_mat<T>::transpose()
 {

#ifdef __COMPILE_TIMERS__
t_transpose.start();
#endif

	if(isZeros || isIdentity)
	{
		resize(dim2,dim1);
		#ifdef __COMPILE_TIMERS__
			t_transpose.stop();
		#endif
		return;
	}

	mat = mat.transpose().eval();
	faust_unsigned_int dim1_copy = dim1;
	dim1 = dim2;
        dim2 = dim1_copy; 
        
	

#ifdef __COMPILE_TIMERS__
t_transpose.stop();
#endif

 }

template<typename T> 
 void faust_mat<T>::multiplyRight(faust_mat<T> const& A)
 { 

#ifdef __COMPILE_TIMERS__
t_mult_right.start();
#endif


	if (dim2 != A.dim1)
	{
		handleError(class_name, "multiplyRight : dimension conflict between matrix");		
	}

	if(A.isIdentity)
	{	
		#ifdef __COMPILE_TIMERS__
			t_mult_right.stop();
		#endif
		return;
	}

	if(isZeros || A.isZeros)
	{	
		//std::cout<<"zero"<<std::endl;
		resize(dim1,A.dim2);
		T *const ptr_data_dst = getData();
		memset(ptr_data_dst, 0, sizeof(T) * dim1*dim2);
		isZeros = true;
		isIdentity = false;
		#ifdef __COMPILE_TIMERS__
			t_mult_right.stop();
		#endif
		return;
	}

	if(isIdentity)
	{	
		this->operator=(A);
		#ifdef __COMPILE_TIMERS__
			t_mult_right.stop();
		#endif
		return;
	}
	
	faust_mat this_copy((*this));	
	gemm<T>(this_copy,A,(*this),1.0,0.0,'N','N');


#ifdef __COMPILE_TIMERS__
t_mult_right.stop();
#endif

 }


template<typename T> 
 void faust_mat<T>::multiplyLeft(faust_mat<T> const& A)
 { 

#ifdef __COMPILE_TIMERS__
 t_mult_left.start();
#endif


	if (dim1 != A.dim2)
	{
		handleError(class_name, "multiplyLeft : dimension conflict between matrix");		
	}

	if(A.isIdentity)
	{
		#ifdef __COMPILE_TIMERS__
			t_mult_left.stop();
		#endif
		return;
	}

	if(isZeros || A.isZeros)
	{
		resize(A.dim1,dim2);
		T *const ptr_data_dst = getData();
		memset(ptr_data_dst, 0, sizeof(T) * dim1*dim2);
		isZeros = true;
		isIdentity = false;
		#ifdef __COMPILE_TIMERS__
			t_mult_left.stop();
		#endif
		return;
	}
	
	if(isIdentity)
	{
		this->operator=(A);
		#ifdef __COMPILE_TIMERS__
			t_mult_left.stop();
		#endif
		return;
	}

		
	faust_mat this_copy((*this));	
	gemm<T>(A,this_copy,(*this),1.0,0.0,'N','N');

#ifdef __COMPILE_TIMERS__
 t_mult_left.stop();
#endif

 }
 
 
 template<typename T>
 T faust_mat<T>::spectralNorm() const 
 {	
	#ifdef __COMPILE_TIMERS__
			t_spectral_norm.start();
	#endif
	
	 T res=mat.operatorNorm();
	 
	#ifdef __COMPILE_TIMERS__
		t_spectral_norm.stop();
	#endif
	return res;
 }
 
 template<typename T>
 T faust_mat<T>::spectralNorm(const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag) const
{	
	#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.start();
	#endif
	if(isZeros)
	{
		flag = -2;
		#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.stop();
		#endif
		return 0;
	}
		
	if(isIdentity)
	{
		flag = -3;
		#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.stop();
		#endif
		return 1;
	}
		
	faust_unsigned_int nb_row = getNbRow();
	faust_unsigned_int nb_col = getNbCol();
		
	
		faust_mat<T> AtA;
		if (nb_row <= nb_col)
		{	
			gemm<T>((*this),(*this),AtA,1.,0,'N','T');
		}else
		{
			gemm<T>((*this),(*this),AtA,1.,0,'T','N');
		}



			T  res=std::sqrt(power_iteration(AtA,nbr_iter_max,threshold,flag));
			

		
		#ifdef __COMPILE_TIMERS__
		t_spectral_norm2.stop();
		#endif
		return res;
			
	}
	
		
		
	
 
 
 
 
 
 
 
 template<typename T>
 void faust_mat<T>::scalarMultiply(T const lambda)
 {
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.start();
#endif
	 mat = lambda * mat;
#ifdef __COMPILE_TIMERS__
t_scalar_multiply.stop();
#endif
 }
 
 
 template<typename T>
 void faust_mat<T>::add(faust_mat<T> const& A)
 {  
#ifdef __COMPILE_TIMERS__
t_add.start();
#endif
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		handleError(class_name, "add : matrix dimension not equal");	
	}
	mat = mat + A.mat;
        isZeros = false;
        isIdentity = false;
#ifdef __COMPILE_TIMERS__
t_add.stop();
#endif
 }

 template<typename T>
 void faust_mat<T>::sub(faust_mat<T> const& A)
 {  
#ifdef __COMPILE_TIMERS__
t_sub.start();
#endif
	if ((getNbCol() != A.getNbCol()) || (getNbRow() != A.getNbRow()))
	{
		handleError(class_name, "sub : matrix dimension not equal");	
	}
	mat = mat - A.mat;

        isZeros = false;
        isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_sub.stop();
#endif
 }












 
  // Affichage
  template<typename T>
  void faust_mat<T>::Display() const
  {     //std::cout << "nb_row=" << getNbRow() << endl;
        //std::cout << "nb_col=" << getNbCol()   <<endl;  
	std::cout<<"DIM1 : "<<dim1<<" DIM2 : "<<dim2<<std::endl;
	std::cout<<"iszeros : "<<isZeros<<" isidentity : "<<isIdentity<<std::endl;
	std::cout << mat <<endl; 
	
  }
  
  template<typename T>
  void faust_mat<T>::print_file(const char* filename)const
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
  template<typename T>
  void faust_mat<T>::operator=(faust_mat<T> const& A)
  {
	  mat = A.mat;
	  dim1 = A.dim1;
	  dim2 = A.dim2;
          isZeros = A.isZeros;
          isIdentity = A.isIdentity;
  }
  
  template<typename T>
  template<typename U>
  void faust_mat<T>::operator=(faust_mat<U> const& A)
  {	
	resize(A.dim1,A.dim2);
	// mat = A.mat.cast<T>();
	for (int i=0;i<dim1*dim2;i++)
			(*this)[i]=(T) A(i);
    isZeros = A.isZeros;
    isIdentity = A.isIdentity;	
  }
  
  template<typename T>
  void faust_mat<T>::operator=(faust_spmat<T> const& S)
  {
          resize(S.getNbRow(),S.getNbCol());
	  setZeros();
          T*const ptr_data = getData();
          for(int i=0 ; i< S.mat.outerSize() ; i++)
             for(typename Eigen::SparseMatrix<T,Eigen::RowMajor>::InnerIterator it(S.mat,i); it; ++it)
                ptr_data[it.col() * dim1 + it.row()] = it.value();
          isZeros = false;
          isIdentity = false;
  }
 

 template<typename T>
void faust_mat<T>::operator*=(const faust_spmat<T>& S)
{
	if(dim2 != S.dim1)
	{
		handleError(class_name, "operator*= : incorrect matrix dimensions");
	}

	if (isIdentity)
	{
		this->operator=(S);
		isIdentity = false;
		isZeros = false;
	}
	else if (isZeros)
	{
		resize(dim1, S.dim2);
		setZeros();
	}
	else
	{
		mat = mat * S.mat;
		dim2 = S.dim2;
	}
	
}


template<typename T>
void faust_mat<T>::operator+=(const faust_spmat<T>& S)
{
	if(dim1!=S.dim1 || dim2!=S.dim2)
	{
		handleError(class_name,"operator+= : incorrect matrix dimensions");
	}
	mat += S.mat;
	isIdentity = false;
	isZeros = false;
}

template<typename T>
void faust_mat<T>::operator-=(const faust_spmat<T>& S)
{
	if(dim1!=S.dim1 || dim2!=S.dim2)
	{
		handleError(class_name,"operator-= : incorrect matrix dimensions");
	}
	mat -= S.mat;
	isIdentity = false;
	isZeros = false;
}


template<typename T>
void faust_mat<T>::scalarMultiply(faust_mat<T> const& A)
{
	if(dim1!=A.dim1 || dim2!=A.dim2)
	{
		handleError(class_name,"scalarMultiply : incorrect matrix dimensions\n");
	}
	mat = (mat.array() * A.mat.array()).matrix();
	isIdentity = false;
	isZeros = false;
}


template<typename T>
void faust_mat<T>::multiplyLeft(const faust_spmat<T>& S)
{
	if(S.dim2 != dim1)
	{
		//std::cerr << "Error in faust_mat<T>::operator*= : incorrect matrix dimensions" << std::endl;
		//exit(EXIT_FAILURE);
		handleError(class_name,"multiplyLeft : incorrect matrix dimensions\n");
	}

	if (isIdentity)
	{
		this->operator=(S);
		isIdentity = false;
		isZeros = false;
	}
	else if (isZeros)
	{
		resize(S.dim1, dim2);
		setZeros();
	}
	else
	{	
		#ifdef __GEMM_WITH_MKL__
			int dim1S = S.getNbRow();
			int dim2S = S.getNbCol();
			int dim1this = dim1;
			int dim2this = dim2;
			int dim1C = dim1S;
			int dim2C = dim2this;
			T alpha = 1.0,beta = 0.0;
			char matdescrS[8];
			matdescrS[0]='G';
			for (int i=1 ; i<=6 ;i++)
			matdescrS[i] = 'C';
			matdescrS[7] = '\0';
			
			faust_spmat<T> Scopy(S);
			faust_mat<T> C(dim1C,dim2C);
			if (!Scopy.isCompressedMode())
			{
				Scopy.makeCompression();
				std::cout<<class_name<<"multiplyLeft : sparse matrix S is not compressed (possible lack of time)"<<std::endl;
			}
			#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.start();
			#endif
			#ifdef FAUST_SINGLE	
				mkl_scscmm("N",&dim1S,&dim2C,&dim2S,&alpha,matdescrS,Scopy.getValuePtr(),Scopy.getInnerIndexPtr(),Scopy.getOuterIndexPtr(),&(Scopy.getOuterIndexPtr()[1]),getData(),&dim1this,&beta,C.getData(),&dim1C);
			#else
				mkl_dcscmm("N",&dim1S,&dim2C,&dim2S,&alpha,matdescrS,Scopy.getValuePtr(),Scopy.getInnerIndexPtr(),Scopy.getOuterIndexPtr(),&(Scopy.getOuterIndexPtr()[1]),getData(),&dim1this,&beta,C.getData(),&dim1C);
			#endif	
			(*this) = C;
			
			#ifdef __COMPILE_TIMERS__
			t_local_multiplyLeft.stop();
			cout << "1 "<<setprecision(10)<<t_local_multiplyLeft.get_time()<<endl;
			t_local_multiplyLeft.reset();
			#endif	
		#else
			
			mat = S.mat * mat;
			dim1 = S.dim1;
		#endif
	}
	
}
 
template<typename T>
void faust_mat<T>::init_from_file(const char* filename)
{
#ifdef __COMPILE_TIMERS__
t_print_file.start();
#endif
  // la premiere ligne contient 2 entiers : dim1 et dim2
  // chacune des autres lignes contient une valeur par ligne
  // suivant la premiere dimension puis la deuxieme 

  ifstream* vec_stream;
  vec_stream = new ifstream(filename);
  if (!vec_stream->is_open())
	  handleError(class_name, "init_from_file : unable to open file");
  istream_iterator<T> start(*vec_stream), eos;
  vector<T> vec(start, eos); 

  if((vec[0]*vec[1]+2) != vec.size())
  {
	  handleError(class_name, "init_from_file : problem with the file");
  }
  resize(vec[0],vec[1]);
  memcpy(getData(), &vec[2], sizeof(T) * dim1 * dim2); 
  
  isZeros = false;
  isIdentity = false;

#ifdef __COMPILE_TIMERS__
t_print_file.stop();
#endif
}



#ifdef __COMPILE_TIMERS__
template<typename T>
faust_timer faust_mat<T>::t_constr;
template<typename T>
faust_timer faust_mat<T>::t_get_coeff;
template<typename T>
faust_timer faust_mat<T>::t_get_coeffs;
template<typename T>
faust_timer faust_mat<T>::t_set_coeff;
template<typename T>
faust_timer faust_mat<T>::t_set_coeffs;
template<typename T>
faust_timer faust_mat<T>::t_set_coeffs2;
template<typename T>
faust_timer faust_mat<T>::t_resize;
template<typename T>
faust_timer faust_mat<T>::t_check_dim;
template<typename T>
faust_timer faust_mat<T>::t_max;
template<typename T>
faust_timer faust_mat<T>::t_transpose;
template<typename T>
faust_timer faust_mat<T>::t_mult_right;
template<typename T>
faust_timer faust_mat<T>::t_mult_left;
template<typename T>
faust_timer faust_mat<T>::t_scalar_multiply;
template<typename T>
faust_timer faust_mat<T>::t_add;
template<typename T>
faust_timer faust_mat<T>::t_sub;
template<typename T>
faust_timer faust_mat<T>::t_print_file;
template<typename T>
faust_timer faust_mat<T>::t_multiply;
template<typename T>
faust_timer faust_mat<T>::t_gemm;
template<typename T>
faust_timer faust_mat<T>::t_add_ext;

template<typename T>
faust_timer faust_mat<T>::t_spectral_norm;
template<typename T>
faust_timer faust_mat<T>::t_spectral_norm2;
template<typename T>
faust_timer faust_mat<T>::t_power_iteration;

template<typename T>
void faust_mat<T>::print_timers()const
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
