#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 

#include <Eigen/Dense>

#include "faust_constant.h"
#include <vector>
#include <iterator>
#include "faust_mat_generic.h"
#include "faust_exception.h"
#include <iostream>

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

/*! \class faust_mat
   * \brief template class representing dense matrix
   *
   *  This class implements basic linear algebra operation (addition, multiplication, frobenius and spectral norm...)
   *
   * The matrix format is ColMajor.
   *
   *\tparam T scalar numeric type, e.g float or double
   */



template<typename T> class faust_mat;
template<typename T> class faust_spmat;
template<typename T> class faust_vec;
template<typename T> class faust_core;

template<typename T>
void multiply(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
	

template<typename T>	
void add(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);

template<typename T>
void gemm(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T& alpha, const T& beta, char  typeA, char  typeB);

template<typename T>
void multiply(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);
  
template<typename T>  
 void gemv(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);









template<typename T>
class faust_mat : public faust_mat_generic
{
	// declarer comme friend toutes les classes derivant de la template class faust_mat
	template<class> friend class faust_mat;
public:
	static const char * name;  

      /*!
     *  \brief Constructor
     *
     *  faust_mat constructor
     *
     *  \tparam data : pointer to the data array of the matrix
		\tparam	   nbRow : number of row of the matrix
		\tparam	   nbCol : number of column of the matrix
     */	
  faust_mat(const T  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol );
  faust_mat() : faust_mat_generic(), mat(0,0), isIdentity(false), isZeros(false) {}
  /*!
     *  \brief Constructor
     *
     *  faust_mat copy constructor
     *
     *  \tparam A : another faust_mat
     */	
  faust_mat(const faust_mat<T> & A) : faust_mat_generic(A.dim1,A.dim2), mat(A.mat), isIdentity(A.isIdentity), isZeros(A.isZeros) {}
	template<typename U>
   faust_mat(const faust_mat<U> & A){this->operator=(A);}
   template<typename U>
   faust_mat(const faust_spmat<U> & A){this->operator=(A);}
  faust_mat(const faust_spmat<T> & A){this->operator=(A);}

  faust_mat(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol) : faust_mat_generic(nbRow,nbCol), mat(nbRow,nbCol), isIdentity(false), isZeros(false){}
  faust_mat(const faust_unsigned_int nbRow) : faust_mat_generic(nbRow,nbRow), mat(nbRow,nbRow), isIdentity(false), isZeros(false){}

 ~faust_mat(){resize(0,0);}
	
  
  /*!
     *  \brief 
	 * resize the faust_mat
		\tparam	   nbRow : new number of row of the matrix
		\tparam	   nbCol : new number of column of the matrix
		\warning nbRow and nbCol must be greater or equal to 0
     */		
  void resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol);
  
  /*!
     *  \brief 
	 * resize the faust_mat
		\tparam	   nbRow : new number of row and column of the matrix
	 \warning nbRow must be greater or equal to 0 	
     */
  void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
  
  /*!
     *  \brief 
	 * check if the dimension of the matrix are consistent,
	 * if not throws an error
     */
  void check_dim_validity();
  
  /*!
     *  \brief 
	 * set the matrix to the zero matrix
     */
  void setZeros();
  
    /*!
     *  \brief 
	 * set the matrix to the one diagonal matrix
     */
  void setEyes();
  
  /*!
     *  \brief 
	 * access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing
	 *\tparam i : position
	 *\return : ith coefficient of the matrix
     */	
  T& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}
  
  /*!
     *  \brief 
	 * access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing
	 *\tparam i : position
	 *\return : read-only ith coefficient of the matrix
     */	
  const T& operator[](faust_unsigned_int i)const{return mat.data()[i];}

   /*!
     *  \brief 
	 * access to the ith coefficient of the matrix pointer, Colmajor format and zero indexing
	 *\tparam i : position
	 *\return : read-only ith coefficient of the matrix
     */	
  const T& operator()(faust_unsigned_int i)const{return mat.data()[i];}
  
  /*!
     *  \brief 
	 * access to (i,j) coefficient of the matrix pointer in zero indexing
	 *\tparam i : row position
	 *\tparam j : col position
	 *\return : read-only (i,j) coefficient of the matrix
     */	
  const T& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*dim1+i];}


   void operator*=(const faust_spmat<T>& M);
   void operator+=(const faust_spmat<T>& M);
   void operator-=(const faust_spmat<T>& M);

   void multiplyLeft(const faust_spmat<T>& M);


  T* getData(){isZeros=false; isIdentity=false;return mat.data();} 
  const T* getData()const{return mat.data();} 

 

  

  bool isEqual(const faust_mat<T> & B) const;
  bool isEqual(const faust_mat<T> & B, T threshold) const;



  /*!
     *  \brief 
	 * initialize faust_mat from text file 
	 *\tparam filename : name of the file
	 * 
	 * the first line of the file contains 2 integer : the number of row and the number of column <br>
	 * all the other line contains one coefficient
	 * in ColMajor access
     */	  
void init_from_file(const char* filename);


  
  

  

  void abs() {mat=mat.cwiseAbs();}
  

  
  
 /*!
     *  \brief 
	 * compute the Frobenius norm of the faust_mat
	 *\return  the Frobenius norm
     */	
  T norm() const {return mat.norm();}
  
  /*!
     *  \brief 
	 * normalize the matrix according to its Frobenius norm
     */
  void normalize() {scalarMultiply(1.0/norm());}
   
   /*!
     *  \brief 
     *
     *  return the spectral norm (maximum singular value in absolute value) using Eigen library, very slow	
	*
	*		   
		\return	the estimated spectral norm   
     */
  T spectralNorm() const;

     /*!
     *  \brief 
     *  return the estimated spectral norm (maximum singular value in absolute value) using power iteration algorithm	
     *
     *  \param nbr_iter_max : maximum number of iteration for the power algo
		\param threshold : threshold until convergence
		\param flag : convergence flag
	*		   
		\return	the estimated spectral norm   
    *
	* See also, template<typename T> 
T power_iteration(const faust_mat<T> & A, const faust_unsigned_int nbr_iter_max,T threshold,faust_int & flag);
	*/	  
  T spectralNorm(const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag) const;
  
  /*!
     *  \brief 
	 * compute the trace of the faust_mat
	 *\return  the trace
     */
  T trace() const {return mat.trace();}
  
  /*!
     *  \brief 
	 * transpose the faust_mat
     */
  void transpose();
  
  /*!
     *  \brief 
	 * replace this by (this) * A  
     */
  void multiplyRight(faust_mat<T> const& A);
  

  /*!
     *  \brief 
	 * replace this by A * (*this)  
     */
  void multiplyLeft(faust_mat<T> const& A);
  
  /*!
     *  \brief 
	 * replace this by lambda * (*this)  
     */	
  void scalarMultiply(T const lambda);
  
  /*!
     *  \brief 
	 * replace this by lambda * (*this) using element by element multiplication  
     */
  void scalarMultiply(faust_mat<T> const& A);
  
  /*!
     *  \brief 
	 * (*this) = (*this) + A   
     */
  void add(faust_mat<T> const& A);

  /*!
     *  \brief 
	 * (*this) = (*this) - A   
     */	
  void sub(faust_mat<T> const& A);

  

  /*!
     *  \brief 
	 * displays the faust_mat  
     */
  void Display() const;
  
  /*!
     *  \brief 
	 * write faust_mat into text file 
	 *\tparam filename : name of the file
	 * 
	 * the first line of the file contains 2 integer : the number of row and the number of column
	 * all the other line contains one coefficient
	 * in ColMajor access of the faust_mat
     */	
  void print_file(const char* filename)const;

  


  void operator=(faust_mat<T> const& A);
  
  template<typename U>
  void operator=(faust_mat<U> const& A);
  template<typename U>
  void operator=(faust_spmat<U> const& A){faust_spmat<T> AT(A);this->operator=(AT);};
  
  void operator=(faust_spmat<T> const& A);


  void operator-=(faust_mat<T> const& A){sub(A);}
  void operator+=(faust_mat<T> const& A){add(A);}

  void operator*=(faust_mat<T> const& A){multiplyRight(A);}
  

  void operator*=(T lambda){scalarMultiply(lambda);}
  void operator/=(T lambda){scalarMultiply(1.0/lambda);}


  

  friend void multiply<>(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
	

	
 
  friend void gemm<>(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T& alpha, const T& beta, char  typeA, char  typeB);
  friend void multiply<>(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);
  friend void gemv<>(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);
  

  bool estIdentite()const{return isIdentity;}
  bool estNulle()const{return isZeros;}
  
  private: 
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;
       bool isIdentity;
       bool isZeros;
	   static const char * class_name;
	   
	   	   

#ifdef __COMPILE_TIMERS__
  public:
	  faust_timer t_local_muliplyLeft;
  //temporary members
      static faust_timer t_constr;
      static faust_timer t_get_coeff;
      static faust_timer t_get_coeffs;
      static faust_timer t_set_coeff;
      static faust_timer t_set_coeffs;
      static faust_timer t_set_coeffs2;
      static faust_timer t_resize;
      static faust_timer t_check_dim;
      static faust_timer t_max;
      static faust_timer t_transpose;
      static faust_timer t_mult_right;
      static faust_timer t_mult_left;
      static faust_timer t_scalar_multiply;
      static faust_timer t_add;
      static faust_timer t_sub;
      static faust_timer t_print_file;
	  static faust_timer t_spectral_norm;
	  static faust_timer t_spectral_norm2;

      static faust_timer t_power_iteration;
	  static faust_timer t_multiply;
      static faust_timer t_gemm;
      static faust_timer t_add_ext;

  void print_timers()const;
#endif
};


#include "faust_mat.hpp"




#endif
