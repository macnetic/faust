#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 

#include <Eigen/Dense>

#include "faust_constant.h"
#include <vector>
#include <iterator>
#include "faust_mat_generic.h"
#include "faust_exception.h"

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

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
 faust_vec<T> solve(const faust_mat<T> & A, const faust_vec<T> & v);







template<typename T>
class faust_mat : public faust_mat_generic
{
public:
	static const char * name;  

  /// Constructeurs ///
  faust_mat(const Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
  faust_mat(const T  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol );	
  faust_mat() : faust_mat_generic(), mat(0,0), isIdentity(false), isZeros(false) {}
  faust_mat(const faust_mat<T> & A) : faust_mat_generic(A.dim1,A.dim2), mat(A.mat), isIdentity(A.isIdentity), isZeros(A.isZeros) {}
  faust_mat(const faust_spmat<T> & A){this->operator=(A);}

  faust_mat(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol) : faust_mat_generic(nbRow,nbCol), mat(nbRow,nbCol), isIdentity(false), isZeros(false){}
  faust_mat(const faust_unsigned_int nbRow) : faust_mat_generic(nbRow,nbRow), mat(nbRow,nbRow), isIdentity(false), isZeros(false){}

 ~faust_mat(){resize(0,0);}
	
	
  /// GETTEUR SETTEUR ///
  // faust_unsigned_int getNbRow() const {return dim1;}
  // faust_unsigned_int getNbCol() const {return dim2;}
  /*T getCoeff(const faust_unsigned_int i,const faust_unsigned_int j) const;
  void getCoeffs(std::vector<T> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const;
  void setCoeff(const T & value,const int id_row, const int id_col);
  void setCoeffs(const T value,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  void setCoeffs(const std::vector<T> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col);*/
  

  void resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol);
  void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
  
  void check_dim_validity();
  
  // (*this) = la matrice nulle
  //void setZeros() {mat.setZero();isZeros=true;}
  void setZeros();
  
  // (*this) = identite, pas forcement carree
  //void setEyes() {mat.setIdentity();if(dim1==dim2)isIdentity=true;}
  void setEyes();

  T& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}

  const T& operator[](faust_unsigned_int i)const{return mat.data()[i];}

  const T& operator()(faust_unsigned_int i)const{return mat.data()[i];}
  const T& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*dim1+i];}


   void operator*=(const faust_spmat<T>& M);
   void operator+=(const faust_spmat<T>& M);
   void operator-=(const faust_spmat<T>& M);

   void multiplyLeft(const faust_spmat<T>& M);


  T* getData(){isZeros=false; isIdentity=false;return mat.data();} 
  const T* getData()const{return mat.data();} 

 

  
  /// EGALITE ///
  //bool isZeros() const {return mat.isZero(FAUST_PRECISION);}
  bool isEqual(const faust_mat<T> & B) const;
  bool isEqual(const faust_mat<T> & B, T threshold) const;
  //bool isEyes() const {return mat.isIdentity(FAUST_PRECISION);}
  
void init_from_file(const char* filename);


  
  
  /// OPERATION BASIQUE ///
  
  //arithmetique
  
  T max() const {return mat.maxCoeff();}
  T min() const {return mat.minCoeff();}
  void abs() {mat=mat.cwiseAbs();}
  
  // return the maximum of all coefficients of this and puts in row_id and col_id its location
  T max(std::vector<faust_unsigned_int> & id_row,std::vector<faust_unsigned_int> & id_col) const;
  T min(std::vector<faust_unsigned_int> & id_row,std::vector<faust_unsigned_int> & id_col) const;
  
  
  // frobenius norm
  T norm() const {return mat.norm();}
  void normalize() {scalarMultiply(1.0/norm());}
  // spectral norm, "norm2", equal to the largest singular value  
  T spectralNorm() const;
  T spectralNorm(const faust_unsigned_int nbr_iter_max,T threshold, faust_int & flag) const;
  
  // trace
  T trace() const {return mat.trace();}
  
  //transposition
  void transpose();
  
  // multiply (*this) = (*this) * A
  void multiplyRight(faust_mat<T> const& A);
  // multiply (*this) =  A * (*this)
  void multiplyLeft(faust_mat<T> const& A);
  
  // scalarMultiply (*this) = (*this) * lambda
  void scalarMultiply(T const lambda);
  // (*this)(i,j)=((*this)(i,j)) * A(i,j)	
  void scalarMultiply(faust_mat<T> const& A);
  // (*this) = (*this) + A
  void add(faust_mat<T> const& A);

  // (*this) = (*this) - A
  void sub(faust_mat<T> const& A);

  
  // Affichage
  void Display() const;
  void print_file(const char* filename)const;

  
  /// SURCHARGE OPERATEUR ///
  // affectation
  void operator=(faust_mat<T> const& A);
  void operator=(faust_spmat<T> const& A);


  void operator-=(faust_mat<T> const& A){sub(A);}
  void operator+=(faust_mat<T> const& A){add(A);}

  void operator*=(faust_mat<T> const& A){multiplyRight(A);}
  

  void operator*=(T lambda){scalarMultiply(lambda);}
  void operator/=(T lambda){scalarMultiply(1.0/lambda);}


  
  ////////////////// friends //////////////////////
  // intra classe//
  // friend void add<>(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
  friend void multiply<>(const faust_mat<T> & A, const faust_mat<T> & B, faust_mat<T> & C);
	

	
 
  friend void gemm<>(const faust_mat<T> & A,const faust_mat<T> & B, faust_mat<T> & C,const T& alpha, const T& beta, char  typeA, char  typeB);
  friend void multiply<>(const faust_core<T> & A, const faust_mat<T> & B, faust_mat<T> & C,const T & alpha, char typeA, char typeMult);
  friend void gemv<>(const faust_mat<T> & A,const faust_vec<T> & x,faust_vec<T> & y,const T & alpha, const T & beta, char typeA);
  friend faust_vec<T> solve<>(const faust_mat<T> & A, const faust_vec<T> & v);
  ///////////friend faust_spmat<T>::operator=(faust_mat<T> const& S);
  bool estIdentite(){return isIdentity;}
  bool estNulle(){return isZeros;}
  
  private: 
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> mat;
  //Eigen::Matrix<T,0,0> mat;
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

 //bool operator==(faust_mat const& A, faust_mat const& B);


#endif
