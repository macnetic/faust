#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 
//#include "faust_spmat.h"
#include <Eigen/Dense>

#include "faust_constant.h"
#include <vector>
#include <iterator>

#ifdef __COMPILE_TIMERS__
  #include "faust_timer.h"
#endif

class faust_vec;//forward declaration of faust_vec class
class faust_spmat;
class faust_core;

class faust_mat
{
public:
	static const char * name;  

  /// Constructeurs ///
  faust_mat(const Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
  faust_mat(const faust_real  *data_,const faust_unsigned_int nbRow, const faust_unsigned_int nbCol );	
  faust_mat() : dim1(0), dim2(0), mat(0,0), isIdentity(false), isZeros(false) {}
  faust_mat(const faust_mat & A) : dim1(A.dim1), dim2(A.dim2), mat(A.mat), isIdentity(A.isIdentity), isZeros(A.isZeros) {}
  faust_mat(const faust_spmat & A){this->operator=(A);}

  faust_mat(const faust_unsigned_int nbRow, const faust_unsigned_int nbCol) : dim1(nbRow), dim2(nbCol), mat(nbRow,nbCol), isIdentity(false), isZeros(false){}
  faust_mat(const faust_unsigned_int nbRow) : dim1(nbRow), dim2(nbRow), mat(nbRow,nbRow), isIdentity(false), isZeros(false){}

 ~faust_mat(){resize(0,0);}
	
	
  /// GETTEUR SETTEUR ///
  faust_unsigned_int getNbRow() const {return dim1;}
  faust_unsigned_int getNbCol() const {return dim2;}
  /*faust_real getCoeff(const faust_unsigned_int i,const faust_unsigned_int j) const;
  void getCoeffs(std::vector<faust_real> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const;
  void setCoeff(const faust_real & value,const int id_row, const int id_col);
  void setCoeffs(const faust_real value,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  void setCoeffs(const std::vector<faust_real> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col);*/
  

  void resize(const faust_unsigned_int nbRow,const faust_unsigned_int nbCol);
  void resize(const faust_unsigned_int nbRow){resize(nbRow,nbRow);}
  
  void check_dim_validity();
  
  // (*this) = la matrice nulle
  //void setZeros() {mat.setZero();isZeros=true;}
  void setZeros();
  
  // (*this) = identite, pas forcement carree
  //void setEyes() {mat.setIdentity();if(dim1==dim2)isIdentity=true;}
  void setEyes();

  faust_real& operator[](faust_unsigned_int i){isZeros=false; isIdentity=false;return mat.data()[i];}

  const faust_real& operator[](faust_unsigned_int i)const{return mat.data()[i];}

  const faust_real& operator()(faust_unsigned_int i)const{return mat.data()[i];}
  const faust_real& operator()(faust_unsigned_int i, faust_unsigned_int j)const{return mat.data()[j*dim1+i];}


   void operator*=(const faust_spmat& M);
   void operator+=(const faust_spmat& M);
   void operator-=(const faust_spmat& M);

   void multiplyLeft(const faust_spmat& M);


  faust_real* getData(){isZeros=false; isIdentity=false;return mat.data();} 
  const faust_real* getData()const{return mat.data();} 

 

  
  /// EGALITE ///
  //bool isZeros() const {return mat.isZero(FAUST_PRECISION);}
  bool isEqual(const faust_mat & B) const;
  bool isEqual(const faust_mat & B, faust_real threshold) const;
  //bool isEyes() const {return mat.isIdentity(FAUST_PRECISION);}
  
void init_from_file(const char* filename);


  
  
  /// OPERATION BASIQUE ///
  
  //arithmetique
  
  faust_real max() const {return mat.maxCoeff();}
  faust_real min() const {return mat.minCoeff();}
  void abs() {mat=mat.cwiseAbs();}
  
  // return the maximum of all coefficients of this and puts in row_id and col_id its location
  faust_real max(std::vector<faust_unsigned_int> & id_row,std::vector<faust_unsigned_int> & id_col) const;
  faust_real min(std::vector<faust_unsigned_int> & id_row,std::vector<faust_unsigned_int> & id_col) const;
  
  
  // frobenius norm
  faust_real norm() const {return mat.norm();}
  void normalize() {scalarMultiply(1.0/norm());}
  // spectral norm, "norm2", equal to the largest singular value  
  faust_real spectralNorm() const;
  faust_real spectralNorm(const faust_unsigned_int nbr_iter_max,faust_real threshold, faust_int & flag) const;
  
  // trace
  faust_real trace() const {return mat.trace();}
  
  //transposition
  void transpose();
  
  // multiply (*this) = (*this) * A
  void multiplyRight(faust_mat const& A);
  // multiply (*this) =  A * (*this)
  void multiplyLeft(faust_mat const& A);
  
  // scalarMultiply (*this) = (*this) * lambda
  void scalarMultiply(faust_real const lambda);
  // (*this)(i,j)=((*this)(i,j)) * A(i,j)	
  void scalarMultiply(faust_mat const& A);
  // (*this) = (*this) + A
  void add(faust_mat const& A);

  // (*this) = (*this) - A
  void sub(faust_mat const& A);

  
  // Affichage
  void Display() const;
  void print_file(const char* filename)const;

  
  /// SURCHARGE OPERATEUR ///
  // affectation
  void operator=(faust_mat const& A);
  void operator=(faust_spmat const& A);


  void operator-=(faust_mat const& A){sub(A);}
  void operator+=(faust_mat const& A){add(A);}

  void operator*=(faust_mat const& A){multiplyRight(A);}
  

  void operator*=(faust_real lambda){scalarMultiply(lambda);}
  void operator/=(faust_real lambda){scalarMultiply(1.0/lambda);}


  
  ////////////////// friends //////////////////////
  // intra classe//
  friend void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C);
  //friend void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta);
  friend void add(const faust_mat & A, const faust_mat & B, faust_mat & C);
  friend void gemm(const faust_mat & A,const faust_mat & B, faust_mat & C,const faust_real& alpha, const faust_real& beta, char  typeA, char  typeB);
  friend void multiply(const faust_core & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, char typeA, char typeMult);
  friend void gemv(const faust_mat & A,const faust_vec & x,faust_vec & y,const faust_real & alpha, const faust_real & beta, char typeA);
  friend faust_vec solve(const faust_mat & A, const faust_vec & v);
  ///////////friend faust_spmat::operator=(faust_mat const& S);
  bool estIdentite(){return isIdentity;}
  bool estNulle(){return isZeros;}
  
  private: 
       faust_unsigned_int dim1;
       faust_unsigned_int dim2;
  Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat;
  //Eigen::Matrix<faust_real,0,0> mat;
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


 //bool operator==(faust_mat const& A, faust_mat const& B);


#endif
