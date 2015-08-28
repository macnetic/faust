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

  /// Constructeurs ///
  faust_mat(const Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
  faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol );	
  faust_mat() : mat(0,0) , dim1(0) , dim2(0), isIdentity(false),isZeros(false) {}
  faust_mat(const faust_mat & A) : dim1(A.dim1),dim2(A.dim2),mat(A.mat),isIdentity(A.isIdentity),isZeros(A.isZeros) {}
  faust_mat(const faust_spmat & A){this->operator=(A);}

  faust_mat(const int nbRow, const int nbCol) : mat(nbRow,nbCol),dim1(nbRow),dim2(nbCol),isIdentity(false),isZeros(false){}
  faust_mat(const int nbRow) : mat(nbRow,nbRow),dim1(nbRow),dim2(nbRow),isIdentity(false),isZeros(false){}

  faust_mat(const int nbRow, const int nbCol, const faust_real* data_);

	
	
  /// GETTEUR SETTEUR ///
  int getNbRow() const {return dim1;}
  int getNbCol() const {return dim2;}
  faust_real getCoeff(const int i,const int j) const;
  void getCoeffs(std::vector<faust_real> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const;
  void setCoeff(const faust_real & value,const int id_row, const int id_col);
  void setCoeffs(const faust_real value,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  void setCoeffs(const std::vector<faust_real> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  

  void resize(const int nbRow,const int nbCol);
  void resize(const int nbRow){resize(nbRow,nbRow);}
  
  void check_dim_validity();
  
  // (*this) = la matrice nulle
  //void setZeros() {mat.setZero();isZeros=true;}
  void setZeros();
  
  // (*this) = identite, pas forcement carree
  //void setEyes() {mat.setIdentity();if(dim1==dim2)isIdentity=true;}
  void setEyes();

  faust_real& operator[](int i){isZeros=false; isIdentity=false;return mat.data()[i];}

  const faust_real& operator[](int i)const{return mat.data()[i];}

  const faust_real& operator()(int i)const{return mat.data()[i];}
  const faust_real& operator()(int i, int j)const{return mat.data()[j*dim1+i];}


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
void write_into_file(const char* filename);

  
  
  /// OPERATION BASIQUE ///
  
  //arithmetique
  
  faust_real max() const {return mat.maxCoeff();}
  faust_real min() const {return mat.minCoeff();}
  void abs() {mat=mat.cwiseAbs();}
  
  // return the maximum of all coefficients of this and puts in row_id and col_id its location
  faust_real max(std::vector<int> & id_row,std::vector<int> & id_col) const;
  faust_real min(std::vector<int> & id_row,std::vector<int> & id_col) const;
  
  
  // frobenius norm
  faust_real norm() const {return mat.norm();}
  void normalize() {scalarMultiply(1.0/norm());}
  // spectral norm, "norm2", equal to the largest singular value  
  faust_real spectralNorm() const;
  faust_real spectralNorm(const int nbr_iter_max,faust_real threshold, int & flag) const;
  
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
  Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat;
  //Eigen::Matrix<faust_real,0,0> mat;
       int dim1;
       int dim2;
       bool isIdentity;
       bool isZeros;

#ifdef __COMPILE_TIMERS__
  public: 
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
