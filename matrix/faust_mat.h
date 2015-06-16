#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 
#include <Eigen/Dense>
#include "faust_constant.h"
#include <vector>
#include <iterator>

class faust_vec;//forward declaration of faust_vec class
class faust_mat
{
public:

  /// Constructeurs ///
  faust_mat(const Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
  faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol );	
  faust_mat() : mat(0,0) , dim1(0) , dim2(0) {}
  faust_mat(const faust_mat & A) : dim1(A.dim1),dim2(A.dim2),mat(A.mat) {}
  faust_mat(const int nbRow, const int nbCol) : mat(nbRow,nbCol),dim1(nbRow),dim2(nbCol){}
  faust_mat(const int nbRow) : mat(nbRow,nbRow),dim1(nbRow),dim2(nbRow){}

	
	
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
  
  // (*this) = la matrice nulle
  void setZeros() {mat.setZero();}
  
  // (*this) = identite, pas forcement carree
  void setEyes() {mat.setIdentity();}

  faust_real* getData(){return mat.data();}
  const faust_real* getData()const{return mat.data();}  

  
  /// EGALITE ///
  bool isZeros() const {return mat.isZero(FAUST_PRECISION);}
  bool isEqual(const faust_mat & B) const;
  bool isEyes() const {return mat.isIdentity(FAUST_PRECISION);}
  
void init_from_file(const char* filename);
  
  
  /// OPERATION BASIQUE ///
  
  //arithmetique
  
  faust_real max() const {return mat.maxCoeff();}
  void abs() {mat=mat.cwiseAbs();}
  
  // return the maximum of all coefficients of this and puts in row_id and col_id its location
  faust_real max(std::vector<int> & id_row,std::vector<int> & id_col) const;
  
  
  // frobenius norm
  faust_real norm() const {return mat.norm();}
  
  // spectral norm, "norm2", equal to the largest singular value  
  faust_real spectralNorm() const {return mat.operatorNorm();}
  
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
  void operator-=(faust_mat const& A){sub(A);}
  void operator+=(faust_mat const& A){add(A);}

  void operator*=(faust_mat const& A){multiplyRight(A);}

  void operator*=(faust_real lambda){scalarMultiply(lambda);}
  void operator/=(faust_real lambda){scalarMultiply(1.0/lambda);}

  faust_real operator()(int i, int j)const{return mat(i,j);}

  
  ////////////////// friends //////////////////////
  // intra classe//
  friend void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C);
  //friend void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta);
  friend void add(const faust_mat & A, const faust_mat & B, faust_mat & C);
  friend void gemm(const faust_mat & A,const faust_mat & B, faust_mat & C,const faust_real& alpha, const faust_real& beta, char  typeA, char  typeB);
  
  
  private: 
  Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat;
       int dim1;
       int dim2;
  
};


 bool operator==(faust_mat const& A, faust_mat const& B);


#endif
