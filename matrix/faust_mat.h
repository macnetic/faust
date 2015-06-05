#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 
#include <Eigen/Dense>
#include "faust_constant.h"


class faust_vec;//forward declaration of faust_vec class
class faust_mat
{
public:
  /// Constructeurs ///

    faust_mat(const Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
	faust_mat(const faust_real  *mat_,const int & nbRow, const int & nbCol );
	
	
    faust_mat(const int & nbRow, const int & nbCol);
	
	
  /// GETTEUR SETTEUR ///
  int getNbRow() const;
  int getNbCol() const;
  
  void resize(const int & nbRow,const int & nbCol);
  
  // (*this) = la matrice nulle
  void setZeros();
  
  // (*this) = identite 
  void setEyes();


  
  /// EGALITE ///
  
  bool isZeros() const;
  bool isEqual(const faust_mat & B) const;
  bool isEyes() const;
  
  
  
  /// OPERATION BASIQUE ///
  
  // frobenius norm
  faust_real norm() const;
  
  // trace
  faust_real trace() const;
  
  
  //transposition
  void transpose();
  
  // multiply (*this) = (*this) * A
  void multiplyRight(faust_mat const& A);
  
  // scalarMultiply (*this) = (*this) * lambda
  void scalarMultiply(faust_real const& lambda);
  
  // (*this) = (*this) + A
  void add(faust_mat const& A);
  
  
  // Affichage
  void Display() const;

  
  /// SURCHARGE OPERATEUR ///
  // affectation
  void operator=(faust_mat const& A);

  
  
  ////////////////// friends //////////////////////
  // intra classe//
  friend void multiply(const faust_mat & A, const faust_mat & B, faust_mat & C);
  friend void gemm(const faust_mat & A, const faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta);
  friend void add(const faust_mat & A, const faust_mat & B, faust_mat & C);
  friend void gemm(faust_mat & A, faust_mat & B, faust_mat & C,const faust_real & alpha, const faust_real & beta, char  typeA, char  typeB);

 
  
  private: 
  Eigen::Matrix<faust_real, Eigen::Dynamic, Eigen::Dynamic> mat;
       int dim1;
       int dim2;
  

};


 bool operator==(faust_mat const& A, faust_mat const& B);


#endif