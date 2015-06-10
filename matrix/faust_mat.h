#ifndef FAUST_MAT_H
#define FAUST_MAT_H
 
#include <Eigen/Dense>
#include "faust_constant.h"
#include <vector>

class faust_vec;//forward declaration of faust_vec class
class faust_mat
{
public:
  /// Constructeurs ///

    faust_mat(const Eigen::Matrix<faust_real, Eigen::Dynamic,Eigen::Dynamic> & mat_);	
	faust_mat(const faust_real  *mat_,const int nbRow, const int nbCol );
	faust_mat();
	
	
    faust_mat(const int nbRow, const int nbCol);
	
	
  /// GETTEUR SETTEUR ///
  int getNbRow() const;
  int getNbCol() const;
  faust_real getCoeff(const int i,const int j) const;
  void getCoeffs(std::vector<faust_real> & valueS,const std::vector<int> & id_row, const std::vector<int>  & id_col) const;
  void setCoeff(const faust_real & value,const int id_row, const int id_col);
  void setCoeffs(const faust_real value,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  void setCoeffs(const std::vector<faust_real> & valueS,const std::vector<int> & id_row,const std::vector<int>  & id_col);
  
  void resize(const int nbRow,const int nbCol);
  
  // (*this) = la matrice nulle
  void setZeros();
  
  // (*this) = identite, pas forcement carree
  void setEyes();

  inline faust_real* getData(){return mat.data();}
  

  

  
  /// EGALITE ///
  
  bool isZeros() const;
  bool isEqual(const faust_mat & B) const;
  bool isEyes() const;
  
  
  
  /// OPERATION BASIQUE ///
  
  //arithmetique
  
  faust_real max() const;
  void abs();
  
  // return the maximum of all coefficients of this and puts in row_id and col_id its location
  faust_real max(std::vector<int> & id_row,std::vector<int> & id_col) const;
  
  
  // frobenius norm
  faust_real norm() const;
  
  
  // spectral norm, "norm2", equal to the largest singular value  
  faust_real spectralNorm() const;
  
  
  
  // trace
  faust_real trace() const;
  
  
  //transposition
  void transpose();
  
  // multiply (*this) = (*this) * A
  void multiplyRight(faust_mat const& A);
  
  // scalarMultiply (*this) = (*this) * lambda
  void scalarMultiply(faust_real const lambda);
  
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

 
 
 
 ///////////////////// fonction INLINE ///////////////////////////
 
 
//setteur 
inline int faust_mat::getNbRow() const
  {
	  return dim1; 
  }
 
inline int faust_mat::getNbCol() const
  {
	  return dim2;
  }


  
  
inline void faust_mat::setEyes()
{
	mat.setIdentity();
}

inline void faust_mat::setZeros()
{
	mat.setZero();
}  




 
 
 // egalite
inline bool faust_mat::isZeros() const
 {
	return mat.isZero(FAUST_PRECISION);

 }
 
inline bool faust_mat::isEyes() const
{
	return mat.isIdentity(FAUST_PRECISION);
}

 /// OPERATION BASIQUE ///
 

//arithmetique

inline void faust_mat::abs()
{
	mat=mat.cwiseAbs();
}


inline faust_real faust_mat::max() const
{
	return mat.maxCoeff();
} 
 
 
// frobenius norm 
inline faust_real faust_mat::norm() const
 {
	 return mat.norm();
 }

// spectral norm, "norm2", equal to the largest singular value  
inline faust_real faust_mat::spectralNorm() const
{
	return mat.operatorNorm();
}
 
 
 
inline faust_real faust_mat::trace() const
 {
	 return mat.trace();
 }

 

#endif
