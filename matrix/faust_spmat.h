#ifndef __FAUST_SPMAT_H__
#define __FAUST_SPMAT_H__

#include "faust_constant.h"
#include "faust_mat.h"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>

#include "faust_vec.h"
//class faust_vec;
//class faust_mat;

class faust_spmat
{
	public:
		faust_spmat();
		faust_spmat(const faust_spmat& M);
		faust_spmat(const faust_mat& M);
		faust_spmat(const int dim1_, const int dim2_);
		faust_spmat(const int nnz_, const int dim1_, const int dim2_);
		faust_spmat(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<faust_real>& values, const int dim1_, const int dim2_);
		faust_spmat(const Eigen::SparseMatrix<faust_real>& mat_);

		void resize(const int nnz_, const int dim1_, const int dim2_);
		void transpose();

		void operator= (const faust_spmat& M);
		void operator*=(const faust_real alpha);
		void operator/=(const faust_real alpha);
		
		int getNbRow()const{return dim1;}
		int getNbCol()const{return dim2;}
		int getNonZeros()const{return nnz;}
		void init_from_file(const char* filename);
		void Display() const;
		faust_real norm(){return mat.norm();}	

		void print_file(const char* filename)const;
		void init_from_file(char* filename);

		~faust_spmat(){}

	private:
		void update_dim(){dim1=mat.rows();dim2=mat.cols();nnz=mat.nonZeros();}

		

	private:
		Eigen::SparseMatrix<faust_real> mat;	
		int dim1;
		int dim2;
		int nnz;

	friend void faust_mat::operator=(faust_spmat const& S);

	// *this = (*this) * S
	friend void faust_mat::operator*=(const faust_spmat& S);
	// *this = (*this) + S
	friend void faust_mat::operator+=(const faust_spmat& S);
	// *this = (*this) - S
	friend void faust_mat::operator-=(const faust_spmat& S);
	friend void solve(const faust_spmat & A,faust_vec & x, const faust_vec & y);

	// *this = S * (*this) 
	friend void faust_mat::multiplyLeft(const faust_spmat& S);

	// *this = S * (*this) 
	friend void  faust_vec::multiplyLeft(faust_spmat const& A);

	


};



#endif
