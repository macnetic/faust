#include <iostream>
#include "faust_spmat.h"
#include "faust_mat.h"
#include "faust_vec.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>

faust_spmat::faust_spmat() : 
	mat(Eigen::SparseMatrix<faust_real>(0,0)),
	dim1(0),
	dim2(0),
	nnz(0){}

faust_spmat::faust_spmat(const faust_spmat& M) :
	mat(M.mat),
	dim1(M.mat.rows()),
	dim2(M.mat.cols()),
	nnz(M.mat.nonZeros()){}


faust_spmat::faust_spmat(const int dim1_, const int dim2_) : 
	mat(Eigen::SparseMatrix<faust_real>(dim1_,dim2_)),
	dim1(dim1_),
	dim2(dim2_),
	nnz(0)
{
	resize(nnz, dim1, dim2);
}

faust_spmat::faust_spmat(const std::vector<int>& rowidx, const std::vector<int>& colidx, const std::vector<faust_real>& values, const int dim1_, const int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{
		std::cerr << "vectors rowidx, colidx and values have not the same size" << std::endl;
		exit(EXIT_FAILURE);
	}
	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}


faust_spmat::faust_spmat(const int nnz_, const int dim1_, const int dim2_) : 
	mat(Eigen::SparseMatrix<faust_real>(dim1_,dim2_)),
	dim1(dim1_),
	dim2(dim2_),
	nnz(nnz_)
{
	resize(nnz, dim1, dim2);
}

faust_spmat::faust_spmat(const Eigen::SparseMatrix<faust_real>& mat_) : 
	mat(mat_),
	dim1(mat_.rows()),
	dim2(mat_.cols()),
	nnz(mat_.nonZeros()){}

void faust_spmat::resize(const int nnz_, const int dim1_, const int dim2_)
{
	mat.resize(dim1_, dim2_);
	mat.reserve(nnz_);
	update_dim();
}

void faust_spmat::transpose()
{
	Eigen::SparseMatrix<faust_real> mat_tmp = mat.transpose();
	mat = mat_tmp;
	update_dim();
}


void faust_spmat::operator=(const faust_spmat& M)
{
	mat = M.mat;
	mat.makeCompressed();
	update_dim();
}

void faust_spmat::operator*=(const faust_real alpha)
{
	if (fabs(alpha) == 0.0)
		resize(0, 0, 0);
	else
	{	
		mat *= alpha;
		update_dim();
	}	
}
void faust_spmat::operator/=(const faust_real alpha)
{

	if(fabs(alpha) == 0.0)
	{
		std::cerr << "Error in faust_spmat::operator/= : dividing by 0" << std::endl;
		exit(EXIT_FAILURE);
	}
	mat /= alpha;
	update_dim();	
}





