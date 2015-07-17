#include "faust_spmat.h"
#include "faust_mat.h"
#include "faust_vec.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

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

faust_spmat::faust_spmat(const faust_mat& M) : 
	mat(Eigen::SparseMatrix<faust_real>(M.getNbRow(),M.getNbCol())),
	dim1(M.getNbRow()),
	dim2(M.getNbCol()),
	nnz(0)
{
   int* rowind = new int[dim1*dim2];
   int* colind = new int[dim1*dim2];
   faust_real* values = new faust_real[dim1*dim2];

   for (int j=0 ; j<dim2 ; j++)
      for (int i=0; i<dim1 ; i++)
         if(M(i,j)!=0.0)
         {
            rowind[nnz] = i;
            colind[nnz] = j;
            values[nnz] = M(i,j);;
	    nnz++; 
         }

   vector<Eigen::Triplet<faust_real> > tripletList;
   tripletList.reserve(nnz);
   for(int i=0 ; i<nnz ; i++)
      tripletList.push_back(Eigen::Triplet<faust_real>(rowind[i], colind[i], values[i]));
   mat.setFromTriplets(tripletList.begin(), tripletList.end());

   delete[] rowind ; rowind=NULL;
   delete[] colind ; colind=NULL;
   delete[] values ; values=NULL;
}

faust_spmat::faust_spmat(const vector<int>& rowidx, const vector<int>& colidx, const vector<faust_real>& values, const int dim1_, const int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{
		cerr << "vectors rowidx, colidx and values have not the same size" << endl;
		exit(EXIT_FAILURE);
	}
	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}



void faust_spmat::init_from_file(const char* filename)
{
	FILE* fp=fopen(filename,"r");
	int nb_row,nb_col,_nnz,id_row,id_col;
	faust_real value;
	fscanf(fp, "%d %d %d\n",&nb_row,&nb_col,&_nnz);
	std::cout<<"INSIDE nb_row : "<<nb_row<<" nb_col : "<<nb_col<<" nnz : "<<_nnz<<std::endl;
	resize(_nnz,nb_row,nb_col);
	typedef Eigen::Triplet<double> T;
	std::vector<T> tripletList;
	while (fscanf(fp,"%d %d %lf\n",&id_row,&id_col,&value)!=EOF)
	{
		tripletList.push_back(T(id_row-1,id_col-1,value));
	}
	mat.setFromTriplets(tripletList.begin(),tripletList.end());
}



void faust_spmat::Display() const
{
	faust_mat mat_tmp(*this);
	mat_tmp.Display();
	
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
	nnz = nnz_;
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
		cerr << "Error in faust_spmat::operator/= : dividing by 0" << endl;
		exit(EXIT_FAILURE);
	}
	mat /= alpha;
	update_dim();	
}


 
void faust_spmat::print_file(const char* filename)const
{
	ofstream fichier;
	fichier.open(filename);
	
	for(int i=0 ; i< mat.outerSize() ; i++)
		for(Eigen::SparseMatrix<faust_real>::InnerIterator it(mat,i); it; ++it)
			fichier << it.row()+1 << " " << it.col()+1 << " " << setprecision(20) << it.value() << endl;
	fichier << dim1 << " " << dim2 << " 0.0" << endl;
	fichier.close();
}



void faust_spmat::init_from_file(char* filename)
{
	// la premiere ligne contient le nombre de lignes et de colonnes de la matrice
	// chacune des autres lignes contient trois valeur par ligne : rowind colind value
	// avec rowind et colind des indices commencant a 1. 
	vector<int> row;
	vector<int> col;
	vector<faust_real> val;
	
	int row_tmp;
	int col_tmp;
	double val_tmp;
	
	int dim1_tmp,dim2_tmp;
	
	FILE* fp=fopen(filename,"r");
	if (fp == NULL)
	{
		cerr << "Error in faust_spmat::init_from_file : unable to open \"" << filename << "\"" << endl;
		exit(EXIT_FAILURE);
	}
	
	fscanf(fp,"%d %d\n", &dim1_tmp,&dim2_tmp);
	while(fscanf(fp,"%d %d %lf\n", &row_tmp,&col_tmp,&val_tmp)!=EOF)
	{
		row.push_back(row_tmp - 1);
		col.push_back(col_tmp - 1);
		val.push_back((faust_real)val_tmp);
	}
	fclose(fp);
	
	if(col.size()!=row.size() || col.size()!=val.size()
		|| dim1_tmp<0 || dim2_tmp<0
		|| *min_element(&row[0],&row[row.size()-1]) <0
		|| *min_element(&col[0],&col[col.size()-1]) <0
		|| *max_element(&row[0],&row[row.size()-1]) > dim1_tmp-1
		|| *max_element(&col[0],&col[col.size()-1]) > dim2_tmp-1)
	{
		cerr << "Error in faust_spmat::init_from_file : Unable to initialize sparse matrix from file "<<filename << endl;
		exit(EXIT_FAILURE);
	}

	resize(row.size(), dim1_tmp, dim2_tmp);
	vector<Eigen::Triplet<faust_real> > tripletList;
	tripletList.reserve(row.size());
	for(int i=0 ; i<row.size() ; i++)
		tripletList.push_back(Eigen::Triplet<faust_real>(row[i], col[i], val[i]));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}

