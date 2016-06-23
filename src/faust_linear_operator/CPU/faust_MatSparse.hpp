#ifndef __FAUST_SP_MAT_HPP
#define __FAUST_SP_MAT_HPP

#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;
template<typename FPP>
const char * Faust::MatSparse<FPP,Cpu>::m_className="Faust::MatSparse<FPP,Cpu>::";

template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse() :
	Faust::MatGeneric<FPP,Cpu>(),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(0,0)),
	nnz(0){}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const Faust::MatSparse<FPP,Cpu>& M) :
	Faust::MatGeneric<FPP,Cpu>(M.getNbRow(),M.getNbCol()),
	mat(M.mat),
	nnz(M.mat.nonZeros()){}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_) :
	Faust::MatGeneric<FPP,Cpu>(dim1_,dim2_),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(dim1_,dim2_)),
	nnz(0)
{
	resize(nnz, this->dim1, this->dim2);
}



template<typename FPP>
template<typename FPP1>
void Faust::MatSparse<FPP,Cpu>::operator=(const Faust::MatSparse<FPP1,Cpu>& M)
{
	resize(M.getNonZeros(),M.getNbRow(),M.getNbCol());

	vector<Eigen::Triplet<FPP> > tripletList;
   	tripletList.reserve(nnz);
   	int nbEltIns = 0;
   	int nb_elt_rowi;

	for (int i=0;i<M.getNbRow();i++)
	{
		nb_elt_rowi = M.getOuterIndexPtr()[i+1]-M.getOuterIndexPtr()[i];
		for (int j = 0;j<nb_elt_rowi;j++)
		{
			tripletList.push_back(Eigen::Triplet<FPP>((int) i,M.getInnerIndexPtr()[j+nbEltIns], (FPP) M.getValuePtr()[j+nbEltIns]));
		}
		nbEltIns += nb_elt_rowi;

	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	makeCompression();
}


template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP* value, const size_t* id_row, const size_t* col_ptr) :
	Faust::MatGeneric<FPP,Cpu>(dim1_,dim2_),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(dim1_,dim2_)),
	nnz(nnz_)
{
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz);
	int nbEltIns = 0;
	int nb_elt_colj;
	//std::cout<<"SPMAT CONSTRUCTOR"<<std::endl;
	//std::cout<<"row "<< dim1_<<" col "<<dim2_<<std::endl;
	for (int j=0;j<dim2_;j++)
	{
		nb_elt_colj = col_ptr[j+1]-col_ptr[j];
		//std::cout<<"nb_elt "<< nb_elt_colj<<" col "<<j<<std::endl;
		for (int i = 0;i<nb_elt_colj;i++)
		{
			//std::cout<<"i : "<<id_row[i+nbEltIns]<<" j :"<<j<<" value : "<<value[i+nbEltIns]<<std::endl;
			tripletList.push_back(Eigen::Triplet<FPP>((int) id_row[i+nbEltIns],j, (FPP) value[i+nbEltIns]));
		}
		nbEltIns += nb_elt_colj;
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}


template<typename FPP>
template<typename FPP1>
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col) :
	Faust::MatGeneric<FPP,Cpu>(dim1_,dim2_),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(dim1_,dim2_)),
	nnz(nnz_)
{
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz);
	int nbEltIns = 0;
	int nb_elt_rowi;
	//std::cout<<"SPMAT CONSTRUCTOR"<<std::endl;
	//std::cout<<"row "<< dim1_<<" col "<<dim2_<<std::endl;
	for (int i=0;i<dim1_;i++)
	{
		nb_elt_rowi = row_ptr[i+1]-row_ptr[i];
		//std::cout<<"nb_elt "<< nb_elt_colj<<" col "<<j<<std::endl;
		for (int j = 0;j<nb_elt_rowi;j++)
		{
			//std::cout<<"i : "<<id_row[i+nbEltIns]<<" j :"<<j<<" value : "<<value[i+nbEltIns]<<std::endl;
			tripletList.push_back(Eigen::Triplet<FPP>(i,(int) id_col[j+nbEltIns], (FPP) value[i+nbEltIns]));
		}
		nbEltIns += nb_elt_rowi;
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}













template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const Faust::MatDense<FPP,Cpu>& M) :
	Faust::MatGeneric<FPP,Cpu>(M.getNbRow(),M.getNbCol()),
	mat(Eigen::SparseMatrix<FPP>(M.getNbRow(),M.getNbCol())),
	// dim1(M.getNbRow()),
	// dim2(M.getNbCol()),
	nnz(0)
{
	int* rowind = new int[this->dim1*this->dim2];
	int* colind = new int[this->dim1*this->dim2];
	FPP* values = new FPP[this->dim1*this->dim2];

	for (int j=0 ; j<this->dim2 ; j++)
		for (int i=0; i<this->dim1 ; i++)
         		if(M(i,j)!=0.0)
         		{
            			rowind[nnz] = i;
            			colind[nnz] = j;
            			values[nnz] = M(i,j);;
	    			nnz++;
         		}

	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz);
	for(int i=0 ; i<nnz ; i++)
		tripletList.push_back(Eigen::Triplet<FPP>(rowind[i], colind[i], values[i]));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());

	delete[] rowind ; rowind=NULL;
	delete[] colind ; colind=NULL;
	delete[] values ; values=NULL;
}


template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const double* value, const size_t* id_row, const size_t* col_ptr)
{
	resize(0,0,0);
	resize(nnz_,dim1_,dim2_);
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz_);
	int nbEltIns = 0;
	int nb_elt_colj;
	std::cout<<"SPMAT SET"<<std::endl;
	std::cout<<"row "<< dim1_<<" col "<<dim2_<<std::endl;
	for (int j=0;j<dim2_;j++)
	{
		nb_elt_colj = col_ptr[j+1]-col_ptr[j];
		//std::cout<<"nb_elt "<< nb_elt_colj<<" col "<<j<<std::endl;
		for (int i = 0;i<nb_elt_colj;i++)
		{
			//std::cout<<"i : "<<id_row[i+nbEltIns]<<" j :"<<j<<" value : "<<value[i+nbEltIns]<<std::endl;
			//mat.insert((int)id_row[i+nbEltIns],j)=value[i+nbEltIns];
			tripletList.push_back(Eigen::Triplet<FPP>((int) id_row[i+nbEltIns],j,(FPP) value[i+nbEltIns]));
		}
		nbEltIns += nb_elt_colj;
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	nnz = nnz_;
}



template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const vector<int>& rowidx, const vector<int>& colidx, const vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{
		//cerr << "vectors rowidx, colidx and values have not the same size" << endl;
		//exit(EXIT_FAILURE);
		handleError(m_className,"::constructor : vectors rowidx, colidx and values have not the same size\n");
	}

	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}


template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_) :
	Faust::MatGeneric<FPP,Cpu>(dim1_,dim2_),
	mat(Eigen::SparseMatrix<FPP>(dim1_,dim2_)),
	nnz(nnz_)
{
	resize(nnz,this->dim1, this->dim2);
}





template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::init(const vector<int>& rowidx, const vector<int>& colidx, const vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{

		handleError(m_className,"init : vectors rowidx, colidx and values have not the same size\n");
	}
	setZeros();
	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}





template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::Display() const
{
	std::cout<<"dim1="<<this->dim1<<" ; dim2="<<this->dim2<<" ; nnz="<<nnz<<std::endl;

	cout << "rowPtr = " << getRowPtr() << " -> [ " ;
	for (int i=0 ; i<this->dim1+1 ; i++)
		cout <<  getRowPtr()[i] << " ";
	cout << " ]"<<endl;
	cout << "colInd = " << getColInd() << " -> [ " ;
	for (int i=0 ; i<nnz ; i++)
		cout <<  getColInd()[i] << " ";
	cout << " ]"<<endl;
	cout << "values = " << getValuePtr() << " -> [ " ;
	for (int i=0 ; i<nnz ; i++)
		cout <<  getValuePtr()[i] << " ";
	cout << " ]"<<endl<<endl;


}


template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::display_support() const
{
	for (int i=0;i<this->getNbRow();i++)
	{
		int ind_line_begin=getRowPtr()[i];
		int ind_line_end=getRowPtr()[i+1];

		int nb_elt_per_row=ind_line_end-ind_line_begin;
		int precedent = 0;
		for (int k=0;k<nb_elt_per_row;k++)
		{
			for (int l=precedent+1;l<getColInd()[ind_line_begin+k];l++)
			{
				std::cout<<"O";
			}
			precedent=getColInd()[ind_line_begin+k];
			std::cout<<"X";
		}
		std::cout<<std::endl;


	}
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	mat.resize(dim1_, dim2_);
	mat.reserve(nnz_);
	update_dim();
	nnz = nnz_;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::check_dim_validity() const
{


	if ( (this->getNbCol() != mat.cols()) ||  (this->getNbRow() != mat.rows()))
	{
		cout<<"nb cols attribute : "<<this->getNbCol()<<endl;
		cout<<"nb cols from eigen : "<<mat.cols()<<endl;
		cout<<"nb rows attribute : "<<this->getNbRow()<<endl;
		cout<<"nb rows from eigen : "<<mat.rows()<<endl;
		handleError(m_className, "check_dim_validity : Size incompatibility in the Faust::MatSparse");
	}

	if (this->nnz != mat.nonZeros())
	{
		cout<<"nnz attribute : "<<mat.nonZeros()<<endl;
		cout<<"nnz from eigen : "<<this->nnz<<endl;
		handleError(m_className, "check_dim_validity : incompatibility in the number of non zeros");
	}
}


template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::transpose()
{
	Eigen::SparseMatrix<FPP> mat_tmp = mat.transpose();
	mat = mat_tmp;
	update_dim();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::operator=(const Faust::MatSparse<FPP,Cpu>& M)
{
	mat = M.mat;
	mat.makeCompressed();
	update_dim();
}





template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::operator= (const Faust::MatDense<FPP,Cpu>& Mdense)
{
	int nbRow,nbCol,new_nnz;
	nbRow = Mdense.getNbRow();
	nbCol = Mdense.getNbCol();
	mat.resize(nbRow,nbCol);
	//mat = Mdense.mat.sparseView();
	typedef Eigen::Triplet<FPP> Tip;
	std::vector<Tip> tripletList;

	for (int i=0;i<nbRow*nbCol;i++)
	{
		if (Mdense[i] != 0)
		{
			tripletList.push_back( Tip(i%nbRow,((int) (i/nbRow)),Mdense[i]) );
		}
	}

	mat.setFromTriplets(tripletList.begin(),tripletList.end());
	mat.makeCompressed();
	update_dim();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::operator*=(const FPP alpha)
{
	if (fabs(alpha) == 0.0)
		resize(0, 0, 0);
	else
	{
		mat *= alpha;
		update_dim();
	}
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::operator/=(const FPP alpha)
{

	if(fabs(alpha) == 0.0)
	{
		handleError(m_className,"operator/= : dividing by 0");
	}
	mat /= alpha;
	update_dim();
}


template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::print_file(const char* filename)const
{print_file(filename,std::fstream::out);}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::print_file(const char* filename,std::ios_base::openmode mode)const
{
	ofstream fichier;
	fichier.open(filename,mode);

	fichier << this->dim1 << " " << this->dim2 <<" "<<getNonZeros() << endl;
	for(int i=0 ; i< mat.outerSize() ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
			fichier << it.row()+1 << " " << it.col()+1 << " " << setprecision(20) << it.value() << endl;

		fichier << endl;

	fichier.close();
}



template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::init_from_file(FILE* fp)
{
	vector<int> row;
	vector<int> col;
	vector<FPP> val;

	int row_tmp;
	int col_tmp;
	double val_tmp;

	int dim1_tmp,dim2_tmp,nnz_tmp;

	fscanf(fp,"%d %d %d\n", &dim1_tmp,&dim2_tmp,&nnz_tmp);
	for (int i=0;i<nnz_tmp;i++)
	{
		if (fscanf(fp,"%d %d %lf\n", &row_tmp,&col_tmp,&val_tmp)!=EOF)
		{
			row.push_back(row_tmp - 1);
			col.push_back(col_tmp - 1);
			val.push_back((FPP)val_tmp);
		}else
		{
			handleError(m_className,"init_from_file : premature end of file");
		}
	}

	if(col.size()!=row.size() || col.size()!=val.size()
		|| dim1_tmp<0 || dim2_tmp<0
		|| *min_element(&row[0],&row[row.size()-1]) <0
		|| *min_element(&col[0],&col[col.size()-1]) <0
		|| *max_element(&row[0],&row[row.size()-1]) > dim1_tmp-1
		|| *max_element(&col[0],&col[col.size()-1]) > dim2_tmp-1)
	{
		handleError(m_className,"init_from_file : Unable to initialize sparse matrix from this file");
	}

	resize(nnz_tmp, dim1_tmp, dim2_tmp);
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(row.size());
	for(int i=0 ; i<row.size() ; i++)
		tripletList.push_back(Eigen::Triplet<FPP>(row[i], col[i], val[i]));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());

	mat.makeCompressed();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::init_from_file(const char* filename)
{
	// la premiere ligne contient le nombre de lignes et de colonnes de la matrice
	// ainsi que le nombre de nonzeros
	// chacune des autres lignes contient trois valeur par ligne : rowind colind value
	// avec rowind et colind des indices commencant a 1.
	FILE* fp=fopen(filename,"r");
	if (fp == NULL)
	{
		handleError(m_className,"init_from_file : unable to open file");
	}
	init_from_file(fp);
	fclose(fp);
}


#endif
