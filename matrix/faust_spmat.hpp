#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;
template<typename T>
const char * faust_spmat<T>::class_name="faust_spmat<T>::";

template<typename T>
faust_spmat<T>::faust_spmat() : 
	faust_mat_generic(),
	mat(Eigen::SparseMatrix<T>(0,0)),
	nnz(0){}

template<typename T>	
faust_spmat<T>::faust_spmat(const faust_spmat<T>& M) :
	faust_mat_generic(M.getNbRow(),M.getNbCol()),
	mat(M.mat),
	nnz(M.mat.nonZeros()){}

template<typename T>
faust_spmat<T>::faust_spmat(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_) :
	faust_mat_generic(dim1_,dim2_),
	mat(Eigen::SparseMatrix<T>(dim1_,dim2_)),
	nnz(0)
{
	resize(nnz, dim1, dim2);
}



template<typename T>
template<typename U>
void faust_spmat<T>::operator=(const faust_spmat<U>& M)
{
resize(M.getNonZeros(),M.getNbRow(),M.getNbCol());	

vector<Eigen::Triplet<T> > tripletList;
   tripletList.reserve(nnz);
   int nbEltIns = 0;
   int nb_elt_rowi;

	for (int i=0;i<M.getNbRow();i++)
	{	
		nb_elt_rowi = M.getOuterIndexPtr()[i+1]-M.getOuterIndexPtr()[i];

		for (int j = 0;j<nb_elt_rowi;j++)
		{	
			tripletList.push_back(Eigen::Triplet<T>((int) i,M.getInnerIndexPtr()[j+nbEltIns], (T) M.getValuePtr()[j+nbEltIns]));
		}
		nbEltIns += nb_elt_rowi;
			
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	makeCompression();
}



template<typename T>
faust_spmat<T>::faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const T* value, const size_t* id_row, const size_t* col_ptr) :
	faust_mat_generic(dim1_,dim2_),
	mat(Eigen::SparseMatrix<T>(dim1_,dim2_)),
	nnz(nnz_)
{	
	vector<Eigen::Triplet<T> > tripletList;
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
			tripletList.push_back(Eigen::Triplet<T>((int) id_row[i+nbEltIns],j, (T) value[i+nbEltIns]));
		}
		nbEltIns += nb_elt_colj;
			
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}


template<typename T>
faust_spmat<T>::faust_spmat(const faust_mat<T>& M) : 
	faust_mat_generic(M.getNbRow(),M.getNbCol()),
	mat(Eigen::SparseMatrix<T>(M.getNbRow(),M.getNbCol())),
	// dim1(M.getNbRow()),
	// dim2(M.getNbCol()),
	nnz(0)
{
   int* rowind = new int[dim1*dim2];
   int* colind = new int[dim1*dim2];
   T* values = new T[dim1*dim2];

   for (int j=0 ; j<dim2 ; j++)
      for (int i=0; i<dim1 ; i++)
         if(M(i,j)!=0.0)
         {
            rowind[nnz] = i;
            colind[nnz] = j;
            values[nnz] = M(i,j);;
	    nnz++; 
         }

   vector<Eigen::Triplet<T> > tripletList;
   tripletList.reserve(nnz);
   for(int i=0 ; i<nnz ; i++)
      tripletList.push_back(Eigen::Triplet<T>(rowind[i], colind[i], values[i]));
   mat.setFromTriplets(tripletList.begin(), tripletList.end());

   delete[] rowind ; rowind=NULL;
   delete[] colind ; colind=NULL;
   delete[] values ; values=NULL;
}





template<typename T>
void faust_spmat<T>::set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const double* value, const size_t* id_row, const size_t* col_ptr) 
{	resize(0,0,0);
	resize(nnz_,dim1_,dim2_);
	vector<Eigen::Triplet<T> > tripletList;
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
			tripletList.push_back(Eigen::Triplet<T>((int) id_row[i+nbEltIns],j,(T) value[i+nbEltIns]));
		}
		nbEltIns += nb_elt_colj;
			
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	nnz = nnz_;
}






template<typename T>
faust_spmat<T>::faust_spmat(const vector<int>& rowidx, const vector<int>& colidx, const vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{
		//cerr << "vectors rowidx, colidx and values have not the same size" << endl;
		//exit(EXIT_FAILURE);
		handleError(class_name,"::constructor : vectors rowidx, colidx and values have not the same size\n");
	}
	
	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}


template<typename T>
faust_spmat<T>::faust_spmat(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_) : 
	faust_mat_generic(dim1_,dim2_),
	mat(Eigen::SparseMatrix<T>(dim1_,dim2_)),
	nnz(nnz_)
{
	resize(nnz, dim1, dim2);
}







template<typename T>
void faust_spmat<T>::init(const vector<int>& rowidx, const vector<int>& colidx, const vector<T>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	if(rowidx.size()!=colidx.size() || rowidx.size()!=values.size())
	{

		handleError(class_name,"init : vectors rowidx, colidx and values have not the same size\n");
	}
	setZeros();
	resize(rowidx.size(), dim1_, dim2_);
	for (int i=0 ; i<rowidx.size() ; i++)
		mat.coeffRef(rowidx[i], colidx[i]) = values[i];
	mat.makeCompressed();
	nnz = mat.nonZeros();
}





template<typename T>
void faust_spmat<T>::Display() const
{
	std::cout<<"dim1 : "<<dim1<<" dim2 : "<<dim2<<" nnz : "<<nnz<<std::endl;
	std::cout<<mat<<std::endl;
	// faust_mat<T> mat_tmp(*this);
	// mat_tmp.Display();
	
}



template<typename T>
void faust_spmat<T>::resize(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	mat.resize(dim1_, dim2_);
	mat.reserve(nnz_);
	update_dim();
	nnz = nnz_;
}


template<typename T>
void faust_spmat<T>::transpose()
{
	Eigen::SparseMatrix<T> mat_tmp = mat.transpose();
	mat = mat_tmp;
	update_dim();
}

template<typename T>
void faust_spmat<T>::operator=(const faust_spmat<T>& M)
{
	mat = M.mat;
	mat.makeCompressed();
	update_dim();
}

template<typename T>
void faust_spmat<T>::operator= (const faust_mat<T>& Mdense)
{
	int nbRow,nbCol,new_nnz;
	nbRow = Mdense.getNbRow();
	nbCol = Mdense.getNbCol();
	mat.resize(nbRow,nbCol);
	//mat = Mdense.mat.sparseView();
	typedef Eigen::Triplet<T> Tip;
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

template<typename T>
void faust_spmat<T>::operator*=(const T alpha)
{
	if (fabs(alpha) == 0.0)
		resize(0, 0, 0);
	else
	{	
		mat *= alpha;
		update_dim();
	}	
}

template<typename T>
void faust_spmat<T>::operator/=(const T alpha)
{

	if(fabs(alpha) == 0.0)
	{
		handleError(class_name,"operator/= : dividing by 0");
	}
	mat /= alpha;
	update_dim();	
}


template<typename T> 
void faust_spmat<T>::print_file(const char* filename)const
{print_file(filename,std::fstream::out);}

template<typename T> 
void faust_spmat<T>::print_file(const char* filename,std::ios_base::openmode mode)const
{
	ofstream fichier;
	fichier.open(filename,mode);
	
	fichier << dim1 << " " << dim2 <<" "<<getNonZeros() << endl;
	for(int i=0 ; i< mat.outerSize() ; i++)
		for(typename Eigen::SparseMatrix<T,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
			fichier << it.row()+1 << " " << it.col()+1 << " " << setprecision(20) << it.value() << endl;

		fichier << endl;

	fichier.close();
}



template<typename T>
void faust_spmat<T>::init_from_file(FILE* fp)
{
	vector<int> row;
	vector<int> col;
	vector<T> val;
	
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
			val.push_back((T)val_tmp);
		}else
		{
			handleError(class_name,"init_from_file : premature end of file");
		}
	}

	if(col.size()!=row.size() || col.size()!=val.size()
		|| dim1_tmp<0 || dim2_tmp<0
		|| *min_element(&row[0],&row[row.size()-1]) <0
		|| *min_element(&col[0],&col[col.size()-1]) <0
		|| *max_element(&row[0],&row[row.size()-1]) > dim1_tmp-1
		|| *max_element(&col[0],&col[col.size()-1]) > dim2_tmp-1)
	{
		handleError(class_name,"init_from_file : Unable to initialize sparse matrix from this file");
	}

	resize(nnz_tmp, dim1_tmp, dim2_tmp);
	vector<Eigen::Triplet<T> > tripletList;
	tripletList.reserve(row.size());
	for(int i=0 ; i<row.size() ; i++)
		tripletList.push_back(Eigen::Triplet<T>(row[i], col[i], val[i]));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	
	mat.makeCompressed();	
}

template<typename T>
void faust_spmat<T>::init_from_file(const char* filename)
{
	// la premiere ligne contient le nombre de lignes et de colonnes de la matrice
	// ainsi que le nombre de nonzeros
	// chacune des autres lignes contient trois valeur par ligne : rowind colind value
	// avec rowind et colind des indices commencant a 1. 

		

	
	FILE* fp=fopen(filename,"r");
	if (fp == NULL)
	{
		handleError(class_name,"init_from_file : unable to open file");	
	}
	init_from_file(fp);
	
	fclose(fp);
	
	
}

