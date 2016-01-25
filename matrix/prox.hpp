
#include <vector>
#include <iostream>
#include "algorithm"


// const char * interface_prox_name="prox : ";

template<typename T>
inline bool partial_sort_comp (const std::pair<int, T>& pair1, const std::pair<int, T>& pair2) 
{ 
   return fabs(pair1.second) > fabs(pair2.second); 
}

template<typename T>
void sort_idx(const std::vector<T> &v, std::vector<int>& idx, int s) 
{
      std::vector<std::pair<int, T> > vec_pair(v.size());
      for (int i=0 ; i<v.size() ; i++)
          vec_pair[i] = std::make_pair(i,v[i]);
     
      std::partial_sort(vec_pair.begin(), vec_pair.begin()+s, vec_pair.end(),partial_sort_comp<T>);
      idx.resize(s);
      for (int i=0 ; i<s ; i++)
          idx[i]=vec_pair[i].first;
}

template<typename T>
void prox_sp(faust_mat<T> & M,faust_unsigned_int k)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;

	if (k<=0)
		M.setZeros();
	else
	{
		if (k<nb_elt_mat)
		{
		const std::vector<T> vec(M.getData(), M.getData()+nb_elt_mat);
		std::vector<int> index;
		sort_idx(vec, index, k);
		index.erase(index.begin()+k, index.end());
		
		M.setZeros();
		for (int i=0 ; i<index.size() ; i++)
			M.getData()[index[i]] = vec[index[i]];
		}
		
	}
	M.normalize();
}




template<typename T>
void prox_spcol(faust_mat<T> & M,faust_unsigned_int k)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;
	
	if (k<=0)
		M.setZeros();
	else
	{
		if (k<dim1)
		{
			std::vector<std::vector<T> > mat(dim2,std::vector<T>(dim1));
			std::vector<std::vector<int> > index(dim2,std::vector<int>(dim1));
			for (int j=0 ; j < dim2 ; j++)
			{
				mat[j].assign(M.getData()+j*dim1, M.getData()+(j+1)*dim1);
				sort_idx(mat[j], index[j], k); 
				index[j].erase(index[j].begin()+k, index[j].end());
			}
			M.setZeros();
			T*const ptr_data = M.getData();
			for (int j=0 ; j<index.size() ; j++)
				for (int i=0 ; i< index[j].size() ; i++)
					ptr_data[j*dim1+index[j][i]] = mat[j][index[j][i]];
		}	
			
	}
	M.normalize();
}



template<typename T>
void prox_splin(faust_mat<T> & M,faust_unsigned_int k)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;
	if (k<=0)
		M.setZeros();
	else
	{
		if (k<dim2)
		{
			std::vector<std::vector<T> > mat(dim1,std::vector<T>(dim2));
			std::vector<std::vector<int> > index(dim1,std::vector<int>(dim2));
			for (int i=0 ; i < dim1 ; i++)
			{
				for (int j=0 ; j<dim2 ; j++)
					mat[i][j] = M.getData()[j*dim1+i];
				sort_idx(mat[i], index[i], k); 
				index[i].erase(index[i].begin()+k, index[i].end());
			}
			M.setZeros();
			T*const ptr_data = M.getData();
			for (int i=0 ; i<index.size() ; i++)
				for (int j=0 ; j< index[i].size() ; j++)
					ptr_data[(index[i][j])*dim1+i] = mat[i][index[i][j]];
		}	
		
	}
	M.normalize();

}



template<typename T>
void prox_normcol(faust_mat<T> & M,T s)
{

	faust_unsigned_int dim1 = M.getNbRow();
	faust_unsigned_int dim2 = M.getNbCol();
	if (s<0)
	{
		handleError("prox : ","prox_normcol : s < 0 ");	
	}
	
	
	faust_mat<T> current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col;
	std::vector<T> values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);
	T norm_col;
	
	if (s == 0)
	{
		M.setZeros();
	}else
	{
		faust_vec<T> current_col(dim1);
		T scalarMultiply;
		for (int j=0;j<dim2;j++)
		{
			memcpy(current_col.getData(),&(M.getData()[j*dim1]),dim1*sizeof(T));
			scalarMultiply = current_col.norm();
			if (scalarMultiply != 0)
			{	
				scalarMultiply = s/scalarMultiply;
			}
			current_col*= scalarMultiply;
			memcpy(&(M.getData()[j*dim1]),current_col.getData(),dim1*sizeof(T));
		}	
	}
	
}

template<typename T>
void prox_normlin(faust_mat<T> & M,T s)
{
	M.transpose();
	prox_normcol(M,s);
	M.transpose();
	
}

template<typename T>
void prox_sp_pos(faust_mat<T> & M,faust_unsigned_int k)
{
	T*const ptr_data = M.getData();
	//treshold de la matrice
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (ptr_data[i] < 0)
			ptr_data[i]=0;
		
	prox_sp(M,k);
	
}



// M needs to be square and k must divide dimension of M
template<typename T>
void prox_blkdiag(faust_mat<T> & M,int k)
{   
	faust_unsigned_int i,j;
	faust_unsigned_int Msize = M.getNbRow();
	if (Msize != M.getNbCol())
	{
		handleError("prox : ","prox_blkdiag : input matrix must be square");	
	}
	
	int sizeblock = Msize/k;
	faust_mat<T> new_M(Msize,Msize);
	new_M.setZeros();
	std::cout<<"sizeblock : "<<sizeblock<<std::endl;
	
	for (int ii=0;ii<k;ii++)
	{
		for (i=sizeblock*ii;i<sizeblock*(ii+1);i++)
		{
			for (j=sizeblock*ii;j<sizeblock*(ii+1);j++)
			{	
				new_M.getData()[i+j*Msize]=M.getData()[i+j*Msize];	
			}
		
		}
		
	}

	M = new_M;

	T normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}	
		
	
}






template<typename T>
void prox_supp(faust_mat<T> & M,const faust_mat<T> & supp)
{
	if ( (supp.getNbRow() != M.getNbRow()) || (supp.getNbCol() != M.getNbCol()) )
	{
		handleError("prox : ","prox_supp : dimensions of the matrix are not equal");	
	}
	M.scalarMultiply(supp);
	M.normalize();
}

	
	








	
	
	





