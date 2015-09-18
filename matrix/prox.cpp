#include "prox.h"
#include <vector>
#include <iostream>
#include "algorithm"
#include "faust_exception.h"
#include "faust_vec.h"


inline bool partial_sort_comp (const std::pair<int, faust_real>& pair1, const std::pair<int, faust_real>& pair2) 
{ 
   return fabs(pair1.second) > fabs(pair2.second); 
}


void sort_idx(const std::vector<faust_real> &v, std::vector<int>& idx, int s) 
{
      std::vector<std::pair<int, faust_real> > vec_pair(v.size());
      for (int i=0 ; i<v.size() ; i++)
          vec_pair[i] = std::make_pair(i,v[i]);
     
      std::partial_sort(vec_pair.begin(), vec_pair.begin()+s, vec_pair.end(),partial_sort_comp);
      idx.resize(s);
      for (int i=0 ; i<s ; i++)
          idx[i]=vec_pair[i].first;
}

void prox_sp_normfree(faust_mat & M,faust_unsigned_int k)
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
		const std::vector<faust_real> vec(M.getData(), M.getData()+nb_elt_mat);
		std::vector<int> index;
		sort_idx(vec, index, k);
		index.erase(index.begin()+k, index.end());
		
		M.setZeros();
		for (int i=0 ; i<index.size() ; i++)
			M.getData()[index[i]] = vec[index[i]];
		}
		
	}
}

void prox_sp(faust_mat & M, faust_unsigned_int k)
{
	prox_sp_normfree(M,k);
	M.normalize();	
}


void prox_spcol_normfree(faust_mat & M,faust_unsigned_int k)
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
			std::vector<std::vector<faust_real> > mat(dim2,std::vector<faust_real>(dim1));
			std::vector<std::vector<int> > index(dim2,std::vector<int>(dim1));
			for (int j=0 ; j < dim2 ; j++)
			{
				mat[j].assign(M.getData()+j*dim1, M.getData()+(j+1)*dim1);
				sort_idx(mat[j], index[j], k); 
				index[j].erase(index[j].begin()+k, index[j].end());
			}
			M.setZeros();
			faust_real*const ptr_data = M.getData();
			for (int j=0 ; j<index.size() ; j++)
				for (int i=0 ; i< index[j].size() ; i++)
					ptr_data[j*dim1+index[j][i]] = mat[j][index[j][i]];
		}	
			
	}
}

	
void prox_spcol(faust_mat & M,faust_unsigned_int k)
{
	prox_spcol_normfree(M,k);
	M.normalize();		
}

void prox_splin_normfree(faust_mat & M,faust_unsigned_int k)
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
			std::vector<std::vector<faust_real> > mat(dim1,std::vector<faust_real>(dim2));
			std::vector<std::vector<int> > index(dim1,std::vector<int>(dim2));
			for (int i=0 ; i < dim1 ; i++)
			{
				for (int j=0 ; j<dim2 ; j++)
					mat[i][j] = M.getData()[j*dim1+i];
				sort_idx(mat[i], index[i], k); 
				index[i].erase(index[i].begin()+k, index[i].end());
			}
			M.setZeros();
			faust_real*const ptr_data = M.getData();
			for (int i=0 ; i<index.size() ; i++)
				for (int j=0 ; j< index[i].size() ; j++)
					ptr_data[(index[i][j])*dim1+i] = mat[i][index[i][j]];
		}	
		
	}

}

void prox_splin(faust_mat & M,faust_unsigned_int k)
{	
	prox_splin_normfree(M,k);
	M.normalize();	
}


void prox_normcol(faust_mat & M,faust_real s)
{

	faust_unsigned_int dim1 = M.getNbRow();
	faust_unsigned_int dim2 = M.getNbCol();
	if (s<0)
	{
		handleError("prox : prox_normcol : s < 0 ");	
	}
	
	
	faust_mat current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col;
	std::vector<faust_real> values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);
	faust_real norm_col;
	
	if (s == 0)
	{
		M.setZeros();
	}else
	{
		faust_vec current_col(dim1);
		faust_real scalarMultiply;
		for (int j=0;j<dim2;j++)
		{
			memcpy(current_col.getData(),&(M.getData()[j*dim1]),dim1*sizeof(faust_real));
			scalarMultiply = current_col.norm();
			if (scalarMultiply != 0)
			{	
				scalarMultiply = s/scalarMultiply;
			}
			current_col*= scalarMultiply;
			memcpy(&(M.getData()[j*dim1]),current_col.getData(),dim1*sizeof(faust_real));
		}	
	}
	
}

void prox_normlin(faust_mat & M,faust_real s)
{
	M.transpose();
	prox_normcol(M,s);
	M.transpose();
	
}

void prox_sp_pos(faust_mat & M,faust_unsigned_int k)
{
	faust_real*const ptr_data = M.getData();
	//treshold de la matrice
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (ptr_data[i] < 0)
			ptr_data[i]=0;
		
	prox_sp(M,k);
	
}

void prox_sp_pos_normfree(faust_mat & M,faust_unsigned_int k)
{
	faust_real*const ptr_data = M.getData();
	//treshold de la matrice
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (ptr_data[i] < 0)
			ptr_data[i]=0;
		
	prox_sp_normfree(M,k);
	
}

// M needs to be square and k must divide dimension of M
void prox_blkdiag(faust_mat & M,int k)
{   
	faust_unsigned_int i,j;
	faust_unsigned_int Msize = M.getNbRow();
	if (Msize != M.getNbCol())
	{
		handleError(" prox : prox_blkdiag : input matrix must be square");	
	}
	
	int sizeblock = Msize/k;
	faust_mat new_M(Msize,Msize);
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

	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}	
		
	
}
// M needs to be square and k must divide dimension of M




void prox_supp(faust_mat & M, const faust_mat & supp)
{
	
	prox_supp(M,supp);
	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}
	
}


void prox_supp_normfree(faust_mat & M,const faust_mat & supp)
{
	if ( (supp.getNbRow() != M.getNbRow()) || (supp.getNbCol() != M.getNbCol()) )
	{
		handleError("prox : prox_supp : dimensions of the matrix are not equal");	
	}
	M.scalarMultiply(supp);
}

/*
void prox_toeplitz(faust_mat & M, int k)
{	
	int nl = M.getNbRow();
	int nc = M.getNbCol();
	faust_mat crit(nl+nc-1,1);
	faust_mat new_M(nl,nc);
	new_M.setZeros();
	std::vector<faust_real> diag_meanS;
	std::vector<int> num_elt_diagS;
	int num_elt_diag;
	int min_dim = std::min(nl,nc);
	int shift_dim = nc-nl;
	faust_real current_mean;
	int id_col,id_row;
	
	diag_meanS.resize(nl+nc-1);
	num_elt_diagS.resize(nl+nc+1);
	
	if (k>nl+nc-1)
	{
		std::cerr << "ERROR prox_toeplitz : k > (nb_Row+nbCol-1) of M " << std::endl;
		exit( EXIT_FAILURE);	
	}
	
	std::cout<<"boucle critere :"<<std::endl;
	for (int i=-nl+1;i<nc;i++)
	{	

		if (shift_dim < 0)
		{
			if (i<=0)
			{
				if (i<shift_dim)
				{	std::cout<<"cas1"<<std::endl;
					num_elt_diag = min_dim-abs(i+abs(shift_dim));
				}else
				{	
					std::cout<<"cas2"<<std::endl;
					num_elt_diag = min_dim;
				}
			}else
			{
				std::cout<<"cas3"<<std::endl;
				num_elt_diag = min_dim-i;
			}
		}else
		{
			if (i>=0)
			{
				if (i>shift_dim)
				{	std::cout<<"cas4"<<std::endl;
					num_elt_diag = min_dim-abs(i-abs(shift_dim));
				}else
				{	
					std::cout<<"cas5"<<std::endl;
					num_elt_diag = min_dim;
				}
			}else
			{
				std::cout<<"cas6"<<std::endl;
				num_elt_diag = min_dim-abs(i);
			}
		}
		
		
		num_elt_diagS[i+nl-1]=num_elt_diag;
		
		if (i < 0)
		{
			id_row = -i;
			id_col = 0;
		}else
		{
			id_row = 0;
			id_col = i;
		}
		std::cout<<"id_row : "<<id_row<<std::endl;
		std::cout<<"id_col : "<<id_col<<std::endl;
		std::cout<<"nombre diagonal : "<<num_elt_diag<<std::endl;
		
		current_mean = 0;
		std::cout<<"coeff diag"<<std::endl;
		for (int l=0;l<num_elt_diag;l++)
		{	
			std::cout<<M.getCoeff(id_row+l,id_col+l)<<std::endl;
			current_mean = current_mean + M.getCoeff(id_row+l,id_col+l);
		}
		
		current_mean = current_mean/((faust_real) num_elt_diag);
		std::cout<<"mean : "<<current_mean<<std::endl;
		
		diag_meanS[i+nl-1]=current_mean;
		crit.setCoeff(current_mean*current_mean*((faust_real)num_elt_diag),i+nl-1,0);
		
	}
	
	std::cout<<" crit : "<<std::endl;
	crit.Display();
	prox_sp(crit,k);
	int ll;
	std::cout<<"boucle sparse"<<std::endl;
	for (int i=0;i<nl+nc-1;i++)
	{
		if (crit.getCoeff(i,0) != 0)
		{
			ll = i-nl+1;
			if (ll < 0)
			{
				id_row = -ll;
				id_col = 0;
			}else
			{
				id_row = 0;
				id_col = ll;
			}
			current_mean = diag_meanS[i];
			for (int l=0;l<num_elt_diagS[i];l++)
			{
				new_M.setCoeff(current_mean,id_row+l,id_col+l);
			}
		}
		
	}
	
	M = new_M;

	faust_real normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}	
	
	
	
	
}*/
	
	








	
	
	





