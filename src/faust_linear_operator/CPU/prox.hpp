#ifndef __FAUST_PROX_HPP__
#define __FAUST_PROX_HPP__

#include <vector>
#include <iostream>
#include "algorithm"


// const char * interface_prox_name="prox : ";

template<typename FPP>
inline bool partial_sort_comp (const std::pair<int, FPP>& pair1, const std::pair<int, FPP>& pair2)
{
   return fabs(pair1.second) > fabs(pair2.second);
}

template<typename FPP>
void sort_idx(const std::vector<FPP> &v, std::vector<int>& idx, int s)
{
      std::vector<std::pair<int, FPP> > vec_pair(v.size());
      for (int i=0 ; i<v.size() ; i++)
          vec_pair[i] = std::make_pair(i,v[i]);

      std::partial_sort(vec_pair.begin(), vec_pair.begin()+s, vec_pair.end(),partial_sort_comp<FPP>);
      idx.resize(s);
      for (int i=0 ; i<s ; i++)
          idx[i]=vec_pair[i].first;
}

template<typename FPP>
void prox_sp(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
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
		const std::vector<FPP> vec(M.getData(), M.getData()+nb_elt_mat);
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

template<typename FPP>
void prox_spcol(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
{
	prox_spcol_normfree(M,k);
	M.normalize();

}


template<typename FPP>
void prox_spcol_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
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
			std::vector<std::vector<FPP> > mat(dim2,std::vector<FPP>(dim1));
			std::vector<std::vector<int> > index(dim2,std::vector<int>(dim1));
			for (int j=0 ; j < dim2 ; j++)
			{
				mat[j].assign(M.getData()+j*dim1, M.getData()+(j+1)*dim1);
				sort_idx(mat[j], index[j], k);
				index[j].erase(index[j].begin()+k, index[j].end());
			}
			M.setZeros();
			FPP*const ptr_data = M.getData();
			for (int j=0 ; j<index.size() ; j++)
				for (int i=0 ; i< index[j].size() ; i++)
					ptr_data[j*dim1+index[j][i]] = mat[j][index[j][i]];
		}

	}

}

template<typename FPP>
void prox_splin_normfree(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
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
			std::vector<vector<FPP> > mat(dim1,vector<FPP>(dim2));
			std::vector<std::vector<int> > index(dim1,std::vector<int>(dim2));
			for (int i=0 ; i < dim1 ; i++)
			{
				for (int j=0 ; j<dim2 ; j++)
					mat[i][j] = M.getData()[j*dim1+i];
				sort_idx(mat[i], index[i], k);
				index[i].erase(index[i].begin()+k, index[i].end());
			}
			M.setZeros();
			FPP*const ptr_data = M.getData();
			for (int i=0 ; i<index.size() ; i++)
				for (int j=0 ; j< index[i].size() ; j++)
					ptr_data[(index[i][j])*dim1+i] = mat[i][index[i][j]];
		}

	}
}

template<typename FPP>
void prox_splin(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
{
	prox_splin_normfree(M,k);
	M.normalize();

}


template<typename FPP>
void prox_splincol(Faust::MatDense<FPP,Cpu> &M,faust_unsigned_int k)
{
	Faust::MatDense<FPP,Cpu> Mspcol = M;
	Faust::MatDense<FPP,Cpu> Msplin = M;


	prox_spcol_normfree(Mspcol,k);
	prox_splin_normfree(Msplin,k);

	//Msplin.transpose();
	//prox_spcol_normfree(Msplin,k);
	//Msplin.transpose();

	for (int i=0;i<M.getNbCol()*M.getNbRow();i++)
	{
		if (Mspcol(i)!=0)
			Msplin.getData()[i]=0;
	}
	Mspcol+=Msplin;
	Mspcol.normalize();
	M=Mspcol;

}


template<typename FPP>
void prox_normcol(Faust::MatDense<FPP,Cpu> & M,FPP s)
{

	faust_unsigned_int dim1 = M.getNbRow();
	faust_unsigned_int dim2 = M.getNbCol();
	if (s<0)
	{
		handleError("prox : ","prox_normcol : s < 0 ");
	}


	Faust::MatDense<FPP,Cpu> current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col;
	vector<FPP> values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);
	FPP norm_col;

	if (s == 0)
	{
		M.setZeros();
	}else
	{
		Faust::Vect<FPP,Cpu> current_col(dim1);
		FPP scalarMultiply;
		for (int j=0;j<dim2;j++)
		{
			memcpy(current_col.getData(),&(M.getData()[j*dim1]),dim1*sizeof(FPP));
			scalarMultiply = current_col.norm();
			if (scalarMultiply != 0)
			{
				scalarMultiply = s/scalarMultiply;
			}
			current_col*= scalarMultiply;
			memcpy(&(M.getData()[j*dim1]),current_col.getData(),dim1*sizeof(FPP));
		}
	}

}

template<typename FPP>
void prox_normlin(Faust::MatDense<FPP,Cpu> & M,FPP s)
{
	M.transpose();
	prox_normcol(M,s);
	M.transpose();

}

template<typename FPP>
void prox_sp_pos(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k)
{
	FPP*const ptr_data = M.getData();
	//treshold de la matrice
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (ptr_data[i] < 0)
			ptr_data[i]=0;

	prox_sp(M,k);

}



// M needs to be square and k must divide dimension of M
template<typename FPP>
void prox_blkdiag(Faust::MatDense<FPP,Cpu> & M,int k)
{
	faust_unsigned_int i,j;
	faust_unsigned_int Msize = M.getNbRow();
	if (Msize != M.getNbCol())
	{
		handleError("prox : ","prox_blkdiag : input matrix must be square");
	}

	int sizeblock = Msize/k;
	Faust::MatDense<FPP,Cpu> new_M(Msize,Msize);
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

	FPP normM = M.norm();
	if (normM != 0)
	{
		M.scalarMultiply(1/normM);
	}


}



template<typename FPP>
void prox_supp(Faust::MatDense<FPP,Cpu> & M,const Faust::MatDense<FPP,Cpu> & supp)
{
	if ( (supp.getNbRow() != M.getNbRow()) || (supp.getNbCol() != M.getNbCol()) )
	{
		handleError("prox : ","prox_supp : dimensions of the matrix are not equal");
	}
	M.scalarMultiply(supp);
	M.normalize();
}



#endif
