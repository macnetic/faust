/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_PROX_HPP__
#define __FAUST_PROX_HPP__

#include <vector>
#include <iostream>
#include <type_traits>


// const char * interface_prox_name="prox : ";

template<typename FPP>
inline bool Faust::partial_sort_comp (const std::pair<int, FPP>& pair1, const std::pair<int, FPP>& pair2)
{
   return Faust::fabs(pair1.second) > Faust::fabs(pair2.second);
}

template<typename FPP>
void Faust::sort_idx(const std::vector<FPP> &v, std::vector<int>& idx, int s)
{
      std::vector<std::pair<int, FPP> > vec_pair(v.size());
      for (int i=0 ; i<v.size() ; i++)
          vec_pair[i] = std::make_pair(i,v[i]);

      std::partial_sort(vec_pair.begin(), vec_pair.begin()+s, vec_pair.end(),Faust::partial_sort_comp<FPP>);
      idx.resize(s);
      for (int i=0 ; i<s ; i++)
          idx[i]=vec_pair[i].first;
}

template<typename FPP>
void Faust::prox_sp(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized /* true by deft */, const bool pos)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;
	if(pos) Faust::pre_prox_pos(M);
	if (k<=0)
		M.setZeros();
	else
	{
		if (k<nb_elt_mat /* && k < M.getNonZeros()*/)
		{
		const std::vector<FPP> vec(M.getData(), M.getData()+nb_elt_mat);
		std::vector<int> index;
		Faust::sort_idx(vec, index, k);
		index.erase(index.begin()+k, index.end());

		M.setZeros();
		for (int i=0 ; i<index.size() ; i++)
			M.getData()[index[i]] = vec[index[i]];
		}

	}
	if(normalized) M.normalize();
}

template<typename FPP>
void Faust::prox_spcol(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized /* deft to true */, const bool pos)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;

	if(pos) Faust::pre_prox_pos(M);
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
				Faust::sort_idx(mat[j], index[j], k);
				index[j].erase(index[j].begin()+k, index[j].end());
			}
			M.setZeros();
			FPP*const ptr_data = M.getData();
			for (int j=0 ; j<index.size() ; j++)
				for (int i=0 ; i< index[j].size() ; i++)
					ptr_data[j*dim1+index[j][i]] = mat[j][index[j][i]];
		}

	}
	if(normalized) M.normalize();
}

template<typename FPP>
void Faust::prox_splin(Faust::MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized, const bool pos)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	const faust_unsigned_int nb_elt_mat = dim1*dim2;

	if(pos) Faust::pre_prox_pos(M);
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
				Faust::sort_idx(mat[i], index[i], k);
				index[i].erase(index[i].begin()+k, index[i].end());
			}
			M.setZeros();
			FPP*const ptr_data = M.getData();
			for (int i=0 ; i<index.size() ; i++)
				for (int j=0 ; j< index[i].size() ; j++)
					ptr_data[(index[i][j])*dim1+i] = mat[i][index[i][j]];
		}

	}
	if(normalized) M.normalize();
}

template<typename FPP>
void Faust::prox_splincol(Faust::MatDense<FPP,Cpu> &M,faust_unsigned_int k, const bool normalized /* deft to true*/, const bool pos)
{
	Faust::MatDense<FPP,Cpu> Mspcol = M;
	Faust::MatDense<FPP,Cpu> Msplin = M;


	if(pos) Faust::pre_prox_pos(M);
	Faust::prox_spcol(Mspcol,k, false);
	Faust::prox_splin(Msplin,k, false);

	//Msplin.transpose();
	//prox_spcol_normfree(Msplin,k);
	//Msplin.transpose();

	for (int i=0;i<M.getNbCol()*M.getNbRow();i++)
	{
		if (Mspcol(i) != FPP(0))
			Msplin.getData()[i]=0;
	}
	Mspcol+=Msplin;
	M=Mspcol;
	if(normalized) M.normalize();
}

template<typename FPP, typename FPP2>
void Faust::prox_normcol(Faust::MatDense<FPP,Cpu> & M, FPP2 s, const bool normalized /* default to false */, const bool pos)
{

	faust_unsigned_int dim1 = M.getNbRow();
	faust_unsigned_int dim2 = M.getNbCol();

	if(pos) Faust::pre_prox_pos(M);

	if (s < 0)
	{
		handleError("prox : ","Faust::prox_normcol : s < 0 ");
	}


	Faust::MatDense<FPP,Cpu> current_col(dim1,1);
	std::vector<int> id_row,id_col_mat,id_col;
	vector<FPP> values_per_Col;
	id_row.resize(dim1);
	id_col.assign(dim1,0);
	values_per_Col.resize(dim1);

	if (s == 0)
	{
		M.setZeros();
	}else
	{
		Faust::Vect<FPP,Cpu> current_col(dim1);
		FPP2 scalarMultiply;
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

	if(normalized) M.normalize();

}

template<typename FPP, typename FPP2>
void Faust::prox_normlin(Faust::MatDense<FPP,Cpu> & M,FPP2 s, const bool normalized /* default to false */, const bool pos)
{
	M.transpose();
	Faust::prox_normcol<FPP,FPP2>(M,s, normalized, pos);
	M.transpose();
}

template<typename FPP>
void
Faust::pre_prox_pos(MatDense<FPP,Cpu> & M)
{
	FPP*const ptr_data = M.getData();
	//treshold de la matrice
	//	bool is_cplx = typeid(ptr_data[0])==typeid(std::complex<double>())||typeid(ptr_data[0])==typeid(std::complex<float>());
	//don't want to duplicate the function for all realizations of template we need
	//so we use a little trick to make the code valid for double/float and complex<double>/complex<float>
	bool is_cplx = std::is_same<FPP,complex<double>>::value || std::is_same<FPP, complex<float>>::value;
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (!is_cplx && std::complex<float>(ptr_data[i]).real() < 0)
			ptr_data[i]=0;
}


template<typename FPP>
void
Faust::prox_sp_pos(MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized /* default to true*/, const bool pos)
{
	Faust::pre_prox_pos(M);
	Faust::prox_sp(M, k, false /*normfree*/);
	if(normalized) M.normalize();
}

/*
// M needs to be square and k must divide dimension of M
template<typename FPP>
void Faust::prox_blkdiag(Faust::MatDense<FPP,Cpu> & M,int k)
{
	faust_unsigned_int i,j;
	faust_unsigned_int Msize = M.getNbRow();
	if (Msize != M.getNbCol())
	{
		handleError("prox : ","Faust::prox_blkdiag : input matrix must be square");
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
*/


template<typename FPP>
void Faust::prox_supp(Faust::MatDense<FPP,Cpu> & M,const Faust::MatDense<FPP,Cpu> & supp, const bool normalized /* deft to false */, const bool pos)
{
	if(pos) pre_prox_pos(M);
	if ( (supp.getNbRow() != M.getNbRow()) || (supp.getNbCol() != M.getNbCol()) )
	{
		handleError("prox : ","Faust::prox_supp : dimensions of the matrix are not equal");
	}
	M.scalarMultiply(supp);
	if(normalized) M.normalize();
}

template<typename FPP>
void Faust::prox_const(Faust::MatDense<FPP,Cpu> & M,const Faust::MatDense<FPP,Cpu> & const_mat, const bool normalized /* deft to false */, const bool pos /* false */)
{
//	if(pos) pre_prox_pos(M); // useless
	if ( (const_mat.getNbRow() != M.getNbRow()) || (const_mat.getNbCol() != M.getNbCol()) )
	{
		handleError("prox : ","Faust::prox_const_mat : dimensions of the matrix are not equal");
	}
	M = const_mat;
	if(normalized) M.normalize();
}


//template<typename FPP>
//void Faust::prox_nonneg_supp_normfree(Faust::MatDense<FPP,Cpu> & M,const Faust::MatDense<FPP,Cpu> & supp)
//{
//
//}


template<typename FPP> void Faust::prox_toeplitz(Faust::MatDense<FPP, Cpu> & M)
{
	FPP m, m_;
	int dl;
	vector<int> I,J;
	Faust::MatDense<FPP,Cpu> P(M.getNbRow(), M.getNbCol());
	for(int i=0;i<M.getNbRow();i++)
	{
		m = M.diagonal(i).mean();
		m_ = M.diagonal(-i).mean();
		dl = M.getNbRow()-i;
		I.clear();
		J.clear();
		for(int j=0;j<dl;j++)
		{
			I.push_back(j);
			J.push_back(i+j);
		}
		for(int k=0;k<dl;k++)
		{
			*(P.getData()+J[k]*M.getNbRow()+I[k]) = m;
			*(P.getData()+I[k]*M.getNbRow()+J[k]) = m_;
		}
	}
	M = P;
}

template<typename FPP> void Faust::prox_circ(Faust::MatDense<FPP, Cpu> & M)
{
//	cout << "Faust::prox_circ" << endl;
	FPP mi, mj, m;
	int dli, dlj, j;
	vector<int> I, J;
	Faust::MatDense<FPP,Cpu> P(M.getNbRow(), M.getNbCol());
	for(int i=0;i<M.getNbRow();i++)
	{
		j = i - M.getNbRow();
		dli = M.getNbRow()-i;
		dlj = M.getNbRow()+j;
		mi =  M.diagonal(i).mean();
		if(i == 0)
			mj = mi;
		else
			mj = M.diagonal(j).mean();
		m = (FPP(dli)*mi+FPP(dlj)*mj)/(FPP(dli)+FPP(dlj));
//		cout << "m=" << m << "mi=" << mi << "mj=" << mj << "i=" << i << "j=" << j << endl;
		I.clear();
		J.clear();
		//TODO: remove vectors and directly use indices to add elts
		//edit diag i
		for(int k=0;k < dli; k++)
			I.push_back(k);
		for(int k=i;k<dli+i;k++)
			J.push_back(k);
		for(int k = 0; k < dli; k++)
		{
			P.getData()[J[k]*M.getNbRow()+I[k]] = m;
		}
		//edit diag j
		I.clear();
		J.clear();
		for(int k = -j; k < dlj-j; k++)
			I.push_back(k);
		for(int k = 0; k < dlj; k++)
			J.push_back(k);
		for(int k = 0; k <dlj; k++)
		{
			P.getData()[J[k]*M.getNbRow()+I[k]] = m;
		}
	}
	M = P;
	for(int i=0;i<M.getNbRow();i++)
	{
		for(int j=0;j<M.getNbCol();j++)
			cout << M(i,j) << " ";
		cout << endl;
	}
}
#endif
