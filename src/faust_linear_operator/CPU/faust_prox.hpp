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
#include <complex>
#include <cmath>

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
	if (k<=0) // it can't be < 0, k is a faust_unsigned_int
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

template<typename FPP> faust_unsigned_int Faust::sparse_size(faust_unsigned_int nnz, faust_unsigned_int nrows)
{
	//TODO: move elsewhere
	auto size = nnz*(sizeof(FPP)+sizeof(int)) + (nrows + 1) * sizeof(int);
	return size;
}

template<typename FPP> faust_unsigned_int Faust::dense_size(faust_unsigned_int nrows, faust_unsigned_int ncols)
{
	//TODO: move elsewhere
	auto size = nrows*ncols*sizeof(FPP);
	return size;
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
			if (scalarMultiply != FPP2(0))
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
void Faust::prox_blockdiag(Faust::MatDense<FPP,Cpu> & M, std::vector<faust_unsigned_int> & m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized /* default to true */, const bool pos)
{
//        if(M.shape != self.shape): raise ValueError('The dimension of the '
//                                                   'projector and matrix must '
//                                                   'agree.')
//        M_ = np.zeros(M.shape)
//        m_ = 0
//        n_ = 0
//        for i,(m,n) in enumerate(zip(self.m_vec, self.n_vec)):
//            print("i=", i, "m=", m, "n=", n)
//            M_[m_:m,n_:n] = M[m_:m,n_:n]
//            m_ = m
//            n_ = n
//        return M_
//
	if(pos) pre_prox_pos(M);
	Faust::MatDense<FPP,Cpu> M_(M.getNbRow(), M.getNbCol());
	M_.setZeros();
	faust_unsigned_int m_ = 0, n_ = 0;
	for(faust_unsigned_int i=0; i< m_vec.size();i++)
	{
		Faust::MatDense<FPP,Cpu> block_mat = M.get_block(m_, n_, m_vec[i]-m_, n_vec[i]-n_);
		M_.set_block(m_, n_, block_mat);
		m_ = m_vec[i];
		n_ = n_vec[i];
	}
	M = M_;
	if(normalized) M.normalize();
}

template<typename FPP>
void Faust::prox_blockdiag(Faust::MatDense<FPP,Cpu> & M, Faust::MatDense<FPP,Cpu> mn_vec, const bool normalized /* default to true */, const bool pos)
{
	std::vector<faust_unsigned_int> m_vec;
	std::vector<faust_unsigned_int> n_vec;
	faust_unsigned_int m, n;

	for(int i=0;i < mn_vec.getNbRow(); i++)
	{
		for(int j=0;j < mn_vec.getNbCol(); j++)
		{
			m = static_cast<faust_unsigned_int>(std::real(std::complex<Real<FPP>>(mn_vec(i,j))));
			n = static_cast<faust_unsigned_int>(std::real(std::complex<Real<FPP>>(mn_vec(i,j))));
			m_vec.push_back(m);
			n_vec.push_back(m);
		}
	}
	Faust::prox_blockdiag(M, m_vec, n_vec, normalized, pos);
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
	bool is_real = std::is_same<FPP, Real<FPP>>::value;
	for (int i=0;i<(M.getNbRow() * M.getNbCol());i++)
		if (is_real && std::complex<float>(ptr_data[i]).real() < 0)
			ptr_data[i]=0;
}

template<typename FPP>
void
Faust::prox_sp_pos(MatDense<FPP,Cpu> & M,faust_unsigned_int k, const bool normalized /* default to true*/, const bool pos)
{
	//TODO: delete this prox, prox_sp with pos argument does the same
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
void Faust::prox_supp(Faust::MatDense<FPP,Cpu> & M,const Faust::MatDense<FPP,Cpu> & supp, const bool normalized /* deft to true */, const bool pos)
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

template<typename FPP> void Faust::prox_hankel(Faust::MatDense<FPP, Cpu> & M, const bool normalized /* default to true*/, const bool pos/* default to false */)
{
//	FPP m, m_;
//	int dl;
//	Faust::MatDense<FPP,Cpu> P(M.getNbRow(), M.getNbCol());
//	for(int i=0;i<M.getNbRow();i++)
//	{
//		m = M.adiagonal(i).mean();
//		m_ = M.adiagonal(-i).mean();
//		dl = M.getNbRow()-i;
//		for(int j=0;j<dl;j++)
//		{
//			*(P.getData()+M.getNbRow()-i-j-1+M.getNbRow()*j) = m_;
//			*(P.getData()+M.getNbRow()-j-1+M.getNbRow()*(j+i)) = m_;
//		}
//	}
//	M = P;
	FPP m;
	if(pos) Faust::pre_prox_pos(M);
	Faust::MatDense<FPP,Cpu> P(M.getNbRow(), M.getNbCol());
	for(int i=1;i<M.getNbCol();i++)
	{
		m = M.adiagonal(-i).mean();
		for(auto ind: M.get_antidiag_indices(-i))
		{
			*(P.getData()+ind.first+M.getNbRow()*ind.second) = m;
		}
	}
	for(int i=0;i<M.getNbRow();i++)
	{
		m = M.adiagonal(i).mean();
		for(auto ind: M.get_antidiag_indices(i))
		{
			*(P.getData()+ind.first+M.getNbRow()*ind.second) = m;
		}
	}
	M = P;
	if(normalized) M.normalize();
}

template<typename FPP> void Faust::prox_toeplitz(Faust::MatDense<FPP, Cpu> & M, const bool normalized /* default to true*/, const bool pos/* default to false */)
{
	FPP m;
	if(pos) Faust::pre_prox_pos(M);
	Faust::MatDense<FPP,Cpu> P(M.getNbRow(), M.getNbCol());
	for(int i=0;i<M.getNbCol();i++)
	{
		m = M.diagonal(i).mean();
		for(auto ind: M.get_diag_indices(i))
		{
			*(P.getData()+ind.first+M.getNbRow()*ind.second) = m;
		}
	}
	for(int i=0;i<M.getNbRow();i++)
	{
		m = M.diagonal(-i).mean();
		for(auto ind: M.get_diag_indices(-i))
		{
			*(P.getData()+ind.first+M.getNbRow()*ind.second) = m;
		}
	}
	M = P;
	if(normalized) M.normalize();
}

template<typename FPP> void Faust::prox_circ(Faust::MatDense<FPP, Cpu> & M, const bool normalized /* default to true*/, const bool pos/* default to false */)
{
	if(M.getNbRow() != M.getNbCol()) throw out_of_range("circulant projector applies only on square matrices");
	FPP mi, mj, m;
	int dli, dlj, j;
	vector<int> I, J;
	if(pos) Faust::pre_prox_pos(M);
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
	if(normalized) M.normalize();
}

	template<typename FPP>
void Faust::prox_skperm(Faust::MatDense<FPP, Cpu> & M_, const unsigned int k, const bool normalized /* default to true*/, const bool pos/* default to false */)
{

	if(pos) Faust::pre_prox_pos(M_);
	//	template<typename FPP> using DenseMat = Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	//	using DenseMatInt = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	// use preproc constants instead of nice C++ template aliases because they are only used here and can't be local to a function
#define DenseMatFPP Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatRealFPP Eigen::Matrix<Real<FPP>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
#define DenseMatInt Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
	// M is a DenseMat<double>&, val is a double
#define set_all_to_self_minus(M, val) \
	for(int i = 0; i < M.rows(); i++) \
	for(int j = 0; j < M.cols(); j++) \
	M(i,j) -= val
	// use a Map to avoid data buffer copying
	Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>> M(M_.getData(), M_.getNbRow(), M_.getNbCol());
	assert(M.rows() == M.cols());
	assert(k > 0);
	unsigned long shape[2] = {(unsigned long)M.rows(), (unsigned long)M.cols()};
	unsigned long size = shape[0];
	int side, vertex, v;
	DenseMatInt matching(M.rows(), M.cols()), visited(2,size);
	DenseMatInt degree(2, size), parent(2, size);
	DenseMatRealFPP potential(2, size), slack(2, size);
	DenseMatRealFPP C(M.rows(), M.cols());

	C = (- M.array().abs().square()).real();

	auto inf = numeric_limits<Real<FPP>>::infinity();
	auto eps = numeric_limits<Real<FPP>>::epsilon();
	eps = 1e-7;
	inf = 1e10;

	matching.setZero();
	potential.setZero();
	potential.block(0,0, 1, size).setConstant(C.minCoeff()); //C.rowwise().minCoeff();
	degree.setZero();
	for(int i=0;i<k*size;i++)
	{
		queue<pair<int,int>> q;
		parent.setConstant(-1);
		visited.setZero();
		slack.setConstant(inf);
		int unmatched_vertex = -1;

		for(int v=0; v < size; v++)
		{
			if(degree(0,v) < k)
			{
				q.push(std::make_pair(0,v));
				visited(0,v) = 1;
			}
		}

		while(true)
		{
			if(q.empty() || q.size() <= 0)
			{
				q = queue<pair<int,int>>();
				Real<FPP> delta = INFINITY;//numeric_limits<Real<FPP>>::infinity();
				for(int side=0;side<2;side++)
				{
					for(int vertex=0; vertex < size; vertex++)
						if(! visited(side, vertex))
							delta = std::min(delta, slack(side, vertex));
				}
				set_all_to_self_minus(slack, delta);
				for(int vertex=0; vertex < size; vertex++)
				{
					if(visited(0,vertex))
						potential(0, vertex) += delta;
					if(visited(1, vertex))
						potential(1, vertex) -= delta;
				}

				for(int side = 0; side < 2; side++)
				{
					for(int vertex = 0; vertex < size; vertex++)
					{
						if(abs(slack(side,vertex)) < eps && ! visited(side, vertex))
						{
							visited(side, vertex) = 1;
							q.push(make_pair(side, vertex));
						}
					}
				}
			}

			side = q.front().first;
			vertex = q.front().second;
			q.pop();

			if(side == 1 && degree(side, vertex) < k)
			{
				unmatched_vertex = vertex;
				if(unmatched_vertex < 0) // simulate python negative index
					unmatched_vertex == degree.cols()+unmatched_vertex;
				degree(1, unmatched_vertex) += 1;
				break;
			}

			DenseMatRealFPP weight;
			DenseMatInt connected;

			if(side == 0)
			{
				weight.resize(1, C.cols());
				weight = C.block(vertex, 0, 1, C.cols());
				connected.resize(1, matching.cols());
				connected = matching.block(vertex, 0, 1, matching.cols());
			}
			else
			{
				weight.resize(C.rows(), 1);
				weight = C.block(0, vertex, C.rows(), 1);
				connected.resize(matching.rows(), 1);
				connected = matching.block(0, vertex, matching.rows(), 1);
				set_all_to_self_minus(connected, 1);
				connected *= -1;
			}

			for(int u=0; u < size; u++)
			{
				Real<FPP> p_diff = weight(u) - potential(side, vertex) - potential(1 - side, u);
				if(std::abs(p_diff) > eps)
				{
					if(side == 0 && ! visited(1 - side, u) && p_diff > 0 && p_diff < (slack(1-side, u)))
					{
						slack(1, u) = p_diff;
						parent(1, u) = vertex;
					}
					if(side == 1 && ! visited(1 - side, u) && p_diff < 0 && - p_diff < (slack(1-side, u)))
					{
						slack(0, u) = - p_diff;
						parent(0, u) = vertex;
					}
					continue;
				}

				if(visited(1 - side, u) || connected(u) == 1)
					continue;

				q.push(make_pair(1 - side, u));
				parent(1 - side, u) = vertex;
				visited(1 - side, u) = 1;

			}


		}

		v = unmatched_vertex;

		while(true)
		{
			if(v < 0)
				v = parent.cols()+v;
			auto u = parent(1, v);
			if(u < 0)
				u = parent.cols()+u;
			auto p = parent(0, u);
			matching(u, v) = 1;
			if(p == -1)
			{
				degree(0, u) += 1;
				break;
			}
			else
			{
				matching(u, p) = 0;
				v = p;
			}
		}
	}
	M = (matching.array() > 0).select(M, 0);
	if(normalized) M_.normalize();
}

template<typename FPP>
void Faust::prox_id(Faust::MatDense<FPP,Cpu> & M, const bool normalized /* deft to false */, const bool pos /* false */)
{
	if(pos) pre_prox_pos(M);
	if(normalized) M.normalize();
}
#endif
