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
#ifndef __FAUST_SP_MAT_HPP
#define __FAUST_SP_MAT_HPP

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cstdlib>
#ifndef NO_MATIO
#include "faust_init_from_matio_mat.h"
#endif
#include "faust_scipy.h"

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
	Faust::MatGeneric<FPP,Cpu>(M.getNbRow(),M.getNbCol(), M.is_ortho, M.is_identity)
{
	*this = M;
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int dim1_, const faust_unsigned_int dim2_) :
	Faust::MatGeneric<FPP,Cpu>(dim1_,dim2_),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(dim1_,dim2_)),
	nnz(0)
{
	resize(nnz, this->dim1, this->dim2);
}


template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::Clone(const bool isOptimize /*default value=false*/) const
{
	if (isOptimize)
	{
		Faust::MatDense<FPP,Cpu> M((*this));
		return optimize(M,(*this));
	} else
	{
		return new Faust::MatSparse<FPP,Cpu>((*this));
	}
}




template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(const Faust::MatDiag<FPP>& D) : Faust::MatGeneric<FPP,Cpu>(D.getNbRow(),D.getNbCol()), mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(D.getNbRow(), D.getNbCol())), nnz(D.getNonZeros())
{
	size_t* id_row = new size_t[nnz];
	size_t* col_ptr = new size_t[this->dim2+1];
	FPP* data = new FPP[nnz];
	FPP c;
	int j = 0;
	col_ptr[0] = 0;
	for(int i=0; i < this->dim2; i++)
		if(i < nnz && (c = D.getData()[i]) != FPP(0))
		{
			id_row[j] = i;
			data[j] = c;
			col_ptr[i+1] = 1+col_ptr[i];
			j++;
		}
		else
			col_ptr[i+1] = col_ptr[i];
	*this = MatSparse<FPP,Cpu>(j, this->dim1, this->dim2, data, id_row, col_ptr);
	delete[] col_ptr;
	delete[] id_row;
	delete[] data;
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
	update_dim();
	this->is_ortho = M.is_ortho;
	this->is_identity = M.is_identity;
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
Faust::MatSparse<FPP,Cpu>::MatSparse(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP1* value, const int* row_ptr, const int* id_col, const bool transpose /*= default to false */) :
	Faust::MatGeneric<FPP,Cpu>(transpose?dim2_:dim1_,transpose?dim1_:dim2_),
	mat(Eigen::SparseMatrix<FPP,Eigen::RowMajor>(transpose?dim2_:dim1_,transpose?dim1_:dim2_)),
	nnz(nnz_)
{
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz);
	int nbEltIns = 0;
	int nb_elt_rowi;
	//std::cout<<"SPMAT CONSTRUCTOR"<<std::endl;
	//std::cout<<"row "<< dim1_<<" col "<<dim2_<<std::endl;
	if(nnz_ > 0)
	{
		for (int i=0;i<dim1_;i++)
		{
			nb_elt_rowi = row_ptr[i+1]-row_ptr[i];
			//		std::cout<<"nb_elt "<< nb_elt_colj<<" col "<<j<<std::endl;
			for (int j = 0;j<nb_elt_rowi;j++)
			{
				//			std::cout<<"i : "<<i <<" j :"<< id_col[j+nbEltIns]<<" value : "<<value[j+nbEltIns]<<std::endl;
				if(transpose)
					tripletList.push_back(Eigen::Triplet<FPP>((int) id_col[j+nbEltIns], i, (FPP) value[j+nbEltIns]));
				else
					tripletList.push_back(Eigen::Triplet<FPP>(i,(int) id_col[j+nbEltIns], (FPP) value[j+nbEltIns]));
			}
			nbEltIns += nb_elt_rowi;
		}
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
	}
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::scalarMultiply(const FPP& lambda)
{
	*this *= lambda;
}


template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::multiply(Faust::MatDense<FPP,Cpu> & M, char opThis) const
{
	faust_unsigned_int nbColOpS,nbRowOpS;
	this->setOp(opThis,nbRowOpS,nbColOpS);

	if(nbColOpS != M.getNbRow())
	{
		handleError(m_className,"multiply : incorrect matrix dimensions\n");
	}

	if (M.is_id())
	{

		M = (*this);
		M.set_id(false);
		M.isZeros = false;

		if (opThis == 'T')
			M.transpose();
		else if (opThis == 'T' || opThis == 'H')
		{
			M.transpose();
			M.conjugate();
		}
	}
	else if (M.isZeros)
	{
		M.resize(nbRowOpS, this->dim2);
		M.setZeros();
	}
	else
	{

		if (opThis == 'N')
			M.mat = this->mat * M.mat;
		else if(opThis == 'T')
			M.mat = this->mat.transpose() * M.mat;
		else if(opThis == 'H')
			M.mat = this->mat.transpose().conjugate() * M.mat;

		M.dim1 = nbRowOpS;
	}


}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::multiply(Faust::MatSparse<FPP,Cpu> & M, char opThis) const
{
	faust_unsigned_int nbColOpS, nbRowOpS;
	this->setOp(opThis,nbRowOpS, nbColOpS);

	if(nbColOpS != M.getNbRow())
	{
		handleError(m_className,"multiply: incorrect matrix dimensions\n");
	}

	if (opThis == 'N')
		M.mat = this->mat * M.mat;
	else if(opThis == 'T')
		M.mat = this->mat.transpose() * M.mat;
	else // if(opThis == 'H')
		M.mat = this->mat.conjugate().transpose() * M.mat;

	M.dim1 = nbRowOpS;
	//M.dim2 doesn't change

	M.nnz = M.mat.nonZeros();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::multiplyRight(Faust::MatSparse<FPP,Cpu> const & M)
{
	this->mat = this->mat*M.mat;
	this->update_dim();
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
	//TODO: opt here. use directly getValuePtr(), getRowPtr(), getColInd() to set triplets
	for (int j=0 ; j<this->dim2 ; j++)
		for (int i=0; i<this->dim1 ; i++)
			if(M(i,j)!=FPP(0.0))
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
	this->set_orthogonal(M.is_ortho);
}


	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::set(const faust_unsigned_int nnz_, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, const FPP* value, const size_t* id_row, const size_t* col_ptr)
{
	resize(0,0,0);
	resize(nnz_,dim1_,dim2_);
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz_);
	int nbEltIns = 0;
	int nb_elt_colj;
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
Faust::MatSparse<FPP,Cpu>::MatSparse(const unsigned int* rowidx, const unsigned int* colidx, const FPP* values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_, faust_unsigned_int nnz)
{
	vector<Eigen::Triplet<FPP> > tripletList;
	tripletList.reserve(nnz);
	for(int i=0 ; i<nnz ; i++)
		tripletList.push_back(Eigen::Triplet<FPP>(rowidx[i], colidx[i], values[i]));
	mat.resize(dim1_, dim2_);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	mat.makeCompressed();
	this->nnz = mat.nonZeros();
	this->update_dim();
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>::MatSparse(std::vector<Eigen::Triplet<FPP>> &tripletList, const faust_unsigned_int dim1, const faust_unsigned_int dim2)
{
	mat.resize(dim1, dim2);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	mat.makeCompressed();
	this->nnz = mat.nonZeros();
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
Faust::MatSparse<FPP,Cpu>::MatSparse(const unsigned int* rowidx, const unsigned int* colidx, const std::vector<FPP>& values, const faust_unsigned_int dim1_, const faust_unsigned_int dim2_)
{
	resize(values.size(), dim1_, dim2_);
	for (int i=0 ; i<values.size() ; i++)
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
std::string Faust::MatSparse<FPP,Cpu>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts /* false by default */) const
{
	std::ostringstream str;
	str << Faust::MatGeneric<FPP, Cpu>::to_string(Sparse, transpose);
	if (displaying_small_mat_elts && this->dim1*this->dim2 < 100)
	{
		str << "rowPtr = " << getRowPtr() << " -> [ " ;
		for (int i=0 ; i<this->dim1+1 ; i++)
			str <<  getRowPtr()[i] << " ";
		str << " ]"<<endl;
		str << "colInd = " << getColInd() << " -> [ " ;
		for (int i=0 ; i<nnz ; i++)
			str <<  getColInd()[i] << " ";
		str << " ]"<<endl;
		str << "values = " << getValuePtr() << " -> [ " ;
		for (int i=0 ; i<nnz ; i++)
			str <<  getValuePtr()[i] << " ";
		str << " ]"<<endl<<endl;

	}
	return str.str();
}

template<typename FPP>
std::string Faust::MatSparse<FPP,Cpu>::to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity)
{
	return Faust::MatGeneric<FPP,Cpu>::to_string(nrows, ncols, transpose, density, nnz, is_identity, Sparse);
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::Display() const
{
	std::cout<<to_string();
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
void Faust::MatSparse<FPP,Cpu>::setEyes()
{
	if(this->getNbRow() == this->getNbCol())
	{
		mat.setIdentity();
		update_dim();
	}
	else
	{
		int nbRow,nbCol,new_nnz;
		nbRow = this->getNbRow();
		nbCol = this->getNbCol();
		new_nnz = nbRow<nbCol?nbRow:nbCol;
		typedef Eigen::Triplet<FPP> Tip;
		std::vector<Tip> tripletList;

		for (int i=0;i<new_nnz;i++)
			tripletList.push_back(Tip(i,i,FPP(1)));

		mat.setFromTriplets(tripletList.begin(),tripletList.end());
		mat.makeCompressed();
		update_dim();
	}
	this->set_id(true);
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::transpose()
{
	Eigen::SparseMatrix<FPP> mat_tmp = mat.transpose();
	mat = mat_tmp;
	update_dim();
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::conjugate(const bool eval /* = true */)
{
	if(eval)
		mat = mat.conjugate().eval();
	else
		mat = mat.conjugate();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::adjoint()
{
	conjugate(false);
	transpose();
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::operator=(const Faust::MatSparse<FPP,Cpu>& M)
{
	mat = M.mat;
	mat.makeCompressed();
	update_dim();

	this->is_ortho = M.is_ortho;
	this->is_identity = M.is_identity;
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
		if (Mdense[i] != FPP(0))
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
	if (alpha == FPP(0.0))
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
void Faust::MatSparse<FPP,Cpu>::setCoeff(const int& i, const int& j, const FPP & val)
{
	if(i > this->getNbRow() || i < 0 || j > this->getNbCol() || j < 0)
		handleError(m_className, "setCoeff() received invalid element indices");
	mat.coeffRef(i,j) = val;
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

template<typename FPP>
matvar_t* Faust::MatSparse<FPP, Cpu>::toMatIOVarDense(bool transpose, bool conjugate) const
{
	matvar_t *var = NULL;
	Faust::MatDense<FPP,Cpu> dense_factor;
	// shouldn't be processed as dense factor, we lose compression of sparse matrix here
	// (see toMatIOVar())
	dense_factor = (*this); //data is copied with operator = redef.
	var = dense_factor.toMatIOVar(transpose, conjugate);
	return var;
}

template<typename FPP>
matvar_t* Faust::MatSparse<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate, const char* var_name/*=nullptr*/) const
{

#ifdef NO_MATIO
	throw std::runtime_error("Sorry but NO_MATIO option was enabled at compiling time, so MAT-IO library wasn't enabled and the matrix can't be saved.");
#else
#if MATIO_MAJOR_VERSION != 1 || MATIO_MINOR_VERSION != 5 // look out MATIO_VERSION_STR for full version string
#error "matio version must be 1.5.x"
#endif
#if MATIO_RELEASE_LEVEL >= 18 // the precise version that switched to unsigned int (1.5.18)
	typedef mat_uint32_t MATIO_INT;
#else
	typedef mat_int32_t MATIO_INT;
#endif
	//TODO: refactor this function because it is a bit too long
	matvar_t* var = NULL;
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat_;
	size_t dims[2];
	mat_sparse_t sparse = {0,};
	//	sparse.njc = (int) this->getNbCol()+1;
	sparse.nzmax = (MATIO_INT) this->nnz;
	sparse.ndata = (MATIO_INT) this->nnz;
	MATIO_INT* jc;
	MATIO_INT* ir = new MATIO_INT[sparse.nzmax];

	Real<FPP>* data;
	mat_complex_split_t z = {0,};
	MATIO_INT nir = 0; //incremented later row by row
	int i = 0;

	matio_types matio_type;
	if(std::is_same<FPP, float>::value)
	{
		matio_type = MAT_T_SINGLE;
	}
	else
	{
		matio_type = MAT_T_DOUBLE;
	}
	int opt = typeid(getValuePtr()[0])==typeid(complex<Real<FPP>>(1.0,1.0))?MAT_F_COMPLEX:0;

	if(opt) {
		z.Re = new Real<FPP>[sparse.nzmax];
		z.Im = new Real<FPP>[sparse.nzmax];
	}
	else
		data = new Real<FPP>[sparse.nzmax];

	if(transpose)
	{
		mat_ = mat.transpose();
		dims[0]=this->getNbCol();
		dims[1]= this->getNbRow();
	}
	else
	{
		mat_ = mat;
		dims[0]=this->getNbRow();
		dims[1]= this->getNbCol();
	}

	sparse.njc = dims[1]+1;
	jc = new MATIO_INT[sparse.njc];

	jc[sparse.njc-1] = this->nnz;
	for(MATIO_INT j=0;j<sparse.njc-1;j++) jc[j] = 0; // TODO: memset/OpenMP ?

	// we use the transpose matrix because we are in row-major order but MatIO wants col-major order
	// and the sparse matrix iterator respects the row-major order
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> st;
	st = mat_.transpose();
	jc[0] = 0; // 1st ele./1st row id of first column, always at id 0 of data/ir
	bool first_col_elt = false;
	MATIO_INT last_col_id = 0;
	for (int k=0; k<st.outerSize(); ++k)
	{
		first_col_elt = true;
		for (typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(st,k); it; ++it)
		{
//			std::cout << "row:" << it.row() << " col:" << it.col() <<  " val:" << it.value() << std::endl;
			// remember we transposed the matrix, so it.row() is the col id, it.col() the row id
			if(it.row() > 0 && first_col_elt)
			{
				jc[last_col_id+1] = nir;
				first_col_elt = false;
				last_col_id = it.row();
			}
			ir[nir] = it.col();
			if(opt)
			{
				((Real<FPP>*)z.Re)[nir] = std::real((complex<Real<FPP>>)it.value());
				if(conjugate)
					((Real<FPP>*)z.Im)[nir] = -std::imag((complex<Real<FPP>>)it.value());
				else
					((Real<FPP>*)z.Im)[nir] = std::imag((complex<Real<FPP>>)it.value());
			}
			else
				data[nir] = std::real(complex<Real<FPP>>(it.value()));
			nir++;
		}
	}
	// count of last nz column
	if(last_col_id > 0) jc[last_col_id+1] = nir;

	// intermediate counts for empty columns
	for(MATIO_INT j=1;j<sparse.njc-1;j++) if(jc[j] == 0) jc[j] = jc[j-1]; // TODO: memset/OpenMP ?

	sparse.ir = ir;
	sparse.jc = jc;
	sparse.nir = nir;
	if(opt) {
		sparse.data = &z;
	}
	else
		sparse.data = data;
	var = Mat_VarCreate(var_name, MAT_C_SPARSE, matio_type, 2, dims, &sparse, opt);
	//	if(var != NULL)
	//		Mat_VarPrint(var,1);
	delete[] jc;
	delete[] ir;
	if(!opt)
		delete[] data;
	return var;
#endif
}
//!  \brief Creates a MatSparse from at matio variable
template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::from_matio_var(matvar_t* var)
{
#ifdef NO_MATIO
	throw std::runtime_error("Sorry but NO_MATIO option was enabled at compiling time, so MAT-IO library wasn't enabled and you can't read this matio variable.");
#else
	init_spmat_from_matvar(*this, var);
#endif
}
//!  \brief Creates a MatSparse from at .mat file
template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::read_from_mat_file(const char *filepath, const char *var_name)
{
#ifdef NO_MATIO
	throw std::runtime_error("Sorry but NO_MATIO option was enabled at compiling time, so MAT-IO library wasn't enabled and you can't read this matio file.");
#else
	init_faust_spmat_from_matio(*this, filepath, var_name);
#endif
}

//!  \brief Saves a MatSparse to a .mat file
template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::save_to_mat_file(const char *filepath, const char *var_name)
{
#ifdef NO_MATIO
	throw std::runtime_error("Sorry but NO_MATIO option was enabled at compiling time, so MAT-IO library wasn't enabled and the matrix can't be saved.");
#else
	//TODO: refactor with MatDense::save_to_mat_file
	int ret;
	matvar_t* matvar = toMatIOVar(false, false, var_name);
	mat_t* matfp = Mat_CreateVer(filepath, NULL, MAT_FT_MAT5);
	if(matfp == NULL)
		handleError("Faust::MatSparse::save_to_mat_file()", "Failed creating file");
//	Mat_VarPrint(matvar, 1);
	ret = Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_NONE); //TODO: enable compression ?
	if(ret)
		handleError("Faust::MatSparse::save_to_mat_file", (std::string("Failed writing the MatSparse to a Matlab file error code: ")+std::to_string(ret)).c_str());
	Mat_VarFree(matvar);
	Mat_Close(matfp);
#endif
}

template<typename FPP>
Real<FPP> Faust::MatSparse<FPP, Cpu>::normL1(faust_unsigned_int& col_id, const bool transpose /*default false*/) const
{
	//TODO: refactor this function with MatDense::normL1()
	faust_unsigned_int i, j, max_j;
	Real<FPP> sum, max_sum;
	Eigen::Matrix<FPP, Eigen::Dynamic,Eigen::Dynamic> vec;
	int dim1, dim2;
	if(transpose){
		dim1 = this->getNbRow();
		dim2 = this->getNbCol();
	}
	{
		dim1 = this->getNbCol();
		dim2 = this->getNbRow();
	}
	for(j=0;j<dim1;j++)
	{
		if(transpose)
			vec=mat.block(j,0,1,this->getNbCol());
		else
			vec=mat.block(0,j,this->getNbRow(),1);
		for(i=0,sum=0;i<dim2;i++)
			sum += std::abs(vec.data()[i]);
		if(j==0 || std::abs(sum) > std::abs(max_sum))
		{
			max_sum = sum;
			max_j = j;
		}
	}
	col_id = max_j;
	return max_sum;
}

template<typename FPP>
Real<FPP> Faust::MatSparse<FPP, Cpu>::normL1(const bool transpose /* default false */) const
{
	faust_unsigned_int id;
	return normL1(id,transpose);
}

template<typename FPP>
Real<FPP> Faust::MatSparse<FPP, Cpu>::normInf(const bool transpose/*=false*/) const
{
	return normL1(!transpose);
}

template<typename FPP>
Real<FPP> Faust::MatSparse<FPP, Cpu>::normInf(faust_unsigned_int& row_id, const bool transpose/*=false*/) const
{
	return normL1(row_id, !transpose);
}

template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, MatDense<FPP, Cpu> & submat) const
{
	if(this->dim1 != row_ids.size() || this->dim2 != col_ids.size())
		submat.resize(row_ids.size(), col_ids.size());
	for(int i=0;i<row_ids.size();i++)
		for(int j=0;j<col_ids.size();j++)
			submat.mat(i,j) = mat.coeff(row_ids[i], col_ids[j]);
}


template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::submatrix(const std::vector<int> &row_ids, const std::vector<int> &col_ids, FPP* submat_data) const
{
	using Map = Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>;
	Map submat_mat(submat_data, row_ids.size(), col_ids.size());
	for(int i=0;i<row_ids.size();i++)
		for(int j=0;j<col_ids.size();j++)
			submat_mat(i,j) = mat.coeff(row_ids[i], col_ids[j]);
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatSparse<FPP,Cpu>::get_col(faust_unsigned_int id) const
{
	Vect<FPP, Cpu> out_vec;
	this->get_col(id, out_vec);
	return out_vec;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::get_col(faust_unsigned_int id, Vect<FPP, Cpu>& out_vec) const
{
	if(id > this->getNbCol())
		handleError("Faust::MatSparse", "Too big column index passed to get_col().");
	if(out_vec.size() != this->getNbRow())
		out_vec.resize(this->getNbRow());
	out_vec.vec = mat.col(id);
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::get_cols(faust_unsigned_int start_col_id, faust_unsigned_int num_cols) const
{
	MatSparse<FPP, Cpu>* subMatrix = new MatSparse<FPP, Cpu>(this->getNbRow(), num_cols);
	get_cols(start_col_id, num_cols, *subMatrix);
	return subMatrix;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::get_cols(const faust_unsigned_int start_col_id, faust_unsigned_int num_cols, MatSparse<FPP, Cpu>& out_cols) const
{
	if(start_col_id + num_cols > this->getNbCol())
		throw std::runtime_error("the column range is not entirely into the matrix dimensions");
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	//	tripletList.reserve((int)(this->getNbRow()*num_cols));
	faust_unsigned_int count = 0;
	for(int i=0 ; i < mat.outerSize() ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
			if(it.col() >= start_col_id && it.col() < start_col_id + num_cols)
			{
				tripletList.push_back(T(it.row(), it.col()-start_col_id, it.value()));
				count++;
			}
	tripletList.resize(count);
//	out_cols.resize(count, this->getNbRow(), num_cols);
	out_cols.mat.setFromTriplets(tripletList.begin(), tripletList.end());
//	out_cols.nnz = out_cols.mat.nonZeros();
	out_cols.update_dim();
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
{
	auto sp_mat = new Faust::MatSparse<FPP, Cpu>(this->getNbRow(), num_cols);
	get_cols(col_ids, num_cols, *sp_mat);
	return sp_mat;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::get_cols(const faust_unsigned_int* orig_col_ids, faust_unsigned_int num_cols, MatSparse<FPP, Cpu>& out_cols) const
{
	//TODO: change the prototype, int for col_ids is consistent with eigen index type
#define SCIPY_GET_COLS
#ifdef SCIPY_GET_COLS // this impl. based on scipy is the fastest according to several benchmarks
					  // https://github.com/scipy/scipy/blob/8a64c938ddf1ae4c02a08d2c5e38daeb8d061d38/scipy/sparse/_compressed.py, function _minor_index_fancy
	auto col_ids = new faust_unsigned_int[num_cols];
	std::iota(col_ids,  col_ids+num_cols, 0);
	std::sort(col_ids, col_ids+num_cols, [orig_col_ids](faust_unsigned_int a, faust_unsigned_int b){return orig_col_ids[a] < orig_col_ids[b];});
	auto col_offsets = new int[this->getNbCol()];
	memset(col_offsets, 0, sizeof(int)*this->getNbCol());
	auto res_indptr = new int[this->getNbRow()+1];
	scipy::csr_column_index1(num_cols, orig_col_ids, this->getNbRow(), this->getNbCol(), getRowPtr(), getColInd(), col_offsets, res_indptr);
	auto nnz = res_indptr[this->getNbRow()];
	auto res_indices = new int[nnz];
	auto res_data = new FPP[nnz];
	scipy::csr_column_index2(col_ids, col_offsets, this->getNonZeros(), getColInd(), getValuePtr(), res_indices, res_data);
	out_cols = MatSparse<FPP, Cpu>(nnz, this->getNbRow(), num_cols, res_data, res_indptr, res_indices);
	delete[] col_ids;
	delete[] col_offsets;
	delete[] res_indptr;
	delete[] res_indices;
	delete[] res_data;
	//TODO: (optimization) initialize out_cols before calling csr_column_index2 and use its buffer directly in csr_column_index2, it would allow to copy only res_indptr
#elif EIGEN_COL_MAJOR_GET_COLS // this implementation is not faster than the two other ones
	Eigen::SparseMatrix<FPP, Eigen::ColMajor> cols_mat(this->getNbRow(), num_cols); // ColMajor to use cols() in write-access (for row-major it is only read-access)
	for(int j=0;j<num_cols;j++)
		cols_mat.col(j) = mat.col(orig_col_ids[j]);
	out_cols.mat = cols_mat;
	out_cols.update_dim();
#else
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	faust_unsigned_int count = 0; // nnz of the out_cols sp mat
	faust_unsigned_int* col_ids = new faust_unsigned_int[num_cols];
	std::map<faust_unsigned_int, std::vector<faust_unsigned_int>> inv_col_ids; // given the col_ids[i] (the key) it returns i (the value)
	// copy the ids because of the sort (don't want to alter orig_col_ids)
	memcpy(col_ids, orig_col_ids, sizeof(faust_unsigned_int)*num_cols);
	// set the map
	for(int i=0;i<num_cols;i++)
		if(col_ids[i] < 0 || col_ids[i] > this->getNbCol())
			throw std::runtime_error("a column index is out of range.");
		else if(inv_col_ids.find(col_ids[i]) != inv_col_ids.end()) // key col_ids[i] is already known
			inv_col_ids[col_ids[i]].push_back(i);
		else // key col_ids[i] is still unknown
			inv_col_ids[col_ids[i]] = std::vector<faust_unsigned_int>(1, i);
	// sort the ids
	std::sort(col_ids, col_ids+num_cols);
	for(int i=0 ; i < mat.outerSize() ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
			if(std::binary_search(col_ids, col_ids+num_cols, it.col()))
			{
				for(auto j: inv_col_ids[it.col()])
				{
					// it.col() is in col_ids, add the entry in the corresponding column (given by the map)
					if(it.value() != FPP(0))
						tripletList.push_back(T(it.row(), j, it.value()));
					count++;
				}
			}
	tripletList.resize(count);
//	out_cols.resize(count, this->getNbRow(), num_cols);
	out_cols.mat.setFromTriplets(tripletList.begin(), tripletList.end());
//	out_cols.nnz = out_cols.mat.nonZeros();
	out_cols.update_dim();
	delete[] col_ids;
#endif
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows, Faust::MatSparse<FPP, Cpu>& out_rows) const
{
	if(start_row_id + num_rows > this->getNbRow())
		throw std::runtime_error("the row range is not entirely into the matrix dimensions");
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	//	tripletList.reserve((int)(this->getNbCol()*num_rows));
	faust_unsigned_int count = 0;
	for(faust_unsigned_int i=start_row_id ; i< start_row_id+num_rows ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
		{
			tripletList.push_back(T(it.row()-start_row_id, it.col(), it.value()));
			count++;
		}
	tripletList.resize(count);
	out_rows.resize(count, num_rows, this->getNbCol());
	out_rows.mat.setFromTriplets(tripletList.begin(), tripletList.end());
	out_rows.nnz = count;
}


template<typename FPP>
Faust::MatSparse<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::get_rows(faust_unsigned_int start_row_id, faust_unsigned_int num_rows) const
{
	MatSparse<FPP, Cpu>* subMatrix = new MatSparse<FPP, Cpu>(num_rows, this->getNbCol());
	get_rows(start_row_id, num_rows, *subMatrix);
	return subMatrix;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows, Faust::MatSparse<FPP, Cpu>& out_rows) const
{
#define SCIPY_GET_ROWS
#ifdef SCIPY_GET_ROWS // the scipy method isn't especially faster (nor slower) than the default one but use it to be consistent with get_cols
	auto res_indptr = new int[num_rows+1];
	// compute number of elements per indexed row
	auto rowptr = getRowPtr();
	for(int i=1;i<=num_rows;i++)
		res_indptr[i] = rowptr[row_ids[i-1]+1] - rowptr[row_ids[i-1]];

	// compute cumulative sum
	res_indptr[0] = 0;
	for(int i=1;i <= num_rows; i++)
		res_indptr[i] += res_indptr[i-1];

	// copy column indices and values for the columns of interest
	auto nnz = res_indptr[num_rows];
	auto res_indices = new int[nnz];
	auto res_data = new FPP[nnz];
	scipy::csr_row_index(num_rows, row_ids, getRowPtr(), getColInd(), getValuePtr(), res_indices, res_data);

	out_rows = MatSparse<FPP, Cpu>(nnz, num_rows, this->getNbCol(), res_data, res_indptr, res_indices);
	delete[] res_indptr;
	delete[] res_indices;
	delete[] res_data;
	//TODO: (optimization) initialize out_rows before calling csr_column_index2 and use its buffer directly in csr_column_index2, it would allow to copy only res_indptr
#else
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	faust_unsigned_int count = 0;
	for(faust_unsigned_int i=0; i < num_rows ; i++)
	{
		if(row_ids[i] < 0 || row_ids[i] > this->getNbRow()) throw std::runtime_error("a row index is out of range.");
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat, row_ids[i]); it; ++it)
		{
			tripletList.push_back(T(i, it.col(), it.value()));
			count++;
		}
	}
	tripletList.resize(count);
	out_rows.resize(count, num_rows, this->getNbCol());
	out_rows.mat.setFromTriplets(tripletList.begin(), tripletList.end());
	out_rows.nnz = count;
#endif
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
{
	MatSparse<FPP,Cpu>* subMatrix = new MatSparse<FPP,Cpu>(num_rows, this->getNbCol());
	get_rows(row_ids, num_rows, *subMatrix);
	return subMatrix;
}

template<typename FPP>
std::vector<int> Faust::MatSparse<FPP, Cpu>::col_nonzero_inds(faust_unsigned_int col_id) const
{
	std::vector<int> ids;
	// if col_id isn't in good range it does nothing (harmless)
	for (int i=0; i < mat.outerSize(); i++)
		for (typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,i); it; ++it)
			if(it.col() == col_id)
			{
				ids.push_back(it.row());
			}
	return ids;
}

template<typename FPP>
std::vector<int> Faust::MatSparse<FPP, Cpu>::row_nonzero_inds(faust_unsigned_int row_id) const
{
	std::vector<int> ids;
	assert(row_id >= 0 && row_id < this->dim1);
	for (typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(mat,row_id); it; ++it)
	{
		ids.push_back(it.col());
	}
	return ids;
}

	template<typename FPP>
Faust::MatSparse<FPP, Cpu>* Faust::MatSparse<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, Real<FPP> density)
{
	std::default_random_engine generator(rand());
	std::uniform_real_distribution<Real<FPP>> distribution(0, 1);
	std::uniform_int_distribution<faust_unsigned_int> int_distribution(0, num_rows*num_cols-1);
	typedef Eigen::Triplet<FPP,faust_unsigned_int> T;
	std::vector<T> tripletList;
	MatSparse<FPP, Cpu>* fsmat = new MatSparse<FPP,Cpu>();
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat(num_rows, num_cols);
	FPP rand_number;
	complex<Real<FPP>> rs;
	try {
		faust_unsigned_int num_elts = (faust_unsigned_int)(num_rows*num_cols*density);
		tripletList.reserve(num_elts);
		faust_unsigned_int i, col_id,row_id, r;
		i = 0;
		while(i < num_elts)
		{
			r = int_distribution(generator);
			row_id = r/num_cols;
			col_id = r - row_id * num_cols;
			assert(col_id >= 0 && col_id < num_cols);
			rs = complex<Real<FPP>>(distribution(generator), distribution(generator));
			memcpy(&rand_number, &rs, sizeof(FPP));
//			rand_number = complex2other<FPP>(rs);
//			cout << row_id << " " << col_id << " " << rand_number << rs << sizeof(FPP) << sizeof(complex<FPP>)<< endl;
			tripletList.push_back(T(row_id,col_id, rand_number));
			i++;
		}
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
		fsmat->mat = mat;
		fsmat->update_dim();
		//	fsmat->Display();
	}
	catch(std::exception e) {
		delete fsmat;
		fsmat = NULL;
		//		std::cerr << "Out of memory." << std::endl;
	}
	return fsmat;
}

	template<typename FPP>
Faust::MatSparse<FPP, Cpu>* Faust::MatSparse<FPP, Cpu>::randMat(faust_unsigned_int num_rows, faust_unsigned_int num_cols, Real<FPP> density, bool per_row)
{
	std::default_random_engine generator(rand());
	std::uniform_real_distribution<Real<FPP>> distribution(0, 1);
	typedef Eigen::Triplet<FPP,faust_unsigned_int> T;
	std::vector<T> tripletList;
	MatSparse<FPP, Cpu>* fsmat = new MatSparse<FPP,Cpu>();
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat(num_rows, num_cols);
	FPP rand_number;
	faust_unsigned_int vec_size, num_vecs, vec_nnz, col_id,row_id, j, k;
	complex<Real<FPP>> rs;
	std::vector<faust_unsigned_int> id_list;
	if(per_row)
	{
		num_vecs = num_rows;
		vec_size = num_cols;
	}
	else
	{
		num_vecs = num_cols;
		vec_size = num_rows;
	}
	vec_nnz = (faust_unsigned_int) (density*vec_size);
	try {
		tripletList.reserve((size_t)(num_vecs*vec_nnz));
		for(faust_unsigned_int i = 0; i < num_vecs; i++)
		{
			if(per_row)
				row_id = i;
			else
				col_id = i;
			//(re)set the indices available for the vec (of weight vec_nnz)
			id_list.resize(0);
			for(j=0;j<vec_size;j++)
				id_list.push_back(j);
			j = vec_nnz;
			while(j > 0 && id_list.size() > 0)
			{
				// set the uniform random law to pick an index among the remaining ones
				std::uniform_int_distribution<faust_unsigned_int> int_distribution(0, id_list.size()-1);
				// pick this index
				k = int_distribution(generator);
				if(per_row)
					col_id = id_list[k];
				else
					row_id = id_list[k];
				// remove the used index (don't want to use it twice to set vector)
				id_list.erase(id_list.begin()+k);
				// random value
				rs = complex<Real<FPP>>(distribution(generator), distribution(generator));
				memcpy(&rand_number, &rs, sizeof(FPP));
				tripletList.push_back(T(row_id, col_id, rand_number));
				j--;
			}
		}
		mat.setFromTriplets(tripletList.begin(), tripletList.end());
		fsmat->mat = mat;
		fsmat->update_dim();
		//	fsmat->Display();
	}
	catch(std::exception e) {
		delete fsmat;
		fsmat = NULL;
		//		std::cerr << "Out of memory." << std::endl;
	}
	return fsmat;
}

	template<typename FPP>
Faust::MatSparse<FPP, Cpu>* Faust::MatSparse<FPP, Cpu>::eye(faust_unsigned_int num_rows, faust_unsigned_int num_cols)
{
	auto eye = new MatSparse<FPP,Cpu>(num_rows, num_cols);
	eye->setEyes();
	return eye;
}

template<typename FPP>
list<pair<int,int>> Faust::MatSparse<FPP,Cpu>::nonzeros_indices(const double& tol/* = 0*/) const
{
	list<pair<int,int>> nz_inds;
	int i, j, k;
	unsigned int rowi_nelts = 0;
	const_cast<Faust::MatSparse<FPP, Cpu>*>(this)->makeCompression(); // can't assume it's already done
	const_cast<Faust::MatSparse<FPP, Cpu>*>(this)->update_dim();
	for(i=0;i<this->dim1;i++)
	{
		rowi_nelts = getOuterIndexPtr()[i+1] - getOuterIndexPtr()[i];
		if(rowi_nelts)
		{
			//non-empty row (nonzero elements)
			for(k=getOuterIndexPtr()[i];k<getOuterIndexPtr()[i+1];k++)
			{
				j = getInnerIndexPtr()[k];
				if(j < this->dim2 /* shouldn't happen */ && std::abs((*this)(i, j)) > tol)
					nz_inds.push_back(make_pair(i,j));
			}
		}
	}
	return nz_inds;
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::delete_col(faust_unsigned_int id)
{
	if(id < 0 || id >= this->getNbCol()) throw std::out_of_range(std::string(m_className)+"::delete_col() index out of bounds");
	Eigen::SparseMatrix<FPP,Eigen::ColMajor> spMat(this->getNbRow(), this->getNbCol()-1);
	if(id > 0)
		spMat.leftCols(id) = this->mat.leftCols(id);
	if(id < this->getNbCol()-1)
		spMat.rightCols(this->getNbCol()-id-1) = this->mat.rightCols(this->getNbCol()-id-1);
	this->mat = spMat;
	this->update_dim();
}

	template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::delete_row(faust_unsigned_int id)
{
	if(id < 0 || id >= this->getNbRow()) throw std::out_of_range(std::string(m_className)+"::delete_row() index out of bounds");
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> spMat(this->getNbRow()-1, this->getNbCol());
	if(id > 0)
		spMat.topRows(id) = this->mat.topRows(id);
	if(id < this->getNbRow()-1)
		spMat.bottomRows(this->getNbRow()-id-1) = this->mat.bottomRows(this->getNbRow()-id-1);
	this->mat = spMat;
	this->update_dim();
}

template<typename FPP>
Faust::MatSparse<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::swap_matrix(faust_unsigned_int order, faust_unsigned_int id1, faust_unsigned_int id2)
{
	unsigned int *colids = new unsigned int[order]; // MS VS refuses the stack because order is not a constexpr, use the heap
	const unsigned int nrows = order;
	const unsigned int ncols = nrows;
	unsigned int *rowptr = new unsigned int[nrows+1];
	int min_id, max_id;
	if(id1<id2)
	{
		min_id = id1;
		max_id = id2;
	}
	else
	{
		min_id = id2;
		max_id = id1;
	}
	std::vector<FPP> values;
	rowptr[0] = 0;
	for(int i=0;i<nrows;i++)
	{
		values.push_back(FPP(1.));
		rowptr[i+1] = rowptr[i]+1;
		colids[i] = i;
	}
	rowptr[nrows] = nrows;
	colids[min_id] = max_id;
	colids[max_id] = min_id;
	auto P = new MatSparse<FPP,Cpu>(rowptr, colids, values, nrows, ncols);
	delete [] colids;
	delete [] rowptr;
	return P;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::swap_rows(faust_unsigned_int id1, faust_unsigned_int id2)
{
	MatDense<FPP,Cpu> dmat(*this);
	dmat.swap_rows(id1, id2);
	*this = dmat;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::swap_cols(faust_unsigned_int id1, faust_unsigned_int id2)
{
	MatDense<FPP,Cpu> dmat(*this);
	dmat.swap_cols(id1, id2);
	*this = dmat;
}

template<typename FPP>
size_t Faust::MatSparse<FPP,Cpu>::getNBytes() const
{
	return this->getNBytes(this->getNonZeros(), this->getNbRow());
}

template<typename FPP>
size_t Faust::MatSparse<FPP,Cpu>::getNBytes(int nnz, int nrows)
{
	return nnz*(sizeof(FPP)+sizeof(int))+(nrows+1)*sizeof(int);
}
template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::copyRowPtr(size_t* out_rowptr) const
{
	auto rowptr = getRowPtr();
	for(int i=0;i<=this->getNbRow();i++)
		out_rowptr[i] = rowptr[i];
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::copyRowPtr(int* out_rowptr) const
{
	memcpy(out_rowptr, getRowPtr(), sizeof(int)*(this->getNbRow()+1));
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::copyColInd(size_t* out_colInd) const
{
	auto colind = getColInd();
	for(int i=0;i<this->getNonZeros();i++)
		out_colInd[i] = colind[i];
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::copyColInd(int* out_colInd) const
{
	memcpy(out_colInd, getColInd(), sizeof(int)*this->getNonZeros());
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::copyValuePtr(FPP* out_values) const
{
	memcpy(out_values, getValuePtr(), sizeof(FPP)*this->getNonZeros());
}

template<typename FPP>
template<typename U>
void Faust::MatSparse<FPP,Cpu>::copyBufs(U* out_rowptr, U* out_colind, FPP* out_values) const
{
	copyValuePtr(out_values);
	copyRowPtr(out_rowptr);
	copyColInd(out_colind);
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::vstack(MatSparse<FPP, Cpu>& top, MatSparse<FPP, Cpu>& bottom)
{
	auto tncols = top.getNbCol();
	auto bncols = bottom.getNbCol();
	auto tnrows = top.getNbRow();
	auto bnrows = bottom.getNbRow();
	auto tnnz = top.getNonZeros();
	auto bnnz = bottom.getNonZeros();
	auto nrows = tnrows+bnrows;
	auto ncols = tncols;
	auto nnz = tnnz + bnnz;
	if(tncols != bncols)
		throw std::runtime_error("vstack error: dimensions must agree.");
	if(this->getNbCol() != ncols || this->getNbRow() != nrows || this->getNonZeros() != nnz)
		resize(nnz, nrows, ncols);
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	tripletList.reserve(tnnz+bnnz);
	for(faust_unsigned_int i=0; i < top.getNbRow() ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(top.mat,i); it; ++it)
		{
			tripletList.push_back(T(it.row(), it.col(), it.value()));
		}

	for(faust_unsigned_int i=0; i < bottom.getNbRow() ; i++)
		for(typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(bottom.mat,i); it; ++it)
		{
			tripletList.push_back(T(it.row()+top.getNbRow(), it.col(), it.value()));
		}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	update_dim();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::hstack(MatSparse<FPP, Cpu>& left, MatSparse<FPP, Cpu>& right)
{
	auto lncols = left.getNbCol();
	auto rncols = right.getNbCol();
	auto lnrows = left.getNbRow();
	auto rnrows = right.getNbRow();
	auto lnnz = left.getNonZeros();
	auto rnnz = right.getNonZeros();
	auto nrows = lnrows;
	auto ncols = lncols + rncols;
	auto nnz = lnnz + rnnz;
	int lrow_count, rrow_count, tot_count;
	int *rowptr, *lrowptr = left.getRowPtr(), *rrowptr = right.getRowPtr();
	if(lnrows != rnrows)
		throw std::runtime_error("hstack error: dimensions must agree.");
	if(this->getNbCol() != ncols || this->getNbRow() != nrows || this->getNonZeros() != nnz)
		resize(nnz, nrows, ncols);
	tot_count = 0;
	int i;
	//TODO : openmp ?
//#pragma omp parallel for
	for(i=0;i<nrows;i++)
	{
		rrow_count = rrowptr[i+1]-rrowptr[i];
		lrow_count = lrowptr[i+1]-lrowptr[i];
		memcpy(getValuePtr()+tot_count, left.getValuePtr()+lrowptr[i], sizeof(FPP)*lrow_count);
		memcpy(getValuePtr()+tot_count+lrow_count, right.getValuePtr()+rrowptr[i], sizeof(FPP)*rrow_count);
		memcpy(this->getColInd()+tot_count, left.getColInd()+lrowptr[i], sizeof(int)*lrow_count);
		for(int j=0;j<rrow_count;j++)
		{
			this->getColInd()[tot_count+lrow_count+j] = right.getColInd()[rrowptr[i]+j]+lncols;
		}
		getRowPtr()[i] = tot_count;
		tot_count += lrow_count + rrow_count;
	}
	getRowPtr()[i] = tot_count;
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::print_bufs(const std::string name/*=""*/)
{
	std::cout << "======" << std::endl;
	std::cout << name << " data: ";
	for(int i=0;i<this->getNonZeros();i++)
		std::cout << *(this->getValuePtr()+i) << " ";
	std::cout <<	std::endl;
	std::cout << name << " colind: ";
	for(int i=0;i<this->getNonZeros();i++)
		std::cout << *(this->getColInd()+i) << " ";
	std::cout << std::endl;
	std::cout << name << " rowptr: ";
	for(int i=0;i<this->getNbRow()+1;i++)
		std::cout << *(this->getRowPtr()+i) << " ";
	std::cout << std::endl << "======" << std::endl;
}

template<typename FPP>
void Faust::MatSparse<FPP, Cpu>::print_asarray(const std::string name/*=""*/)
{
	MatDense<FPP,Cpu> A(*this);
	std::cout << "======" << std::endl;
	std::cout << name << ": " << std::endl;
	for(int i=0;i <this->getNbRow();i++)
	{
		for(int j=0;j < this->getNbCol();j++)
			std:cout << *(A.getData()+j*A.getNbRow()+i) << " ";
		std::cout << std::endl;
	}
	std::cout << "======" << std::endl;
}

	template<typename FPP>
void Faust::copy_sp_mat(Faust::MatSparse<FPP,Cpu>& src, Faust::MatSparse<FPP, Cpu>& dst)
{
	dst.resize(src.getNonZeros(),src.getNbRow(),src.getNbCol());
	memcpy(dst.getValuePtr(), src.getValuePtr(), src.getNonZeros()*sizeof(FPP));
	memcpy(dst.getColInd(), src.getColInd(), src.getNonZeros()*sizeof(int));
	memcpy(dst.getRowPtr(), src.getRowPtr(), (src.getNbRow()+1)*sizeof(int));
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::real(MatSparse<Real<FPP>, Cpu> &real_mat) const
{
	real_mat.resize(this->nnz, this->getNbRow(), this->getNbCol());
	real_mat.mat = mat.real().eval().template cast<Real<FPP>>();
}

template<typename FPP>
template<typename FPP2>
Faust::MatSparse<Real<FPP2>, Cpu> Faust::MatSparse<FPP, Cpu>::to_real() const
{
	Faust::MatSparse<FPP2, Cpu> smat;
	smat.resize(this->getNonZeros(), this->getNbRow(), this->getNbCol());
	smat.mat = mat.real().eval().template cast<Real<FPP2>>();
	return smat;
}


template<typename FPP>
template<typename FPP2>
Faust::MatSparse<FPP2, Cpu> Faust::MatSparse<FPP, Cpu>::cast() const
{
	Faust::MatSparse<FPP2, Cpu> smat;
	smat.resize(this->getNonZeros(), this->getNbRow(), this->getNbCol());
	smat.mat = mat.eval().template cast<FPP2>();
	return smat;
}

template<typename FPP>
bool Faust::MatSparse<FPP,Cpu>::containsNaN() const
{
	return Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>(const_cast<FPP*>(getValuePtr()) /* no worry, just a read access */, getNonZeros(), 1).hasNaN();
}

namespace Faust {
	template<typename FPP>
		template<typename MatType1, typename MatType2>
		void MatSparse<FPP, Cpu>::eigenIndexMul(const faust_unsigned_int* row_ids, const faust_unsigned_int* col_ids, size_t nrows, size_t ncols, const MatType1 &in_mat, MatType2 &out_mat, bool transpose/* = false*/, bool conjugate /*= false*/)
		{

#if (EIGEN_WORLD_VERSION >= 3 && EIGEN_MAJOR_VERSION >= 4)

			if(row_ids == nullptr && col_ids == nullptr)
			{
				if(transpose && conjugate)
					out_mat = mat.adjoint() * in_mat;
				else if(transpose)
					out_mat = mat.adjoint() * in_mat;
				else if(conjugate)
					out_mat = mat.conjugate() * in_mat;
				else
					out_mat = mat * in_mat;
			}
			else if(row_ids == nullptr && col_ids != nullptr)
			{
				out_mat = MatType1::Zero(transpose?mat.cols():mat.rows(), in_mat.cols());
				for(int i=0;i < ncols; i++)
				{
					if(transpose && conjugate)
						out_mat += mat.row(col_ids[i]).adjoint() * in_mat.row(i);
					else if(transpose)
					{
						out_mat += mat.row(col_ids[i]).transpose() * in_mat.row(i);
					}
					else if(conjugate)
						out_mat += mat.col(col_ids[i]).conjugate() * in_mat.row(i);
					else
						out_mat += mat.col(col_ids[i]) * in_mat.row(i);
				}
			}
			else if(row_ids != nullptr && col_ids == nullptr)
			{
				SpMat tmp(nrows,transpose?mat.rows():mat.cols());
				//#pragma omp parallel for // can't use OpenMP because even if the rows are independent the whole rowptr buffer is modified when a row is
				for(int i=0;i<nrows;i++)
					if(transpose && conjugate)
						tmp.row(i) = mat.col(row_ids[i]).adjoint();
					else if(transpose)
						tmp.row(i) = mat.col(row_ids[i]).transpose();
					else if(conjugate)
						tmp.row(i) = mat.row(row_ids[i]).conjugate();
					else
						tmp.row(i) = mat.row(row_ids[i]);
				out_mat = tmp * in_mat;
			}
			else // if(row_ids != nullptr && col_ids != nullptr)
			{
				SpMat tmp(nrows, ncols);
				for(int i = 0; i < nrows; i++)
				{
					for(int j = 0; j < ncols; j++)
					{
					  if(transpose)
						 tmp.coeffRef(i, j) = mat.coeff(col_ids[j], row_ids[i]);
					  else
						  tmp.coeffRef(i, j) = mat.coeff(row_ids[i], col_ids[j]);
					}
				}
				if(conjugate)
					tmp = tmp.conjugate();
				out_mat = tmp * in_mat;
			}
#else
			throw std::runtime_error("MatSparse::eigenIndexMul is not supported with eigen version < 3.9");
#endif
		}


	template <typename FPP>
		MatDense<FPP, Cpu> MatSparse<FPP, Cpu>::to_dense() const
		{
			return MatDense<FPP, Cpu>(*this);
		}
}
#endif
