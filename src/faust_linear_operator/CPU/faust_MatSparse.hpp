/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
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
Faust::MatGeneric<FPP,Cpu>* Faust::MatSparse<FPP,Cpu>::Clone(const bool isOptimize /*default value=false*/) const
{
	if (isOptimize)
	{	
		Faust::MatDense<FPP,Cpu> M((*this));
		return optimize(M,(*this));
	}else
	{
		return new Faust::MatSparse<FPP,Cpu>((*this));
	}
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
			tripletList.push_back(Eigen::Triplet<FPP>(i,(int) id_col[j+nbEltIns], (FPP) value[j+nbEltIns]));
		}
		nbEltIns += nb_elt_rowi;
	}
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
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

	if (M.isIdentity)
	{
				
		M = (*this);
		M.isIdentity = false;
		M.isZeros = false;

		if (opThis == 'T')
			M.transpose();
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
			else
				M.mat = this->mat.transpose() * M.mat;

			M.dim1 = nbRowOpS;
	}


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
std::string Faust::MatSparse<FPP,Cpu>::to_string() const
{
	std::ostringstream str;
	str<<" type : SPARSE";
	str << Faust::MatGeneric<FPP,Cpu>::to_string();
	if (this->dim1*this->dim2 < 100)
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
void Faust::MatSparse<FPP,Cpu>::transpose()
{
	Eigen::SparseMatrix<FPP> mat_tmp = mat.transpose();
	mat = mat_tmp;
	update_dim();
}

template<typename FPP>
void Faust::MatSparse<FPP,Cpu>::conjugate()
{
	mat = mat.conjugate().eval();
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
matvar_t* Faust::MatSparse<FPP, Cpu>::toMatIOVar(bool transpose, bool conjugate) const {
	//TODO: refactor this function because it is a bit too long
	matvar_t* var = NULL;
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat_;
	size_t dims[2];
	mat_sparse_t sparse = {0,};
	//	sparse.njc = (int) this->getNbCol()+1;
	sparse.nzmax = (int) this->nnz;
	sparse.ndata = (int) this->nnz;
	int* jc;
	int* ir = new int[sparse.nzmax];
	double* data;
	mat_complex_split_t z = {0,};
	int nir = 0; //incremented later row by row
	int i = 0;

	int opt = typeid(getValuePtr()[0])==typeid(complex<double>(1.0,1.0))?MAT_F_COMPLEX:0;

	if(opt) {
		z.Re = new double[sparse.nzmax];
		z.Im = new double[sparse.nzmax];
	}
	else
		data = new double[sparse.nzmax];

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
	jc = new int[sparse.njc];

	jc[sparse.njc-1] = this->nnz;
	for(int j=0;j<sparse.njc-1;j++) jc[j] = -1;

	// we use the transpose matrix because we are in row-major order but MatIO wants col-major order
	// and the sparse matrix iterator respects the row-major order
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> st;
	st = mat_.transpose();
	for (int k=0; k<st.outerSize(); ++k)
		for (typename Eigen::SparseMatrix<FPP,Eigen::RowMajor>::InnerIterator it(st,k); it; ++it)
		{
//			std::cout << "row:" << it.row() << " col:" << it.col() <<  " val:" << it.value() << std::endl;
			if(it.row() > 0){
				i=1;
				while(it.row()>=i && jc[it.row()-i] < 0) {
					jc[it.row()-i] = nir;
					i++;
				}
			}
			if(jc[it.row()] < 0) jc[it.row()] = nir;
			ir[nir] = it.col();
			if(opt) {
				((double*)z.Re)[nir] = std::real((complex<double>)it.value());
				if(conjugate)
					((double*)z.Im)[nir] = -std::imag((complex<double>)it.value());
				else
					((double*)z.Im)[nir] = std::imag((complex<double>)it.value());
			}
			else
				data[nir] = std::real(complex<double>(it.value()));
			nir++;
		}
	i=1;
	while(i<=st.rows() && jc[st.rows()-i] < 0) {
		jc[st.rows()-i] = nir;
		i++;
	}
	sparse.ir = ir;
	sparse.jc = jc;
	sparse.nir = nir;
	if(opt) {
		sparse.data = &z;
	}
	else
		sparse.data = data;
	var = Mat_VarCreate(NULL /* no-name */, MAT_C_SPARSE, MAT_T_DOUBLE, 2, dims, &sparse, opt);
//	if(var != NULL)
//		Mat_VarPrint(var,1);
	delete[] jc;
	delete[] ir;
	delete[] data;
	return var;
}

template<typename FPP>
FPP Faust::MatSparse<FPP, Cpu>::normL1(faust_unsigned_int& col_id) const
{
	faust_unsigned_int i, j, max_j;
	FPP sum, max_sum;
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	for(j=0;j<this->getNbCol();j++)
	{
		vec=mat.block(0,j,this->getNbRow(),1);
		for(i=0,sum=0;i<this->getNbRow();i++)
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
FPP Faust::MatSparse<FPP, Cpu>::normL1() const
{
	faust_unsigned_int id;
	return normL1(id);
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatSparse<FPP,Cpu>::get_col(faust_unsigned_int id) const
{
	if(id > this->getNbCol())
		handleError("Faust::MatSparse", "Too big column index passed to get_col().");
	Eigen::Matrix<FPP, Eigen::Dynamic,1> vec;
	vec = mat.col(id);
	return Vect<FPP,Cpu>(this->getNbRow(),vec.data());
}

template<typename FPP>
Faust::MatSparse<FPP, Cpu>* Faust::MatSparse<FPP, Cpu>::randMat(unsigned int num_rows, unsigned int num_cols, double density)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0, 1);
	typedef Eigen::Triplet<FPP> T;
	std::vector<T> tripletList;
	tripletList.reserve((int)(num_rows*num_cols*density));
	MatDense<FPP,Cpu>* dense = Faust::MatDense<FPP,Cpu>::randMat(num_rows, num_cols);
	for(int i=0;i<num_rows;i++)
		for(int j=0;j<num_cols;j++)
				if(distribution(generator) < density)
				{
					tripletList.push_back(T(i,j,(*dense)(i,j)));
				}
	MatSparse<FPP, Cpu>* fsmat = new MatSparse<FPP,Cpu>();
	Eigen::SparseMatrix<FPP,Eigen::RowMajor> mat(num_rows, num_cols);
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	fsmat->mat = mat;
	fsmat->update_dim();
//	fsmat->Display();
//	fsmat->mat.setFromTriplets(tripletList.begin(), tripletList.end());
	delete dense;
	return fsmat;
}


	
#endif
