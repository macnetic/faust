template<typename FPP>
void Faust::MatDiag<FPP>::Display() const
{
	std::cout<<this->to_string();
 }

template<typename FPP>
std::string Faust::MatDiag<FPP>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts /* false by default */) const
{
	std::ostringstream  str;
	str << Faust::MatGeneric<FPP,Cpu>::to_string(Diag, transpose);
	if(isZeros)
		str <<"zeros matrix flag" << std::endl;
	else if (displaying_small_mat_elts && this->dim1 < 100)
	{
		for (int i=0 ; i<this->dim1 ; i++)
			str << mat.diagonal()(i) << " " ;
		str << std::endl;
	}
	return str.str();
}

template<typename FPP>
void Faust::MatDiag<FPP>::operator*=(const FPP alpha)
{
	mat = mat*alpha;
}

template<typename FPP>
void Faust::MatDiag<FPP>::faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB) const
{
	Faust::spgemm(MatSparse<FPP,Cpu>(*this), B, C, alpha, beta, typeA, typeB);
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDiag<FPP>::multiply(const Vect<FPP,Cpu> &v) const
{
	Faust::Vect<FPP,Cpu> v_(this->getNbRow());
	v_.vec = mat * v.vec;
	return v_;
}

template<typename FPP>
void  Faust::MatDiag<FPP>::multiply(Faust::Vect<FPP,Cpu> & vec, char opThis) const
{
	if(opThis = 'H')
		vec.vec = mat.diagonal().conjugate().array() * vec.vec.array();
	else //if (opThis == 'N' || opThis == 'T')
		vec.vec = mat * vec.vec;
}

template<typename FPP>
void  Faust::MatDiag<FPP>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
{
	if(opThis = 'H')
		M.mat = mat.diagonal().conjugate().asDiagonal() * M.mat;
	else //if (opThis == 'N' || opThis == 'T')
		M.mat = mat * M.mat;
}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDiag<FPP>::Clone(const bool isOptimize) const
{
	return new Faust::MatDiag<FPP>(*this);
}

template<typename FPP>
matvar_t* Faust::MatDiag<FPP>::toMatIOVar(bool transpose, bool conjugate, const char* var_name/*=nullptr*/) const
{
	//TODO
	return nullptr;
}

template<typename FPP>
Real<FPP> Faust::MatDiag<FPP>::normL1(const bool transpose) const
{
	faust_unsigned_int i;
	return normL1(i,transpose);
}

template<typename FPP>
Real<FPP> Faust::MatDiag<FPP>::norm() const
{
	return mat.diagonal().norm();
}

template<typename FPP>
Real<FPP> Faust::MatDiag<FPP>::normL1(faust_unsigned_int& col_id, const bool transpose) const
{
	const FPP* data = getData();
	Real<FPP> max = 0, a;
	for(faust_unsigned_int i=0;i<min(this->dim1, this->dim2);i++)
		if(Faust::fabs(a=Faust::fabs(data[i])) > Faust::fabs(max))
		{
			max = a;
			col_id = i;
		}
	return max;
}

template<typename FPP>
Faust::Vect<FPP,Cpu> Faust::MatDiag<FPP>::get_col(faust_unsigned_int id) const
{
	Faust::Vect<FPP,Cpu> v(this->getNbRow());
	memset(v.getData(),0, sizeof(FPP)*this->getNbRow());
	if(id>=this->getNbCol())
		handleError("Matdiag", "column index is out of dimension size.");
	if(id < min(this->getNbRow(), this->getNbCol()))
		v.getData()[id] = getData()[id];
	return v;
}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDiag<FPP>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
{
	if(num_cols+col_id_start>=this->getNbCol())
		handleError("Matdiag", "column index is out of dimension size.");
	faust_unsigned_int dmin = min(this->getNbRow(), this->getNbCol());
	FPP * data = new FPP[num_cols];
	if(col_id_start > dmin)
		memset(data, 0, sizeof(FPP)*num_cols);
	else {
		if(col_id_start+num_cols > dmin)
		{
			faust_unsigned_int overflow = col_id_start+num_cols-dmin;
			memset(data+num_cols-overflow, 0, sizeof(FPP)*overflow);
		}
		memcpy(data, getData()+col_id_start, num_cols*sizeof(FPP));
	}
	Faust::MatDiag<FPP> * ret = new Faust::MatDiag<FPP>(this->getNbRow(), num_cols, data);
	delete[] data;
	return ret;
}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDiag<FPP>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
{
	if(num_rows+row_id_start>=this->getNbRow())
		handleError("Matdiag", "row index is out of dimension size.");
	faust_unsigned_int dmin = min(this->getNbRow(), this->getNbCol());
	FPP * data = new FPP[num_rows];
	if(row_id_start > dmin)
		memset(data, 0, sizeof(FPP)*num_rows);
	else {
		if(row_id_start+num_rows > dmin)
		{
			faust_unsigned_int overflow = row_id_start+num_rows-dmin;
			memset(data+num_rows-overflow, 0, sizeof(FPP)*overflow);
		}
		memcpy(data, getData()+row_id_start, num_rows*sizeof(FPP));
	}
	Faust::MatDiag<FPP> * ret = new Faust::MatDiag<FPP>(num_rows, this->getNbCol(), data);
	delete[] data;
	return ret;
}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDiag<FPP>::get_cols(const faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
{
	return MatSparse<FPP,Cpu>(*this).get_cols(col_ids, num_cols);
}

template<typename FPP>
Faust::MatGeneric<FPP,Cpu>* Faust::MatDiag<FPP>::get_rows(const faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
{
	return MatSparse<FPP,Cpu>(*this).get_rows(row_ids, num_rows);
}
template<typename FPP>
std::list<std::pair<int,int>> Faust::MatDiag<FPP>::nonzeros_indices(const double& tol/*=0*/) const
{
	std::list<std::pair<int,int>> nz_inds;
	for(int i=0;i<this->getNbRow();i++)
		if(std::abs(mat.diagonal()(i, i)) > tol)
			nz_inds.push_back(std::make_pair(i,i));
	return nz_inds;
}

template<typename FPP>
size_t Faust::MatDiag<FPP>::getNBytes() const
{
	return sizeof(FPP)*this->dim1*this->dim2;
}

template<typename FPP>
void Faust::MatDiag<FPP>::setZeros()
{
	memset((void*)getData(), 0, sizeof(FPP) * this->dim1*this->dim2);
	isZeros = true;
	this->is_identity = false;
}

template<typename FPP>
bool Faust::MatDiag<FPP>::containsNaN() const
{
	return Eigen::Map<Eigen::Matrix<FPP, Eigen::Dynamic, Eigen::Dynamic>>(const_cast<FPP*>(getData()) /* no worry, just a read access */, min(this->dim1, this->dim2), 1).hasNaN();
}

namespace Faust
{
	template<typename FPP>
		MatDense<FPP, Cpu> MatDiag<FPP>::to_dense() const
		{
			MatDense<FPP, Cpu> dense(this->getNbCol(), this->getNbCol());
			dense.setOnes();
			multiply(dense, 'N');
			return dense;
		}
}
