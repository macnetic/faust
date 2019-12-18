template<typename FPP>
void Faust::MatDiag<FPP>::Display() const
{
	std::cout<<this->to_string();
 }

template<typename FPP>
std::string Faust::MatDiag<FPP>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts /* false by default */) const
{
	//using ostringstream because it's easier for concatenation (of complex and any number)
	std::ostringstream  str;

	str << " (" << ((typeid(mat.diagonal()(0)) == typeid(complex<double>) || typeid(mat.diagonal()(0)) == typeid(complex<float>))?"complex":"real") << ")";
	str<<" DENSE,";
	str << Faust::MatGeneric<FPP,Cpu>::to_string(transpose);
	if(isZeros)
		str <<"zeros matrix flag" << endl;
	else if (displaying_small_mat_elts && this->dim1 < 100)
	{
		for (int i=0 ; i<this->dim1 ; i++)
			str << mat.diagonal()(i) << " " ;
		str << endl;
	}
	return str.str();
}

template<typename FPP>
void MatDiag<FPP>::operator*=(const FPP alpha)
{
	mat = mat*alpha;
}

template<typename FPP>
void MatDiag<FPP>::faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB) const
{
	Faust::spgemm(MatSparse<FPP,Cpu>(*this), B, C, alpha, beta, typeA, typeB);
}

template<typename FPP>
Vect<FPP,Cpu> MatDiag<FPP>::multiply(const Vect<FPP,Cpu> &v) const
{
	Vect<FPP,Cpu> v_(this->getNbRow());
	v_.vec = mat * v.vec;
	return v_;
}

template<typename FPP>
void  MatDiag<FPP>::multiply(Vect<FPP,Cpu> & vec, char opThis) const
{
	if(opThis == 'T')
	{
		MatDiag<FPP> tmat(this->getNbCol(), this->getNbRow(), this->getData());
		vec.vec = tmat.mat * vec.vec;
	}
	else vec.vec = mat * vec.vec;
}

template<typename FPP>
void  MatDiag<FPP>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
{
	if (opThis == 'N')
		M.mat = this->mat * M.mat;
	else {
		MatDiag<FPP> tmat(this->getNbCol(), this->getNbRow(), this->getData());
		M.mat = tmat.mat * M.mat;
	}
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::Clone(const bool isOptimize) const
{
	return new MatDiag<FPP>(*this);
}

template<typename FPP>
matvar_t* MatDiag<FPP>::toMatIOVar(bool transpose, bool conjugate) const
{
	//TODO
	return nullptr;
}

template<typename FPP>
Real<FPP> MatDiag<FPP>::normL1(const bool transpose) const
{
	faust_unsigned_int i;
	return normL1(i,transpose);
}

template<typename FPP>
Real<FPP> MatDiag<FPP>::norm() const
{
	return mat.diagonal().norm();
}

template<typename FPP>
Real<FPP> MatDiag<FPP>::normL1(faust_unsigned_int& col_id, const bool transpose) const
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
Vect<FPP,Cpu> MatDiag<FPP>::get_col(faust_unsigned_int id) const
{
	Vect<FPP,Cpu> v(this->getNbRow());
	memset(v.getData(),0, sizeof(FPP)*this->getNbRow());
	if(id>=this->getNbCol())
		handleError("Matdiag", "column index is out of dimension size.");
	if(id < min(this->getNbRow(), this->getNbCol()))
		v.getData()[id] = getData()[id];
	return v;
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
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
	MatDiag<FPP> * ret = new MatDiag<FPP>(this->getNbRow(), num_cols, data);
	delete[] data;
	return ret;
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
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
	MatDiag<FPP> * ret = new MatDiag<FPP>(num_rows, this->getNbCol(), data);
	delete[] data;
	return ret;
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
{
	return MatSparse<FPP,Cpu>(*this).get_cols(col_ids, num_cols);
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
{
	return MatSparse<FPP,Cpu>(*this).get_rows(row_ids, num_rows);
}
template<typename FPP>
list<pair<int,int>> MatDiag<FPP>::nonzeros_indices() const
{
	list<pair<int,int>> nz_inds;
	for(int i=0;i<this->getNbRow();i++)
		nz_inds.push_back(make_pair(i,i));
	return nz_inds;
}

