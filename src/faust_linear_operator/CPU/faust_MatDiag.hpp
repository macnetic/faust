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
	//TODO
	return nullptr;
}

template<typename FPP>
matvar_t* MatDiag<FPP>::toMatIOVar(bool transpose, bool conjugate) const
{
	//TODO
	return nullptr;
}

template<typename FPP>
FPP MatDiag<FPP>::normL1(const bool transpose) const
{
	//TODO
	return FPP(0);
}

template<typename FPP>
FPP MatDiag<FPP>::norm() const
{
	//TODO
	return FPP(0);
}

template<typename FPP>
FPP MatDiag<FPP>::normL1(faust_unsigned_int& col_id, const bool transpose) const
{
	//TODO
	return FPP(0);
}

template<typename FPP>
Vect<FPP,Cpu> MatDiag<FPP>::get_col(faust_unsigned_int id) const
{
	//TODO
	Vect<FPP,Cpu> v;
	return v;
}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_cols(faust_unsigned_int col_id_start, faust_unsigned_int num_cols) const
{
	//TODO
	return nullptr;

}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_rows(faust_unsigned_int row_id_start, faust_unsigned_int num_rows) const
{
	//TODO
	return nullptr;

}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_cols(faust_unsigned_int* col_ids, faust_unsigned_int num_cols) const
{
	//TODO
	return nullptr;

}

template<typename FPP>
MatGeneric<FPP,Cpu>* MatDiag<FPP>::get_rows(faust_unsigned_int* row_ids, faust_unsigned_int num_rows) const
{
	//TODO
	return nullptr;
}

