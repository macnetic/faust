template<typename FPP>
void MatDiag<FPP>::operator*=(const FPP alpha)
{

}

template<typename FPP>
void MatDiag<FPP>::faust_gemm(const Faust::MatDense<FPP,Cpu> & B, Faust::MatDense<FPP,Cpu> & C,const FPP & alpha, const FPP & beta, char typeA, char typeB)const
{
}

template<typename FPP>
Vect<FPP,Cpu> MatDiag<FPP>::multiply(const Vect<FPP,Cpu> &v) const
{
	Vect<FPP,Cpu> v_(v);
	return v;
}

template<typename FPP>
void  MatDiag<FPP>::multiply(Vect<FPP,Cpu> & vec, char opThis) const
{
	//TODO
}

template<typename FPP>
void  MatDiag<FPP>::multiply(MatDense<FPP,Cpu> & M, char opThis) const
{
	//TODO
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

