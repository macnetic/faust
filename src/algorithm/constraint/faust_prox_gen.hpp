	template<typename FPP, FDevice DEV, typename FPP2>
		Faust::MatGeneric<FPP,DEV>* Faust::prox_normcol_gen(Faust::MatDense<FPP,DEV> & M,FPP2 s, const bool normalized/*=false*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_normcol(M, s, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);

}

	template<typename FPP, FDevice DEV, typename FPP2>
Faust::MatGeneric<FPP,DEV>* Faust::prox_normlin_gen(Faust::MatDense<FPP,DEV> & M,FPP2 s, const bool normalized/*=false*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_normlin(M, s, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_blockdiag_gen(Faust::MatDense<FPP,DEV> & M, std::vector<faust_unsigned_int>& m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized /*= false*/, const bool pos /*= false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_blockdiag(M, m_vec, n_vec, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);

}

template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_blockdiag_gen(Faust::MatDense<FPP,DEV> & M, Faust::MatDense<FPP,DEV> mn_vec, const bool normalized /* default to false */, const bool pos, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_blockdiag(M, mn_vec, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV> Faust::MatGeneric<FPP,DEV>* Faust::prox_hankel_gen(Faust::MatDense<FPP, DEV> & M, const bool normalized /*= true*/, const bool pos /*= false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_hankel(M, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV> Faust::MatGeneric<FPP,DEV>* Faust::prox_toeplitz_gen(Faust::MatDense<FPP, DEV> & M, const bool normalized /*= true*/, const bool pos /*= false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_toeplitz(M, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV> Faust::MatGeneric<FPP,DEV>* Faust::prox_circ_gen(Faust::MatDense<FPP, DEV> & M, const bool normalized /*= true*/, const bool pos /*= false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_circ(M, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV> Faust::MatGeneric<FPP,DEV>* Faust::prox_anticirc_gen(Faust::MatDense<FPP, DEV> & M, const bool normalized /*= true*/, const bool pos /*= false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_anticirc(M, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_skperm_gen(Faust::MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{

	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k*dim1, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_skperm(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_splincol_gen(Faust::MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{

	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k*dim1, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_splincol(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_spcol_gen(Faust::MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{

	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k*dim1, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_spcol(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
		Faust::MatGeneric<FPP,DEV>* Faust::prox_splin_gen(Faust::MatDense<FPP, DEV> & M, const unsigned int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k*dim2, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_splin(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_supp_gen(Faust::MatDense<FPP,DEV> & M, const Faust::MatDense<FPP,DEV> & supp, const bool normalized/*=true*/, const bool pos/*=false*/,  const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_supp(M, supp, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);

}
	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_const_gen(Faust::MatDense<FPP,DEV> & M, const Faust::MatDense<FPP,DEV> & c, const bool normalized/*=true*/, const bool pos/*=false*/,  const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_const(M, c, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_id_gen(Faust::MatDense<FPP,DEV> & M, const bool normalized/*=false*/, const bool pos/*=false*/,  const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	prox_id(M, normalized, pos);
	auto out_is_dense = Faust::sparse_size<FPP>(M.getNonZeros(), dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

	template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_sp_gen(Faust::MatDense<FPP, DEV> & M, faust_unsigned_int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_sp(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_triu_sp_gen(Faust::MatDense<FPP, DEV> & M, faust_unsigned_int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_triu_sp(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

template<typename FPP, FDevice DEV>
Faust::MatGeneric<FPP,DEV>* Faust::prox_tril_sp_gen(Faust::MatDense<FPP, DEV> & M, faust_unsigned_int k,  const bool normalized/*=true*/, const bool pos/*=false*/, const MatType forcedType/*=None*/)
{
	const faust_unsigned_int dim1 = M.getNbRow();
	const faust_unsigned_int dim2 = M.getNbCol();
	auto out_is_dense = Faust::sparse_size<FPP>(k, dim1) > Faust::dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
	prox_tril_sp(M, k, normalized, pos);
	if(out_is_dense)
		return new Faust::MatDense<FPP,DEV>(M);
	else
		return new Faust::MatSparse<FPP,DEV>(M);
}

//TODO: this function is Cpu specialized, it needs to be generalized to any FDevice or explicitly specialized for GPU2
//template<typename FPP, FDevice DEV>
//Faust::MatGeneric<FPP,DEV>* Faust::prox_sp_gen(Faust::MatDense<FPP,DEV> & M, faust_unsigned_int k, const bool normalized /* true by deft */, const bool pos, const MatType forcedType/*=None*/)
//{
//	// M is the Dense matrix on which to compute projection
//	// M is also the output matrix if k is high enough
//	// otherwise spM is used as output
//	// The output matrix is returned as a MatGeneric (as explained it could be MatSparse or MatDense)
//	const faust_unsigned_int dim1 = M.getNbRow();
//	const faust_unsigned_int dim2 = M.getNbCol();
//	const faust_unsigned_int nnz = M.getNonZeros();
//
//	auto out_is_dense = sparse_size<FPP>(k, dim1) > dense_size<FPP>(dim1, dim2) && forcedType == None || forcedType == Dense;
//
//	if(pos) Faust::pre_prox_pos(M);
//
//	const faust_unsigned_int numel = dim1*dim2;
//
//	if (k < numel /* && k < M.getNonZeros()*/)
//	{
//		const std::vector<FPP> vec(M.getData(), M.getData()+numel);
//		std::vector<int> index;
//		Faust::sort_idx(vec, index, k);
//		index.erase(index.begin()+k, index.end());
//
//		if(out_is_dense)
//		{
//			M.setZeros();
//			for (int i=0 ; i<index.size() ; i++)
//				M.getData()[index[i]] = vec[index[i]];
//
//			if(normalized) M.normalize();
//			return new Faust::MatDense<FPP,DEV>(M);
//		}
//		else
//		{
//			FPP *values = new FPP[index.size()];
//			unsigned int *row_ids = new unsigned int[index.size()];
//			unsigned int *col_ids = new unsigned int[index.size()];
//
//
//
//			std::vector<Eigen::Triplet<FPP> > tripletList;
//			for (int i=0 ; i<index.size() ; i++)
//			{
//				values[i] = vec[index[i]];
//				row_ids[i] = index[i]%dim1;
//				col_ids[i] = index[i]/dim1;
//			}
//
//			if(normalized)
//			{
//				Faust::Vect<FPP,DEV> nvalues(index.size(), values);
//				nvalues.normalize();
//				delete [] values;
//				values = nvalues.getData();
//			}
//			auto spM = new MatSparse<FPP,DEV>(row_ids, col_ids, values, dim1, dim2, index.size());
//			delete[] row_ids;
//			delete[] col_ids;
//			if(! normalized)
//				delete[] values;
//			return spM;
//		}
//	}
//	else
//	{
//		if(sparse_size<FPP>(M.getNonZeros(), dim1) > dense_size<FPP>(dim1, dim2))
//			return new Faust::MatDense<FPP,DEV>(M);
//		else
//		{
//			auto spM = new MatSparse<FPP,DEV>(M);
//			return spM;
//		}
//	}
//}


