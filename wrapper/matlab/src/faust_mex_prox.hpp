#include "class_handle.hpp"
template <typename SCALAR, FDevice DEV>
void faust_prox(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	if(nlhs!=1)
		mexErrMsgTxt("prox(): 1 output variable is expected.");
	if(nrhs < 4)
		mexErrMsgTxt("prox(): wrong number of arguments (must be 4)");

	unsigned int constraint_type = (unsigned int) mxGetScalar(prhs[2]);
	bool const is_param_scalar = (1u == mxGetNumberOfElements(prhs[3]));
	double scal_param;
	const mxArray* mxMat_in = prhs[1];

	const size_t nrows_mat = mxGetM(mxMat_in);
	const size_t ncols_mat = mxGetN(mxMat_in);
	SCALAR* mat_data = NULL;
	mxArray2Ptr(mxMat_in, mat_data);
	Faust::MatDense<SCALAR, DEV> mat(nrows_mat, ncols_mat, mat_data);
	std::unique_ptr<Faust::MatDense<SCALAR, DEV>> mat_param = nullptr;
	bool normalized = true, pos = false;

	if(is_param_scalar)
	{
		scal_param = mxGetScalar(prhs[3]);
	}
	else
	{
		mxArray2Ptr(prhs[3], mat_data);
		// constraint matrix has normally same dim as constrainted matrix (TODO: must be checked from matlab code)
		// it's verified in C++
		mat_param = std::unique_ptr<Faust::MatDense<SCALAR, DEV>>(new Faust::MatDense<SCALAR, DEV>(nrows_mat, ncols_mat, mat_data));
	}

	if(nrhs > 4)
	{
		normalized = (bool) mxGetScalar(prhs[4]);
	}

	if(nrhs > 5)
	{
		pos = (bool) mxGetScalar(prhs[5]);
	}


	try
	{
		switch(constraint_type)
		{
			case CONSTRAINT_NAME_SP:
				Faust::prox_sp(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_SPCOL:
				Faust::prox_spcol(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_SPLIN:
				Faust::prox_splin(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_SPLINCOL:
				Faust::prox_splincol(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_SP_POS:
				Faust::prox_sp_pos(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_NORMLIN:
				Faust::prox_normlin(mat, scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_NORMCOL:
				Faust::prox_normcol(mat, scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_CONST:
				Faust::prox_const(mat, *mat_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_BLKDIAG:
				//not impl. yet in cpp core (but the prox is)
				break;
			case CONSTRAINT_NAME_SUPP:
				Faust::prox_supp(mat, *mat_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_CIRC:
				Faust::prox_circ(mat, normalized, pos);
				break;
			case CONSTRAINT_NAME_ANTICIRC:
				Faust::prox_anticirc(mat, normalized, pos);
				break;
			case CONSTRAINT_NAME_TOEPLITZ:
				Faust::prox_toeplitz(mat, normalized, pos);
				break;
			case CONSTRAINT_NAME_HANKEL:
				Faust::prox_hankel(mat, normalized, pos);
				break;
			case CONSTRAINT_NAME_SKPERM:
				Faust::prox_skperm(mat, (faust_unsigned_int) scal_param, normalized, pos);
				break;
			case CONSTRAINT_NAME_ID:
				Faust::prox_id(mat, normalized, pos);
				break;
            case CONSTRAINT_NAME_TRIL_SP:
				Faust::prox_tril_sp(mat, (faust_unsigned_int) scal_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_TRIU_SP:
				Faust::prox_triu_sp(mat, (faust_unsigned_int) scal_param, normalized, pos);
                break;
            case CONSTRAINT_NAME_SYMM_SP:
				Faust::prox_symm_sp(mat, (faust_unsigned_int) scal_param, normalized, pos);
                break;
			default:
				mexErrMsgTxt("Unknown constraint name/type.");
		}
		plhs[0] = FaustMat2mxArray(mat);
	}
	catch(std::domain_error& e)
	{
		mexErrMsgTxt("Can't normalize because norm is zero.");
	}
	// if(mat_param != nullptr)
	//	delete mat_param;

}

template <typename SCALAR, FDevice DEV>
void faust_prox_blockdiag(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{
	if(nlhs!=1)
		mexErrMsgTxt("prox_blockdiag(): 1 output variable is expected.");
	if(nrhs < 4)
		mexErrMsgTxt("prox_blockdiag(): wrong number of arguments (must be 4)");

	SCALAR* mat_data = NULL;
	const mxArray* mxMat_in = prhs[1];
	const mxArray* mx_m_ptr = prhs[2], *mx_n_ptr = prhs[3];
	double * m_ptr = NULL, *n_ptr = NULL;
	size_t nblocks;

	const size_t nrows_mat = mxGetM(mxMat_in);
	const size_t ncols_mat = mxGetN(mxMat_in);
	mxArray2Ptr(mxMat_in, mat_data);
	Faust::MatDense<SCALAR, DEV> mat(nrows_mat, ncols_mat, mat_data);

	mxArray2Ptr<double>(mx_n_ptr, n_ptr);
	mxArray2Ptr<double>(mx_m_ptr, m_ptr);

	std::vector<faust_unsigned_int> m_vec;
	std::vector<faust_unsigned_int> n_vec;
	bool normalized = (bool) mxGetScalar(prhs[2]);
	bool pos = (bool) mxGetScalar(prhs[3]);

	nblocks = mxGetM(mx_m_ptr);
	if(nblocks == 1 && mxGetN(mx_n_ptr) > 1)
		nblocks = mxGetN(mx_n_ptr);
	if(mxGetM(mx_m_ptr) != nblocks && mxGetN(mx_n_ptr) != nblocks)
		mexErrMsgTxt("the number of row and col offsets must agree.");

	for(int i = 0; i < nblocks; i++)
	{
		m_vec.push_back((faust_unsigned_int)(m_ptr[i]));
		n_vec.push_back((faust_unsigned_int)(n_ptr[i]));
	}
	try
	{
		Faust::prox_blockdiag(mat, m_vec, n_vec, normalized, pos);
	}
	catch(std::domain_error& e)
	{
		mexErrMsgTxt("Can't normalize because norm is zero.");
	}

	plhs[0] = FaustMat2mxArray(mat);
}
