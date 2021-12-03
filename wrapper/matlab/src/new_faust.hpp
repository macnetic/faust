template <typename SCALAR, FDevice DEV>
void new_faust(const mxArray **prhs, const int nrhs, mxArray **plhs, const int nlhs)
{

	if(nlhs!=1)
		mexErrMsgTxt("1 output is expected.");
	if((nrhs<2) || (nrhs>4))
		mexErrMsgTxt("between 1 and 3 inputs are expected.");

	if(!mxIsCell(prhs[1]))
		mexErrMsgTxt("1st arg input must be a cell-array");

	std::vector<Faust::MatGeneric<SCALAR, DEV>*> list_factor;
	mwSize nb_element = mxGetNumberOfElements(prhs[1]);

	Faust::TransformHelper<SCALAR, DEV> tmp_bsr_th;

	if (nb_element != 0)
	{
		mxArray * mxMat, *mx_tmp;
		for (int i=0; i < nb_element; i++)
		{
			mxMat=mxGetCell(prhs[1],i);
			if(mxIsCell(mxMat))
			{
				// adding at MatBSR matrix
				// the BSR cell is assumed well-formed (checked in the matfaust matlab code)
				mx_tmp = mxGetCell(mxMat, 1);
				auto nrows = (int) mxGetScalar(mx_tmp);
				mx_tmp = mxGetCell(mxMat, 2);
				auto ncols = (int) mxGetScalar(mx_tmp);
				mx_tmp = mxGetCell(mxMat, 3);
				auto bnnz = (int) mxGetScalar(mx_tmp);
				SCALAR* bdata = nullptr;
				mx_tmp = mxGetCell(mxMat, 4);
				auto bnrows = mxGetM(mx_tmp); // block number of rows
				auto bncols = mxGetN(mx_tmp) / bnnz; // block number of cols
				mxArray2Ptr(mx_tmp, bdata);
				mx_tmp = mxGetCell(mxMat, 5);
				int* bcolinds_ = nullptr;
				mxArray2Ptr(mx_tmp, bcolinds_);
				auto bcolinds = new int[bnnz];
				for(int i=0;i<bnnz;i++)
					bcolinds[i] = bcolinds_[i]-1;
				mx_tmp = mxGetCell(mxMat, 6);
				int* browptr_ = nullptr;
				mxArray2Ptr(mx_tmp, browptr_);
				int* browptr = new int[nrows/bnrows+1];
				browptr[0] = 0;
				for(int i=1;i<bnrows+1;i++)
					browptr[i] = browptr_[i-1] + browptr[i-1];
//				auto bsr_mat = new Faust::MatBSR<SCALAR, DEV>(nrows, ncols, bnrows, bncols, bnnz, bdata, browptr, bcolinds);
//				list_factor.push_back(bsr_mat);
				// hack to use the same code for GPU2 which is not MatBSR compatible yet
				tmp_bsr_th.push_back(bdata, browptr, bcolinds, nrows, ncols, bnnz, bnrows, bncols, /*optimizedCopy*/ false, false, false);
				auto gen_fac = tmp_bsr_th.get_gen_fact_nonconst(tmp_bsr_th.size()-1);
				list_factor.push_back(gen_fac);
				delete[] bdata;
				delete[] browptr_;
				delete[] browptr;
				delete[] bcolinds;
				delete[] bcolinds_;
				continue;
			}
			// add a MatDense or MatSparse factor
			//TODO: test if the i-th element of the cell is a matrix or a cell (BSR encoding)
			//TODO: if the element is a valid BSR encoding create the BSRMatrix here and add it to list_factor
			concatMatGeneric<SCALAR, DEV>(mxMat, list_factor);

		}
	}

	// 2nd input (optional) : multiplicative scalar
	SCALAR lambda = 1.0;
	if (nrhs > 2)
		lambda = (SCALAR) mxGetScalar(prhs[2]);


	// 3rd input (optional) : boolean to determine the type of copy
	bool optimizedCopy = false;
	if (nrhs > 3)
	{
		double optimizedCopy_inter = mxGetScalar(prhs[3]);
		if ((optimizedCopy_inter != 1.0) && (optimizedCopy_inter != 0.0))
			mexErrMsgTxt("invalid boolean argument.");

		optimizedCopy = (bool) optimizedCopy_inter;

	}

	Faust::TransformHelper<SCALAR,DEV>* F = new Faust::TransformHelper<SCALAR,DEV>(list_factor, lambda, optimizedCopy, /* cloning */ false);

	plhs[0]=convertPtr2Mat<Faust::TransformHelper<SCALAR,DEV> >(F);

	return;
}

