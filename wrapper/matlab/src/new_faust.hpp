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

	if (nb_element != 0)
	{
		mxArray * mxMat;
		for (int i=0; i < nb_element; i++)
		{
			mxMat=mxGetCell(prhs[1],i);

			concatMatGeneric<SCALAR, DEV>(mxMat,list_factor);

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

