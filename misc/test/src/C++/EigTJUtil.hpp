
#define OK_ASSERT(msg) \
	cout << "VALIDATED ASSERTION: " << msg <<endl;

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_err_against_Laplacian(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file)
{
	Faust::MatDense<SCALAR,DEVICE> Lap, tmp;
	Faust::MatDense<SCALAR,DEVICE> ordered_Uhat;
	SCALAR2 err0;
	SCALAR err_lap, err_lap2;
	SCALAR2 ERR_OK_TRESHOLD = 1e-2;
	init_faust_mat_from_matio(Lap,conf_file,"Lap");
	err0 = Faust::init_double_from_matio(conf_file, "err0");

	Faust::MatDense<SCALAR,DEVICE> ordered_D(algo->get_Dspm(true));

	cout << "Lap fro norm: " << Lap.norm() << endl;

	// Compute squared relative error according to the Laplacian
	ordered_Uhat = algo->compute_fourier(true);
	cout << "ordered Uhat fro norm: " << ordered_Uhat.norm() << endl;
	tmp = ordered_Uhat;
	tmp.multiplyRight(ordered_D);
	ordered_Uhat.transpose();
	ordered_Uhat.conjugate();
	tmp.multiplyRight(ordered_Uhat);
	ordered_Uhat.transpose();
	ordered_Uhat.conjugate();
	tmp -= Lap;
	cout << "Error relatively to the Laplacian: " << endl;
	err_lap = SCALAR(tmp.norm()/Lap.norm());
	cout << err_lap << endl;
	cout << "matlab error =" << err0 << endl;
	assert(Faust::fabs(err_lap-err0) < ERR_OK_TRESHOLD);
	OK_ASSERT("found same value as Matlab for relative error of full(Uhat)*Dhat*full(Uhat)' against the Laplacian (precision of "<< ERR_OK_TRESHOLD << ")");
	

	//do the same but simpler with Transform object
	Faust::Transform<SCALAR,DEVICE> t = std::move(algo->get_transform(true));
	cout << "ordered Uhat norm: " << t.normFro() << endl;
	Faust::MatDense<SCALAR,DEVICE> P = t.multiply(ordered_D);//t.get_product();

//	P.multiplyRight(t.get_product('T', true /* isConj */));
	P.multiplyRight(t.get_product('H'));
	P -= Lap;
	cout << "Error relatively to the Laplacian: (computed from Transform object)" << endl;
	err_lap2 = SCALAR(P.norm()/Lap.norm());
	cout << err_lap2 << endl;

	assert(Faust::fabs(err_lap2-err0) < ERR_OK_TRESHOLD);
	OK_ASSERT("found same value as Matlab for relative error of Uhat*Dhat*Uhat' against the Laplacian (precision of "<< ERR_OK_TRESHOLD << ", Uhat is a Faust::Transform)");
	assert(Faust::fabs(err_lap2-err_lap) < 1e-5);
}

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_eigentransform(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file)
{
	Faust::MatDense<SCALAR,DEVICE> U, tmp;
	SCALAR2 err1;
	double ERR_OK_TRESHOLD = 1e-3;

	init_faust_mat_from_matio(U,conf_file,"U");
	err1 = Faust::init_double_from_matio(conf_file, "err1");

	// Compute manually the dense matrix of the eigenvectors (from the Givens transform / factors)
	Faust::MatDense<SCALAR,DEVICE> Uhat(U.getNbRow(), U.getNbCol());
	Uhat.setEyes();
	const vector<Faust::MatSparse<SCALAR,DEVICE>>& givens_facts = algo->get_facts();
	for(int i=givens_facts.size()-1;i>=0;i--)
		givens_facts[i].multiply(Uhat, 'N');

	cout << "Uhat fro norm: " << Uhat.norm() << endl;

	// Compute again the dense matrix but this time using the provided internal function
	Faust::MatDense<SCALAR,DEVICE> Uhat2 = algo->compute_fourier();
	// ordered according to corresponding eigenvalues in ascending order
	Faust::MatDense<SCALAR,DEVICE> ordered_Uhat2 = algo->compute_fourier(true);

	// construct the ordered matrix manually to verify it works
	Faust::MatDense<SCALAR,DEVICE> * ordered_Uhat = Uhat.get_cols(algo->get_ord_indices());

	Uhat2 -= Uhat;
	ordered_Uhat2 -= *ordered_Uhat;
	cout << "norm(Uhat2-Uhat): " << Uhat2.norm() << endl;
	cout << "norm(ordered_Uhat2-ordered_Uhat): " << ordered_Uhat2.norm() << endl;

	//verify Uhat is properly computed (when ordered or not)
	assert(Uhat2.norm()==SCALAR(0));
	assert(ordered_Uhat2.norm()==SCALAR(0));
	OK_ASSERT("Uhat and non-ordered Uhat have the same Fro. norm.");

	ordered_Uhat2 = algo->compute_fourier(true);
	tmp = ordered_Uhat2;
	tmp -= U;
	cout << "Error relatively to Fourier matrix (computed with matlab eig()):" << endl;
	cout << tmp.norm()/U.norm() << endl;
	cout << "matlab err: " << err1 << endl;
	assert(Faust::fabs(tmp.norm()/U.norm()-err1)/err1 < ERR_OK_TRESHOLD*60);
	OK_ASSERT("Relative error ||Uhat-U||/||U|| is the same as Matlab's (relative precision: " << ERR_OK_TRESHOLD*60 << ").");
}

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_pivot_choices(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file, const float same_pivot_target_rate /* default to .99 */, const int stop_count_ite /*= 57*/)
{
	vector<pair<int,int>> coord_choices = algo->get_coord_choices();
	int J = Faust::init_int_from_matio(conf_file, "J");

	mat_t *matfp;

	matfp = Mat_Open(conf_file,MAT_ACC_RDONLY);
	if ( NULL == matfp ) {
		fprintf(stderr,"Error opening MAT file %s", conf_file);
		exit(EXIT_FAILURE);
	}
	matvar_t* cell_choices = Mat_VarRead(matfp, "choices");
	int* p_choices = new int[J];
	int* q_choices = new int[J];
	cout << "=======================================" << endl;
	if ( NULL != cell_choices) {
		for (int i = 0; i < J*2; i+=2 )
		{
			p_choices[i/2] = (int)((double*) (cell_choices->data))[i];
			q_choices[i/2] = (int)((double*) (cell_choices->data))[i+1];
		}
		free(cell_choices);
	}
	Mat_Close(matfp);

	bool eq; // are ref and test pivots equal ?
	bool all_eq = true; //all past pivots are eq

	int count_eq_pivots = 0; //before iteration stop_count_ite

	for(int j=0;j<J;j++)
	{
		cout << "ite=" << j;
		cout << " (1-base index) ref. p=" << p_choices[j] << ", q=" << q_choices[j] << ", algo-> p=" << coord_choices[j].first+1 << " q=" << coord_choices[j].second+1;
		eq = p_choices[j] == coord_choices[j].first+1 && q_choices[j] == coord_choices[j].second+1;
		all_eq = all_eq && eq;
		cout << " equal=" << eq;
		cout << " all ites eq=" << all_eq << endl;
		if(j < stop_count_ite && eq) count_eq_pivots++;
	}
	cout << "rate of equal pivot choices (C++/matlab) from iteration #0 to #"<< stop_count_ite << "/" << J << ": " << count_eq_pivots/(float)stop_count_ite << endl;
	assert(count_eq_pivots/(float)stop_count_ite > same_pivot_target_rate);
	OK_ASSERT("More than "<< same_pivot_target_rate*100 <<"% of pivots selected in the iterations (0 to " << stop_count_ite << ") are the same as the Matlab ones");
	delete []p_choices;
	delete []q_choices;
}

template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_eigenvalues(Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file)
{
	SCALAR2 ERR_OK_TRESHOLD = 1e-2;
	Faust::MatDense<SCALAR,DEVICE> tmp;
	// Get computed eigenvalues without order
	const Faust::MatSparse<SCALAR,DEVICE> D = algo->get_Dspm();
	Faust::MatDense<SCALAR,DEVICE> full_D = Faust::MatDense<SCALAR,DEVICE>(D);
	cout << "D fro norm:" << D.norm() << endl;

	// Get sparse matrix of ordered eigenvalues
	Faust::MatDense<SCALAR,DEVICE> ordered_D(algo->get_Dspm(true));

	cout << "orderded D fro norm: " << ordered_D.norm() << endl;
	cout << "ordered_D:" << endl;
	ordered_D.Display();

	cout << Faust::fabs(D.norm() - ordered_D.norm()) << endl;
	assert(Faust::fabs(D.norm() - ordered_D.norm()) < 1e-4);
	OK_ASSERT("The ordered eigenvalue vector Dhat has the same euclidean norm as the non-ordered one.");

	cout << "ordered D eles" << endl;
	for(int i=0;i<ordered_D.getNbRow();i++)
	{
		cout << ordered_D.getData()[i*ordered_D.getNbRow()+i] << " " << full_D(i,i) << endl;
		assert(i >= ordered_D.getNbRow()-1 || Faust::fabs(ordered_D(i,i)) <= Faust::fabs(ordered_D(i+1,i+1)));
	}
	OK_ASSERT("abs(Dhat(i)) <= abs(Dhat(i+1)) for i from 0 to n-1 (n the the Laplacian order).");

	// Measure error between matlab Dhat and C++ Dhat
	Faust::MatDense<SCALAR,DEVICE> ref_Dhat;
	init_faust_mat_from_matio(ref_Dhat,conf_file,"Dhat");
	tmp = ref_Dhat;
	tmp -= ordered_D;
	cout << "Relative difference on Dhat between matlab ref script and C++ impl:" << endl;
	SCALAR test_D_err = tmp.norm()/ref_Dhat.norm();
	cout << test_D_err << endl;
	SCALAR2 ref_D_err = Faust::init_double_from_matio(conf_file, "err0"); // not used because it's recomputed
	assert(Faust::fabs(test_D_err) < ERR_OK_TRESHOLD);
	OK_ASSERT("Relative error norm(Dhat-Dhat_matlab)/norm(Dhat) is lower or equal to "<< ERR_OK_TRESHOLD << ".");

	cout << "ordered D eles" << endl;
	for(int i=0;i<ordered_D.getNbRow();i++)
	{
		cout << "ref_D(i,i) = "<< ref_Dhat(i,i) << " testD(i,i) =" << ordered_D(i,i) << endl;

		cout << Faust::fabs(ref_Dhat(i,i)-ordered_D(i,i)) << endl;
		assert(Faust::fabs(ref_Dhat(i,i)-ordered_D(i,i)) < ERR_OK_TRESHOLD*10);
	}
	OK_ASSERT("For each i abs(Dhat(i)-Dhat_matlab(i)) <= "<< ERR_OK_TRESHOLD*10 << ".");
}


template<typename SCALAR, FDevice DEVICE, typename SCALAR2>
void test_ite_errors(const Faust::@GIVENS_CLASS@ /*EigTJ<SCALAR,DEVICE,SCALAR2> or EigTJComplex<SCALAR,DEVICE,SCALAR2> */ * algo, const char* conf_file, int ref_ite_period /* default to 0 to use algo ERROR_CALC_PERIOD */)
{
	//
	assert(algo->ERROR_CALC_PERIOD == 100); /** C++ impl. must calculate the error at the same rate than matlab */
	Faust::MatDense<double,Cpu> ref_errs;
	init_faust_mat_from_matio(ref_errs,conf_file,"err");
	const vector<SCALAR2>& errs = algo->get_errs();
	SCALAR2 d;
	int i = 0;
	if(ref_ite_period == 0)
	{
		ref_ite_period = algo->ERROR_CALC_PERIOD; //matlab period of error storage
		i = ref_ite_period-1;
	}
	SCALAR2 ERR_OK_TRESHOLD = 1e-2;
	SCALAR2 ref_err, abs_diff;
	for(SCALAR2 test_err: errs)
	{
		ref_err = ref_errs.getData()[i];
		abs_diff = std::abs(test_err - ref_err);
		cout << "ref test_err:" << ref_err << " test test_err:" << test_err << " abs. diff:" << abs_diff << endl;
		assert(abs_diff < ERR_OK_TRESHOLD);
		i += ref_ite_period; // matlab Givens code stores errors at a period which is not necessarily the same as C++
	}
	OK_ASSERT("Errors of diagonalization measured every " << ref_ite_period << " iteration(s) are equal to Matlab's (abs. precision " << ERR_OK_TRESHOLD <<", error: ||diag(L)-offdiag(L)||^2_fro/||Lap||_fro^2)");
}

