

namespace Faust {

	//TODO: move these functions to Faust::utils ?
	template<typename FPP>
		void conjugate(complex<FPP>* elts, faust_unsigned_int n)
		{
			for(faust_unsigned_int i=0; i< n; i++)
				elts[i] = complex<FPP>(elts[i].real(), - elts[i].imag());
		}

	template<typename FPP>
		void conjugate(FPP* elts, faust_unsigned_int n)
		{
		//nothing to do for real numbers
		}

	template<typename FPP>
		void fft_factors(unsigned int n, vector<MatGeneric<FPP,Cpu>*>&  v)
		{
			handleError("Faust::TransformHelper<FPP,Cpu>", "Can't get fft factors for real matrices.");
		}

	template<typename FPP>
		void fft_factors(unsigned int n, vector<MatGeneric<complex<FPP>,Cpu>*>&  v)
		{
			//Cooley-Tukey
			// number of factors to set : n+1
			v.resize(n+1);
			//TODO: rename index and new_index
			unsigned int dim_size = 1u << n, L, r, L_over_2, L_times_2;
			unsigned int* index = new unsigned int[dim_size];
			unsigned int* new_index = new unsigned int[dim_size];
			complex<FPP>* fac_data = new complex<FPP>[dim_size*2];
			int* fac_row_ptr = new int[dim_size+1]; //int and not uint because of MatSparse
			// we should change it later when MatSparse will be updated
			int* fac_col_ind = new int[dim_size*2];
			complex<FPP> wL;
			const FPP pi = std::acos(-1);
			const std::complex<double> cplx_i(0, 1); //TODO: rename to i (np for loop vars because they have their own scope)
			for(unsigned int i = 0; i < dim_size; i++)
				index[i] = i;
			memcpy(new_index, index, sizeof(unsigned int)*dim_size);
			bit_rev_permu(n, new_index, false);
			vector<complex<FPP>> ones(dim_size);
			for(typename vector<complex<FPP>>::iterator it=ones.begin(); it != ones.end(); it++)
				*it = complex<FPP>(1.0);
			MatSparse<complex<FPP>,Cpu> *P = new MatSparse<complex<FPP>,Cpu>(index, new_index, ones, dim_size, dim_size);
			MatSparse<complex<FPP>,Cpu> * factor;
//			cout << "fft_factors() P:" << endl;
//			P->Display();
			v[n] = P;
			for(unsigned int q=0; q < n; q++) //n+1 factors if counting P
			{
				L_over_2 = 1u << q;
				L = L_over_2 << 1;
				L_times_2 = L << 1;
				r = dim_size >> (q+1); // dim_size/L
				wL = std::exp(FPP(-2.0)*cplx_i*pi/FPP(L));
				//TODO: infer j from i and delete j
				//TODO: there is a potential for OpenMP opt. here and in the next loops
				for(unsigned int i=0, j=0; i < L_over_2; i++,j+=2)
				{
					fac_data[j] = 1; //identity diag
					fac_data[j+1] = std::pow(wL, index[i]);
					fac_data[j+L] = 1;
					fac_data[j+L+1] = -fac_data[j+1];
//					cout << "fac_data[j+1]: " << fac_data[j+1] << endl;
				}
				//Rhs == fac_data[0], ..., fac_data[2*L-1]
				for(unsigned int i=1;i<r;i++)
					//TODO: var L_times_2
					memcpy(fac_data+i*(L<<1), fac_data, sizeof(complex<FPP>)*(L<<1));
				//ok for the factor data, now the indices (CSR format)
//				cout << "fact_data[i]:" << endl;
//				for(int i = 0; i < dim_size*2; i++)
//					cout << fac_data[i] << " ";
//				cout << endl;
				//row_ptr
				// Rhs (see CooleyTukeyFact.m)
				for(int i=0,j=0;i<L+1;i+=2,j++)
					fac_row_ptr[j] = i;
				//the first L_over_2+1 eles of fac_row_ptr are set
				// set now the second part of Rhs
				// (ignoring the first 0 when copying)
				memcpy(fac_row_ptr+L_over_2+1, fac_row_ptr+1, sizeof(int)*L_over_2);
				//shift by L the second part
				for(int i=L_over_2+1;i<L+1;i++)
					fac_row_ptr[i] += L;
//				for(int i=0;i<L+1;i++)
//					cout << "fac_row_ptr[i]=" << fac_row_ptr[i] << " ";
//				cout << endl;
				// ok we have the first L+1 elements of fac_row_ptr
				// we need r*L+1 elements for Aq (whole factor)
				for(int i=1;i<r;i++)
				{
					//ignoring again the first ele. == 0
					memcpy(fac_row_ptr+i*L+1, fac_row_ptr+1, sizeof(int)*L);
					for(int j=i*L+1;j<(i+1)*L+1;j++)
						fac_row_ptr[j] += fac_row_ptr[i*L];
				}
//				for(int i=0;i<dim_size+1;i++)
//					cout << "fac_row_ptr[i]=" << fac_row_ptr[i] << " ";
//				cout << endl;				//fac_row_ptr is fully ready to set the factor
				// now set the column indices
				// starting from Rhs indices
				for(int i=0;i<L_over_2;i++)
				{
					fac_col_ind[(i<<1)] = i;
					fac_col_ind[(i<<1)+1] = L_over_2+i;
				}
				//copy the first half height of the matrix columns into the second buf part
				// because they are the same
				memcpy(fac_col_ind+L, fac_col_ind, sizeof(int)*L);
				// we need r times the Rhs columns
				for(int i=1;i<r;i++)
				{
					memcpy(fac_col_ind+i*L_times_2,fac_col_ind,sizeof(int)*L_times_2);
					// shift because the Rhs are diagonal blocks in factor (Aq)
					for(int j=i*L_times_2;j<(i+1)*L_times_2;j++)
						fac_col_ind[j] += i*L;
				}
//				for(int i=0;i<dim_size*2;i++)
//					cout << "fac_col_ind[i]=" << fac_col_ind[i] << " ";
//				cout << endl;
				// init the MatSparse with its buffers
				factor = new MatSparse<complex<FPP>,Cpu>(dim_size*2, dim_size, dim_size, fac_data,
						fac_row_ptr, fac_col_ind);
//				cout << "fft_factors() factor:" <<endl;
//				MatDense<complex<FPP>,Cpu>(*factor).Display();
//				cout << "id: " << n-q-1 << endl;
//				cout << factor << endl;
//				factor->Display();
				v[n-q-1] = factor; //TODO: directly affects v with the MatSparse, delete factor

			}
			delete[] index;
			delete[] new_index;
			delete[] fac_data;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
				const FPP lambda_, const bool optimizedCopy, const bool cloning_fact) : is_transposed(false),
																						is_conjugate(false),
																						is_sliced(false)
	{
		if(lambda_ != FPP(1.0))
			cerr << "WARNING: the constructor argument for multiplying the Faust by a scalar is DEPRECATED and might not be supported in next versions of FAµST." << endl;
		this->transform = make_shared<Transform<FPP,Cpu>>(facts, lambda_, optimizedCopy, cloning_fact);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : is_transposed(false), is_conjugate(false), is_sliced(false)
	{
		this->transform = make_shared<Transform<FPP,Cpu>>();
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(Transform<FPP,Cpu> &t) : is_transposed(false), is_conjugate(false), is_sliced(false)
	{
		this->transform = make_shared<Transform<FPP,Cpu>>(t);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th_left, TransformHelper<FPP,Cpu>* th_right)
		: is_transposed(false), is_conjugate(false), is_sliced(false)
		{
			this->transform = make_shared<Transform<FPP,Cpu>>(th_left->eval_sliced_Transform(), th_left->is_transposed, th_left->is_conjugate,
					th_right->eval_sliced_Transform(), th_right->is_transposed, th_right->is_conjugate);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, bool transpose, bool conjugate)
		{
			this->transform = th->transform;
			this->is_transposed = transpose?!th->is_transposed:th->is_transposed;
			this->is_sliced = th->is_sliced;
			if(th->is_sliced)
				copy_slices(th);
			this->is_conjugate = conjugate?!th->is_conjugate:th->is_conjugate;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th)
		{
			this->transform = th->transform;
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
			this->is_sliced = th->is_sliced;
			if(th->is_sliced)
				copy_slices(th);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(TransformHelper<FPP,Cpu>* th, Slice s[2])
		{
			this->transform = th->transform;
			this->is_transposed = th->is_transposed;
			this->is_conjugate = th->is_conjugate;
			if(! (s[0].belong_to(0, th->getNbRow()) || s[1].belong_to(0, th->getNbCol())))
				handleError("Faust::TransformHelper::TransformHelper(TransformHelper,Slice)", "Slice overflows a Faust dimension.");
			this->slices[0] = s[0];
			this->slices[1] = s[1];
			this->is_sliced = true;
			this->transform = make_shared<Transform<FPP,Cpu>>(*eval_sliced_Transform());
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A) const
		{
			MatDense<FPP,Cpu> M = this->transform->multiply(A, isTransposed2char());
			if(is_conjugate) M.conjugate();
			return M;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x) const
		{
			Vect<FPP,Cpu> v = this->transform->multiply(x, isTransposed2char());
			if(is_conjugate) v.conjugate();
			return v;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x, const bool transpose)
		{
			is_transposed ^= transpose;
			Vect<FPP,Cpu> v = this->transform->multiply(x, isTransposed2char());
			is_transposed ^= transpose;
			return v;
		}


	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A, const bool transpose)
		{
			is_transposed ^= transpose;
			MatDense<FPP,Cpu> M = this->transform->multiply(A, isTransposed2char());
			is_transposed ^= transpose;
			return M;
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(TransformHelper<FPP, Cpu>* th_right)
		{
			return new TransformHelper<FPP,Cpu>(this, th_right);
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(FPP& scalar)
		{
			const vector<MatGeneric<FPP,Cpu>*>& vec = transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
			//the point here is to minimize the number of copies (with direct access)
			// the constructor then will copy the factors from the vector
			Transform<FPP,Cpu>* t = new Transform<FPP,Cpu>(vec, scalar, false, true); //optimizedCopy == false, cloning_fact == true
			TransformHelper<FPP,Cpu>* th  = new TransformHelper<FPP,Cpu>(*t);
			//TODO: refactor the attributes copy ? 
			th->is_transposed = is_transposed;
			th->is_conjugate = is_conjugate;
			th->is_sliced = is_sliced;
			if(is_sliced)
			{
				th->slices[0] = Slice(slices[0]);
				th->slices[1] = Slice(slices[1]);
			}
			return th;
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::push_back(const MatGeneric<FPP,Cpu>* M)
		{
			//warning: should not be called after initialization of factors (to respect the immutable property)
			//this function is here only for python wrapper (TODO: see how to modify that wrapper in order to delete this function after)
			this->transform->push_back(M, true, is_conjugate); //2nd argument is for opt. (transforming dense matrix in sparse if possible)
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbRow() const
		{
			if(is_sliced){
				faust_unsigned_int id = is_transposed?1:0;
				return	this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return is_transposed?transform->getNbCol():transform->getNbRow();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::getNbCol() const
		{
			if(is_sliced)
			{
				faust_unsigned_int id = is_transposed?0:1;
				return this->slices[id].end_id-this->slices[id].start_id;
			}
			else
				return is_transposed?transform->getNbRow():transform->getNbCol();
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isReal() const
		{
			return this->transform->isReal();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_total_nnz() const
		{
			return this->transform->get_total_nnz();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::display() const
		{
			this->transform->Display(is_transposed);
		}

	template<typename FPP>
		string TransformHelper<FPP,Cpu>::to_string() const
		{
			return this->transform->to_string(is_transposed);
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::get_fact_nnz(const faust_unsigned_int id) const
		{
			if(id == 0 || id == this->size()-1)
				return this->transform->get_fact_nnz(is_transposed?size()-id-1:id);
			else
				return this->transform->get_fact_nnz(is_transposed?size()-id-1:id);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_nb_rows(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,0);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_nb_cols(const faust_unsigned_int id) const
		{
			return get_fact_dim_size(id,1);
		}

	template<typename FPP>
		unsigned int TransformHelper<FPP,Cpu>::get_fact_dim_size(const faust_unsigned_int id, unsigned short dim) const
		{
			//dim == 0 to get num of cols, >0 otherwise
			faust_unsigned_int rid; //real id
			if(is_transposed) {
				rid = size()-id-1;
				dim = !dim;
			}
			else {
				rid = id;
				//dim = dim;
			}
			Faust::MatGeneric<FPP,Cpu>* mat;
			if(rid == 0 || rid == this->size()-1)
				mat = this->transform->get_fact(rid, false);
			else
				mat = this->transform->get_fact(rid, false);
			if(dim)
				return mat->getNbCol();
			else
				return mat->getNbRow();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				const int** rowptr,
				const int** col_ids,
				const FPP** elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols);

		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				int* rowptr,
				int* col_ids,
				FPP* elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* = false*/) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, is_transposed ^ transpose);
			if(is_conjugate)
				Faust::conjugate(elts, *nnz);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				const FPP** elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, elts, num_rows, num_cols);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::get_fact(const faust_unsigned_int id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /* default to false */) const
		{
			this->transform->get_fact(is_transposed?this->size()-id-1:id, elts, num_rows, num_cols, is_transposed ^ transpose);
			if(is_conjugate)
				Faust::conjugate(elts,*num_cols*(*num_rows));
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
		{
			MatDense<FPP,Cpu> dense_factor;
			MatGeneric<FPP,Cpu>* factor_generic;
			factor_generic = this->transform->get_fact(is_transposed?size()-id-1:id);

			switch (factor_generic->getType())
			{
				case Dense :
					{
						MatDense<FPP,Cpu>* factor_dense_ptr = dynamic_cast<MatDense<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_dense_ptr); //copy
					}
					break;

				case Sparse :
					{
						MatSparse<FPP,Cpu>* factor_sparse_ptr = dynamic_cast<MatSparse<FPP,Cpu>* > (factor_generic);
						dense_factor = (*factor_sparse_ptr); //copy
					}
					break;

				default:
					handleError("Faust::TransformHelper", "get_fact : unknown type of the factor matrix");
			}
			delete factor_generic;
			if(is_transposed) dense_factor.transpose();
			if(is_conjugate) dense_factor.conjugate();
			return dense_factor;
		}

	template<typename FPP>
		bool Faust::TransformHelper<FPP,Cpu>::is_fact_sparse(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_sparse(is_transposed?size()-id-1:id);
		}


	template<typename FPP>
		bool Faust::TransformHelper<FPP,Cpu>::is_fact_dense(const faust_unsigned_int id) const
		{
			return this->transform->is_fact_dense(is_transposed?size()-id-1:id);
		}

	template<typename FPP>
	void TransformHelper<FPP, Cpu>::copy_slices(TransformHelper<FPP, Cpu> *th, const bool transpose /* default to false */)
	{
		this->slices[0].copy(th->slices[0]);
		this->slices[1].copy(th->slices[1]);
	}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP, Cpu>::slice(faust_unsigned_int start_row_id, faust_unsigned_int end_row_id,
				faust_unsigned_int start_col_id, faust_unsigned_int end_col_id)
		{
			Slice sr(start_row_id, end_row_id);
			Slice sc(start_col_id, end_col_id);
			Slice s[2];
			if(is_transposed)
			{
				s[0] = sc;
				s[1] = sr;
			}
			else
			{
				s[0] = sr;
				s[1] = sc;
			}
			return new TransformHelper<FPP, Cpu>(this, s);
		}

	template<typename FPP>
	const Transform<FPP, Cpu>* TransformHelper<FPP, Cpu>::eval_sliced_Transform()
		{
			if(is_sliced) {
				bool cloning_fact = false;
				std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) this->size());
				faust_unsigned_int size = this->size();
				MatGeneric<FPP,Cpu>* gen_fac, *first_sub_fac, *last_sub_fac;
				MatGeneric<FPP,Cpu>* gen_facs[size-2];
				gen_fac = this->transform->get_fact(0, cloning_fact);
				first_sub_fac = gen_fac->get_rows(slices[0].start_id, slices[0].end_id-slices[0].start_id);
				//		first_sub_fac->Display();
				if(cloning_fact)
					delete gen_fac;
				if(size > 1) {
					gen_fac = this->transform->get_fact(size-1, cloning_fact);
					last_sub_fac = gen_fac->get_cols(slices[1].start_id, slices[1].end_id-slices[1].start_id);
					//		std::cout << "---" << std::endl;
					//		last_sub_fac->Display();
					if(cloning_fact)
						delete gen_fac;
					factors.insert(factors.begin(), first_sub_fac);
					auto it = factors.begin();
					for(faust_unsigned_int i = 1; i < size-1; i++)
					{
						gen_facs[i-1] = this->transform->get_fact(i, cloning_fact);
						factors.insert(++it, gen_facs[i-1]);
					}
					factors.insert(factors.begin()+(size-1), last_sub_fac);
					factors.resize(size);
				}
				else { //only one factor
					last_sub_fac = first_sub_fac->get_cols(slices[1].start_id, slices[1].end_id-slices[1].start_id);
					delete first_sub_fac;
					factors[0] = last_sub_fac;
					factors.resize(1);
				}
				Transform<FPP, Cpu>* th = new Transform<FPP, Cpu>(factors, 1.0, false, cloning_fact);
				if(cloning_fact) {
					for(faust_unsigned_int i = 1; i < size-1; i++)
						delete gen_facs[i-1];
					delete first_sub_fac;
					delete last_sub_fac;
				}
				return th;
			}
			return this->transform.get(); //needed to return the stored Transform object ptr
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product() const {
			return this->transform->get_product(isTransposed2char(), is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::save_mat_file(const char* filename) const
		{
			this->transform->save_mat_file(filename, is_transposed, is_conjugate);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
		{
			return this->transform->spectralNorm(nbr_iter_max, threshold, flag);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::transpose()
		{
			return new TransformHelper<FPP,Cpu>(this, true, false);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::conjugate()
		{
			return new TransformHelper<FPP,Cpu>(this, false, true);
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::adjoint()
		{
			return new TransformHelper<FPP,Cpu>(this, true, true);
		}


	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isTransposed() const
		{
			return is_transposed;
		}

	template<typename FPP>
		bool TransformHelper<FPP,Cpu>::isConjugate() const
		{
			return is_conjugate;
		}

	template<typename FPP>
		const char TransformHelper<FPP,Cpu>::isTransposed2char() const
		{
			return is_transposed?'T':'N';
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normL1() const {
			return this->transform->normL1();
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normFro() const {
			return this->transform->normFro();
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::~TransformHelper() {
			// transform is deleted auto. when no TransformHelper uses it (no more weak refs left)
#ifdef FAUST_VERBOSE
			cout << "Destroying Faust::TransformHelper object." << endl;
#endif
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density) {
			if(!TransformHelper<FPP,Cpu>::seed_init) {
				srand(time(NULL)); //seed init needed for MatDense rand generation
				TransformHelper<FPP,Cpu>::seed_init = true;

			}
			// pick randomly the number of factors into {min_num_factors, ..., max_num_factors}
			std::uniform_int_distribution<int> num_fac_distr(min_num_factors, max_num_factors);
			std::uniform_int_distribution<int> dim_distr(min_dim_size, max_dim_size);
			std::uniform_int_distribution<int> bin_distr(0,1);
			unsigned int num_factors = num_fac_distr(generator);
			// create factors randomly respecting the RandFaustType asked and dims interval
			std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) num_factors);
			unsigned int num_rows, num_cols = dim_distr(generator);
			for(unsigned int i=0;i<num_factors;i++) {
				num_rows = num_cols;
				num_cols = dim_distr(generator);
				switch(t){
					case DENSE:
						factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, density);
						break;
					case SPARSE:
						factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, density);
						break;
					case MIXTE:
						if(bin_distr(generator))
							factors[i] = MatDense<FPP,Cpu>::randMat(num_rows, num_cols, density);
						else
							factors[i] = MatSparse<FPP,Cpu>::randMat(num_rows, num_cols, density);
						break;
					default:
						handleError("Faust::TransformHelper", "randFaust(): Unknown RandFaustType");
						break;

				}
				if(factors[i] == NULL) return NULL;
			}
			TransformHelper<FPP,Cpu>* randFaust = new TransformHelper<FPP, Cpu>(factors,1.0,false);
			return randFaust;
		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::hadamardFaust(unsigned int n)
		{
			if(n == 0) {
				MatDense<FPP,Cpu>* mat = new MatDense<FPP,Cpu>(1,1);
				(*mat)[0] = 1;
				vector<MatGeneric<FPP,Cpu>*> factors(1);
				factors[0] = mat;
				return new TransformHelper<FPP, Cpu>(factors,1.0,false);
			}
			unsigned int order = 1u << n;
			vector<int> col_ids(order), row_ids(order);
			vector<FPP> vals(order);
			unsigned int order_over_2 = order >> 1;
			unsigned int order_times_2 = order << 1;
			unsigned int i_times_2;

//			// init the permutation matrix

			for(unsigned int i=0; i < order_over_2; i++)
			{
				i_times_2 = i << 1;
				row_ids[i] = i_times_2;
				row_ids[i+order_over_2] = i_times_2+1;
				col_ids[i] = i;
				col_ids[i+order_over_2] = i + order_over_2;
				vals[i] = vals[i+order_over_2] = 1;
			}
//			cout << row_ids.size() << endl;
//			cout << col_ids.size() << endl;
//			cout << vals.size() << endl;
			MatSparse<FPP,Cpu> P(row_ids, col_ids, vals, order, order);
			P.update_dim();

			// init the base matrix
			int row_ptr[order+1];
			row_ptr[0] = 0;
			int bcol_ids[order_times_2];
			FPP bvalues[order_times_2];

//			cout << "row_ptr: ";
			for(int i = 1; i < order+1;i++)
			{
				row_ptr[i] = 2+row_ptr[i-1];
//				cout << row_ptr[i] << ",";
			}
//			cout << endl;

			bool parity = true; //row parity

			int col_id = 0;
//			cout << "bvalues: ";
//			cout << "bcol_ids: ";
			for(unsigned int i=0; i < order_times_2; i+=2)
			{

				if(parity) //row index is pair
					bvalues[i] = bvalues[i+1] = 1;
				else
				{
					bvalues[i+1] = -1;
					bvalues[i] = 1;
				}
//				cout << bvalues[i] << " " << bvalues[i+1];
				parity = ! parity;
				bcol_ids[i] = col_id;
				bcol_ids[i+1] = col_id+1;
//				cout << bcol_ids[i] << " " << bcol_ids[i+1] << " ";
				if(((i + 1) & 3u) == 3u) col_id+=2; // i+1 mod 4 == 0
			}
//			cout << endl;
			MatSparse<FPP,Cpu> B(order_times_2, order, order, bvalues, row_ptr, bcol_ids, false);
//			cout << "Faust::TransformHelper::hadamardFaust(), B:" << endl;
//			B.Display();

			std::vector<MatGeneric<FPP,Cpu>*> factors((size_t) n);
			MatSparse<FPP,Cpu>* factor = new MatSparse<FPP,Cpu>(order,order);
			factor->mat = B.mat*P.mat;
			factor->update_dim();

			factors[0] = factor;
			for(int i=1; i < n; i++)
				factors[i] = factor->Clone();

			TransformHelper<FPP,Cpu>* hadamardFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
			return hadamardFaust;

		}

	template<typename FPP>
		TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::fourierFaust(unsigned int n)
		{

			vector<MatGeneric<FPP,Cpu>*> factors(n+1);
			fft_factors(n, factors);
//			for(int i=0;i<n+1;i++)
//				factors[i]->Display();
			TransformHelper<FPP,Cpu>* fourierFaust = new TransformHelper<FPP, Cpu>(factors, 1.0, false, false);
			return fourierFaust;
//			return nullptr;

		}

	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));
}
