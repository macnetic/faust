

namespace Faust {

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
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
				const FPP lambda_, const bool optimizedCopy, const bool cloning_fact) : is_transposed(false),
																						is_conjugate(false),
																						is_sliced(false)
	{
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
			if(th->is_sliced) { // slice of slice
				this->slices[0].start_id = th->slices[0].start_id+s[0].start_id;
				this->slices[0].end_id = th->slices[0].start_id+s[0].end_id;
				this->slices[1].start_id = th->slices[1].start_id+s[1].start_id;
				this->slices[1].end_id = th->slices[1].start_id+s[1].end_id;
			}
			else {
				this->slices[0] = s[0];
				this->slices[1] = s[1];
			}
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

			TransformHelper<FPP,Cpu>* hadamardFaust = new TransformHelper<FPP, Cpu>(factors,1.0,false);
			return hadamardFaust;

		}

	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));
}
