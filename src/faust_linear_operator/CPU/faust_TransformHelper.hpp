

namespace Faust {

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(const std::vector<MatGeneric<FPP,Cpu> *>& facts,
				const FPP lambda_, const bool optimizedCopy, const bool cloning_fact) : is_transposed(false),
																						is_conjugate(false),
																						is_sliced(false)
	{
		//TODO: transform is useless
		//		Transform<FPP,Cpu>* transform = new Transform<FPP,Cpu>(facts, lambda_, optimizedCopy, cloning_fact);
		this->transform = make_shared<Transform<FPP,Cpu>>(facts, lambda_, optimizedCopy, cloning_fact);
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper() : is_transposed(false), is_conjugate(false), is_sliced(false)
	{
		//TODO: transform is useless
		//		Transform<FPP,Cpu>* transform = new Transform<FPP,Cpu>();
		this->transform = make_shared<Transform<FPP,Cpu>>();
	}

	template<typename FPP>
		TransformHelper<FPP,Cpu>::TransformHelper(Transform<FPP,Cpu> &t) : is_transposed(false), is_conjugate(false), is_sliced(false)
	{
		//TODO: transform is useless
		//		Transform<FPP,Cpu>* transform = new Transform<FPP,Cpu>();
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
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A) const
		{
			MatDense<FPP,Cpu> M = (this->eval_sliced_Transform())->multiply(A, isTransposed2char());
			if(is_conjugate) M.conjugate();
			return M;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x) const
		{
			Vect<FPP,Cpu> v = (this->eval_sliced_Transform())->multiply(x, isTransposed2char());
			if(is_conjugate) v.conjugate();
			return v;
		}

	template<typename FPP>
		Vect<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const Vect<FPP,Cpu> x, const bool transpose)
		{
			is_transposed ^= transpose;
			Vect<FPP,Cpu> v = (this->eval_sliced_Transform())->multiply(x, isTransposed2char());
			is_transposed ^= transpose;
			return v;
		}


	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::multiply(const MatDense<FPP,Cpu> A, const bool transpose)
		{
			is_transposed ^= transpose;
			MatDense<FPP,Cpu> M = (this->eval_sliced_Transform())->multiply(A, isTransposed2char());
			is_transposed ^= transpose;
			return M;
		}

	template<typename FPP>
		TransformHelper<FPP, Cpu>* TransformHelper<FPP,Cpu>::multiply(TransformHelper<FPP, Cpu>* th_right)
		{
			return new TransformHelper<FPP,Cpu>(this, th_right);
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
			return (this->eval_sliced_Transform())->get_total_nnz();
		}

	template<typename FPP>
		faust_unsigned_int TransformHelper<FPP,Cpu>::size() const
		{
			return this->transform->size();
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::display() const
		{
			(this->eval_sliced_Transform())->Display(is_transposed);
		}

	template<typename FPP>
		string TransformHelper<FPP,Cpu>::to_string() const
		{
			return (this->eval_sliced_Transform())->to_string(is_transposed);
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_fact(faust_unsigned_int id) const
		{
			MatDense<FPP,Cpu> dense_factor;
			MatGeneric<FPP,Cpu>* factor_generic;
			if(id == 0 || id == this->size()-1)
				factor_generic = (this->eval_sliced_Transform())->get_fact(is_transposed?size()-id-1:id);
			else
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
		const Transform<FPP, Cpu>* TransformHelper<FPP, Cpu>::eval_sliced_Transform() const
		{
			if(is_sliced) {
//				if(this->sliced_transform) //TODO: possible opt.
//					return this->sliced_transform;
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
			return this->transform.get(); //need to return the stored Transform object ptr
		}

	template<typename FPP>
		MatDense<FPP,Cpu> TransformHelper<FPP,Cpu>::get_product() const {
			return (this->eval_sliced_Transform())->get_product(isTransposed2char(), is_conjugate);
		}

	template<typename FPP>
		void TransformHelper<FPP,Cpu>::save_mat_file(const char* filename) const
		{
			(this->eval_sliced_Transform())->save_mat_file(filename, is_transposed, is_conjugate);
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::spectralNorm(const int nbr_iter_max, double threshold, int &flag) const
		{
			return (this->eval_sliced_Transform())->spectralNorm(nbr_iter_max, threshold, flag);
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
			return (this->eval_sliced_Transform())->normL1();
		}

	template<typename FPP>
		double TransformHelper<FPP,Cpu>::normFro() const {
			return (this->eval_sliced_Transform())->normFro();
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

	template<typename FPP> bool TransformHelper<FPP,Cpu>::seed_init = false;
	template<typename FPP> std::default_random_engine TransformHelper<FPP,Cpu>::generator(time(NULL));

}
