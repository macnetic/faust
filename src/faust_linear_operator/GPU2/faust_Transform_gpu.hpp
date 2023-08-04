namespace Faust
{
	template<typename FPP>
		Transform<FPP,GPU2>::Transform() : dtor_delete_data(false), dtor_disabled(false), data(std::vector<MatGeneric<FPP, GPU2>*>()), total_nnz(0)
	{
	}

	template<typename FPP>
		faust_unsigned_int Transform<FPP,GPU2>::getNbRow() const // TODO: factorize with CPU code
		{
			if (size() != 0)
				return data[0]->getNbRow();
			throw std::runtime_error("Empty Transform");
		}

	template<typename FPP>
		faust_unsigned_int Transform<FPP,GPU2>::getNbCol() const // TODO: factorize with CPU code
		{
			if (size() != 0)
				return data[size()-1]->getNbCol();
			throw std::runtime_error("Empty Transform");
		}


	template<typename FPP>
		faust_unsigned_int Transform<FPP,GPU2>::size() const
		{
			return data.size();
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::push_back(const MatGeneric<FPP,GPU2>* M, bool copying/*=true*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
		{
			auto pushed_M = const_cast<MatGeneric<FPP,GPU2>*>(M);
			if((transpose || conjugate) && !copying)
				throw std::runtime_error("Transform<FPP,GPU2>::push_back(): copying argument must be true if any of transpose or conjugate argument is true.");
			if(copying)
			{
				pushed_M = M->clone();
				if(transpose && conjugate)
					pushed_M->adjoint();
				else if(transpose)
					pushed_M->transpose();
				else if(conjugate)
					pushed_M->conjugate();
			}
			data.push_back(const_cast<MatGeneric<FPP,GPU2>*>(pushed_M));
			if(!dtor_delete_data) ref_man.acquire(const_cast<MatGeneric<FPP,GPU2>*>(pushed_M));
			total_nnz += M->getNonZeros();
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::replace(const MatGeneric<FPP, GPU2>* new_mat, const faust_unsigned_int id)
		{
			// update the underlying gpu_mod array
			total_nnz -= data[id]->getNonZeros();
			// update local (wrapper) data
			if(dtor_delete_data)
				delete data[id];
			else
				ref_man.release(data[id]);
			data[id] = const_cast<MatGeneric<FPP,GPU2>*>(new_mat);
			total_nnz += new_mat->getNonZeros();
			if(! dtor_delete_data)
				ref_man.acquire(data[id]);
			//this->update_total_nnz();
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* Transform<FPP,GPU2>::get_fact(int32_t id, bool cloning_fact) const
		{
			if(cloning_fact)
				return data[id]->clone();
			else
				return data[id];
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::update(const MatGeneric<FPP, GPU2>& M, const faust_unsigned_int id)
		{
			auto fact = get_fact(id, false);
			auto fact_type = fact->getType();
			if(M.getType() != fact_type)
				throw std::runtime_error("The factor matrix to update is not of the same type (dense or sparse) as the input matrix.");
			total_nnz -= fact->getNonZeros();
			if(fact_type == Dense)
			{
				// fact to update is dense
				auto dfact = dynamic_cast<MatDense<FPP,GPU2>*>(fact);
				auto dM = dynamic_cast<const MatDense<FPP,GPU2>*>(&M);
				*dfact = *dM;
			}
			else
			{
				// fact to update is sparse
				auto sfact = dynamic_cast<MatSparse<FPP,GPU2>*>(fact);
				auto sM = dynamic_cast<const MatSparse<FPP,GPU2>*>(&M);
				*sfact = *sM;
			}
			total_nnz += M.getNonZeros();
			//TODO: other types
			//			fact->set_id(M.is_id()); //TODO: in MatGeneric<FPP,GPU2> add is_id/set_id
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::update_total_nnz()
		{
			total_nnz = 0;
			for(auto fac: data)
			{
				total_nnz += fac->getNonZeros();
			}
		}

	template<typename FPP>
		faust_unsigned_int Transform<FPP, GPU2>::get_total_nnz() const
		{
			return total_nnz;
		}

	template<typename FPP>
		void Transform<FPP, GPU2>::get_facts(std::vector<MatGeneric<FPP,GPU2>*> &factors, bool cloning_facts/*=true*/) const
		{
			for(int i=0;i < size(); i++)
			{
				factors.push_back(get_fact(i, cloning_facts));
			}
		}

	template<typename FPP>
		void Transform<FPP, GPU2>::operator=(const Transform<FPP,GPU2>& t)
		{
			this->clear();
			for(int i=0;i<t.size();i++)
			{
				this->push_back(t.get_fact(i, /*cloning*/false), /*copying*/ true);
			}
		}

	template<typename FPP>
		Transform<FPP,GPU2>::Transform(const Transform<FPP,GPU2> & t) : Transform<FPP,GPU2>()
	{
		*this = t;
	}

	template<typename FPP>
		void Transform<FPP,GPU2>::push_first(const MatGeneric<FPP,GPU2>* M, bool copying/*=true*/)
		{
			this->insert(0, M, copying);
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::pop_front()
		{
			this->erase(0);
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::pop_back()
		{
			this->erase(size()-1);
		}

	template<typename FPP>
		MatDense<FPP,GPU2> Transform<FPP,GPU2>::get_product(const char opThis/*='N'*/, const bool isConj/*=false*/) const
		{
			MatDense<FPP, GPU2> M;
			this->get_product(M, opThis, isConj);
			return M;
		}


	template<typename FPP>
		void Faust::Transform<FPP,GPU2>::get_product(Faust::MatDense<FPP,GPU2> &mat, const char opThis, const bool isConj)const //TODO: factorize with CPU code
		{

			if (size() == 0)
				throw std::runtime_error("get_product : empty Faust::Transform");

			if (opThis == 'N')
			{
				if(data.size() == 1)
				{
					auto end = this->size()-1;
					// just one matrix in the Faust, return a copy as dense matrix
					mat = data[end]->to_dense();
					if(isConj)
						mat.conjugate();
					return;
				}
				else
				{
					// at least two factors, compute the first product (of the last two matrices)
					// it avoids making a copy of the last factor
					gemm_gen(*data[this->size()-2], *data[this->size()-1], mat, FPP(1.0), FPP(0.0), 'N', 'N');
				}
				for (int i=this->size()-3; i >= 0; i--)
				{
					data[i]->multiply(mat,opThis);
				}
			}
			else
			{
				if(data.size() == 1)
				{
					// just one matrix in the Faust, return a transpose or transconjugate copy as dense matrix
					mat = data[0]->to_dense();
					if(opThis == 'H' || opThis == 'T' && isConj)
						mat.adjoint();
					else if(opThis == 'T')
						mat.transpose();
					else if(isConj)
						mat.conjugate();
					return;
				}
				else
				{
					// at least two factors, compute the first product (of the last two matrices)
					// it avoids making a copy of the first factor
					gemm_gen(*data[1], *data[0], mat, FPP(1.0), FPP(0.0), opThis, opThis);
				}
				for (int i=2; i < this->size(); i++)
				{
					data[i]->multiply(mat, opThis);
				}

			}

			if(isConj && opThis != 'H') mat.conjugate();

		}


	template<typename FPP>
		void Transform<FPP,GPU2>::multiply(const Transform<FPP,GPU2> & A)
		{
			if (A.size() != 0)
			{
				if (size() == 0)
				{
					(*this)=A;
				}
				else
				{
					if (getNbCol() != A.getNbRow())
					{
						throw std::runtime_error("Dimensions must agree");
					}
					for (int i=0;i<A.size();i++)
					{
						this->push_back(A.get_fact(i, /*cloning*/ false)); // push_back clones it afterward
					}
				}
			}
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::multiplyLeft(const Transform<FPP,GPU2> & A)
		{
			if (A.size() != 0)
			{
				if (size() == 0)
				{
					(*this)=A;
				}
				else
				{
					if (A.getNbCol() != getNbRow())
					{
						throw std::runtime_error("Dimensions must agree");
					}
					for (int i=A.size()-1;i>-1;i--)
					{
						this->push_first(A.get_fact(i, /*cloning*/ false)); // push_back clones it afterward
					}
				}
			}
		}

	template<typename FPP>
		bool Transform<FPP, GPU2>::is_fact_sparse(int id) const
		{
			return get_fact(id, /*cloning*/ false)->getType() == Sparse;
		}

	template<typename FPP>
		bool Transform<FPP, GPU2>::is_fact_dense(int id) const
		{
			return get_fact(id, /*cloning*/ false)->getType() == Dense;
		}

	template<typename FPP>
		bool Transform<FPP, GPU2>::is_fact_bsr(int id) const
		{
			return get_fact(id, /*cloning*/ false)->getType() == BSR;
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::get_fact(const faust_unsigned_int &id,
				FPP* elts,
				faust_unsigned_int* num_rows,
				faust_unsigned_int* num_cols,
				const bool transpose /*=false*/) const
		{
			if(! is_fact_dense(id))
				throw std::runtime_error("faust_Transform_gpu: this get_fact function signature is for MatDense only.");
			auto gen_mat = get_fact(id, /*cloning*/ false);
			auto dmat = dynamic_cast<MatDense<FPP,GPU2>*>(gen_mat);
			FPP tmp;
			*num_cols = gen_mat->getNbCol();
			*num_rows = gen_mat->getNbRow();
			/******* failsafe copy **/
			auto cpu_mdense = new MatDense<FPP,Cpu>(dmat->getNbRow(), dmat->getNbCol());
			dmat->tocpu(*cpu_mdense);
			memcpy(elts, cpu_mdense->getData(),cpu_mdense->getNbCol()*cpu_mdense->getNbRow()*sizeof(FPP));
			delete cpu_mdense;
			/********* efficient copy that should work but doesn't */
			// TODO: normally it must be possible to copy directly to elts but I noticed a segfault (and I couldn't figure out how to fix it)
			//			dmat->tocpu(elts);
			/***************************/
			if(transpose)
			{
				// transpose in-place
				// the matrix is in column-major order
				for(int j=0;j<*num_cols;j++)
					for(int i=0; i<*num_rows;i++)
					{
						tmp = elts[i**num_cols+j];
						elts[i**num_cols+j] = elts[j**num_rows+i];
						elts[j**num_rows+i] = tmp;
					}
				// swap num_cols and num_rows
				// (with only these 2 variables -- F2 arithmetic trick)
				*num_cols = *num_cols^*num_rows;
				*num_rows = *num_cols^*num_rows;
				*num_cols = *num_cols^*num_rows;
			}
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::get_fact(const faust_unsigned_int id,
				int* d_outer_count_ptr, int* d_inner_ptr, FPP* d_elts,
				faust_unsigned_int* nnz,
				faust_unsigned_int* num_rows, faust_unsigned_int* num_cols,
				bool transpose /* default to false */) const
		{
			if(! is_fact_sparse(id))
				throw std::runtime_error("faust_Transform_gpu: this get_fact function signature is for MatSparse only.");
			auto gen_mat = get_fact(id, /*cloning*/ false);
			auto smat = dynamic_cast<MatSparse<FPP,GPU2>*>(gen_mat);
			MatSparse<FPP,Cpu> cpu_smat;
			// warning: arguments are of type faust_unsigned_int in this function but tocpu below uses int32_t
			int32_t nr, nc, nz;
			if(transpose)
			{
				auto t_smat = smat->clone();
				t_smat->transpose();
				t_smat->tocpu(d_outer_count_ptr, d_inner_ptr, d_elts, &nr, &nc, &nz);
				delete t_smat;
			}
			else
			{
				smat->tocpu(d_outer_count_ptr, d_inner_ptr, d_elts, &nr, &nc, &nz);
			}
			*num_rows =  nr;
			*num_cols =  nc;
			*nnz = nz;
		}


	template<typename FPP>
		Vect<FPP,GPU2> Transform<FPP,GPU2>::multiply(const Vect<FPP,GPU2>& x, const char opThis/*='N'*/)
		{
			MatDense<FPP, GPU2> out = this->multiply(*dynamic_cast<const MatDense<FPP, GPU2>*>(&x), opThis);
			Vect<FPP,GPU2> v_out;
			v_out.gpu_mat = out.gpu_mat;
			out.gpu_mat = nullptr; // avoid freeing v_out.gpu_mat when out of scope
			return v_out;
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::scalarMultiply(const FPP& scalar, long int sid/*=-1*/) //TODO: factorize with CPU code (scalarMultiply)
		{
			// find smallest factor
			if (size() == 0)
				throw std::runtime_error("Empty Transform");

			if(sid < 0)
			{
				auto ssize = data[0]->getNbRow() * data[0]->getNbCol();
				for(auto i = 0;i < data.size(); i++)
				{
					auto fac = data[i];
					auto s = fac->getNbRow() * fac->getNbCol();
					if(s < ssize)
					{
						ssize = s;
						sid = i;
					}
				}
			}
			*(data[sid]) *= scalar;
		}

	template<typename FPP>
		MatGeneric<FPP,GPU2>* Transform<FPP,GPU2>::iterator::operator*() const
		{
			return container.get_fact(index, /*cloning_fact*/ false);
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator& Transform<FPP,GPU2>::iterator::operator++()
		{
			index++;
			return *this;
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator Transform<FPP,GPU2>::iterator::operator+(int arg)
		{
			iterator copy(*this);
			copy.index = this->index+arg;
			return copy;
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator Transform<FPP,GPU2>::iterator::operator-(int arg)
		{
			iterator copy(*this);
			copy.index = this->index-arg;
			return copy;
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator Transform<FPP,GPU2>::iterator::operator++(int)
		{
			iterator copy(*this);
			this->index++;
			return copy;
		}

	template<typename FPP>
		bool Transform<FPP,GPU2>::iterator::operator<(const Transform<FPP,GPU2>::iterator& it)
		{
			return this->index < it.index;
		}

	template<typename FPP>
		bool Transform<FPP,GPU2>::iterator::operator!=(const Transform<FPP,GPU2>::iterator& it)
		{
			return this->index != it.index;
		}

	template<typename FPP>
		Transform<FPP,GPU2>::iterator::iterator(const Transform<FPP, GPU2>& container, size_t index) : index(index), container(container)
	{
	}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator Transform<FPP,GPU2>::begin() const
		{
			return Transform<FPP,GPU2>::iterator(*this, 0);
		}

	template<typename FPP>
		typename Transform<FPP,GPU2>::iterator Transform<FPP,GPU2>::end() const
		{
			return Transform<FPP,GPU2>::iterator(*this, size());
		}

	template<typename FPP>
		void Transform<FPP, GPU2>::tocpu(Transform<FPP, Cpu>& cpu_transf) const
		{

			MatDense<FPP, GPU2>* gpu_mdense = nullptr;
			MatSparse<FPP, GPU2>* gpu_msparse = nullptr;
			MatBSR<FPP, GPU2>* gpu_mbsr = nullptr;
			for(auto gpu_mat: data)
			{
				if(gpu_mdense = dynamic_cast<MatDense<FPP, GPU2>*>(gpu_mat))
				{
					auto cpu_mdense = new MatDense<FPP,Cpu>(gpu_mdense->getNbRow(), gpu_mdense->getNbCol());
					gpu_mdense->tocpu(*cpu_mdense);
					cpu_transf.push_back(cpu_mdense, false, false);
				}
				else if(gpu_msparse = dynamic_cast<MatSparse<FPP, GPU2>*>(gpu_mat))
				{
					auto cpu_msparse = new MatSparse<FPP,Cpu>();
					cpu_msparse->resize(gpu_msparse->getNonZeros(), gpu_msparse->getNbRow(), gpu_msparse->getNbCol());
					gpu_msparse->tocpu(*cpu_msparse);
					cpu_transf.push_back(cpu_msparse, false, false, false);
				}
				else if(gpu_mbsr = dynamic_cast<MatBSR<FPP, GPU2>*>(gpu_mat))
				{
					auto cpu_mbsr = new MatBSR<FPP,Cpu>();
					gpu_mbsr->tocpu(*cpu_mbsr);
					cpu_transf.push_back(cpu_mbsr, false, false, false);
				}
				else
					throw std::runtime_error("Invalid matrix pointer");
			}
		}

	//TODO: refactor to generic CPU/GPU code (using if needed a non-member function on Transform<FPP, DEV>)
	template<typename FPP>
		MatDense<FPP,GPU2> Transform<FPP,GPU2>::multiply(const MatDense<FPP,GPU2> &A, const char opThis) /*const*/ //TODO: should be const
		{
			if (size() == 0)
				handleWarning("Transform<FPP,GPU2> : multiply : empty Transform<FPP,GPU2>");

			MatDense<FPP,GPU2> mat(A);


			if (opThis == 'N')
				for (int i=this->size()-1; i >= 0; i--)
					data[i]->multiply(mat, opThis);
			else
				for (int i=0; i < this->size(); i++)
					data[i]->multiply(mat, opThis);

			return mat;
		}

	template<typename FPP>
		Transform<FPP, Cpu> Transform<FPP, GPU2>::tocpu() const
		{
			Transform<FPP, Cpu> cpu_transf;
			tocpu(cpu_transf);
			return cpu_transf;
		}

	template<typename FPP>
		void Transform<FPP, GPU2>::save_mat_file(const char* filename, const bool transpose, const bool conjugate) const
		{
			Transform<FPP,Cpu> cpu_transf;
			this->tocpu(cpu_transf);
			cpu_transf.save_mat_file(filename, transpose, conjugate);
		}

	template<typename FPP>
		Real<FPP> Transform<FPP, GPU2>::normL1(const bool transpose /* = false */, const bool full_array/*=true*/, const int batch_size/*=1*/) const
		{
			double norm;
			MatDense<FPP, GPU2> full = get_product(transpose?'T':'N');
			norm = std::abs(full.normL1(/*transpose*/)); //transpose not necessary because full is already transposed if needed
			return norm;
		}


	template<typename FPP>
		faust_unsigned_int Transform<FPP, GPU2>::get_fact_nnz(const faust_unsigned_int id) const //TODO: factorize with CPU code
		{
			return this->get_fact(id, false)->getNonZeros();
		}


	template<typename FPP>
		void Faust::Transform<FPP, GPU2>::transpose() //TODO: factorize with CPU code
		{
			std::reverse(data.begin(), data.end());
			for (int i=0;i<size();i++)
				data[i]->transpose();
		}


	template<typename FPP>
		void Faust::Transform<FPP,GPU2>::adjoint()
		{
			transpose();
			for(auto it = data.begin(); it != data.end(); it++){
				(*it)->conjugate();
			}
		}

	// compute the largest eigenvalue of A, A must be positive semi-definite
	template<typename FPP>
		template<typename FPP2>
		FPP Transform<FPP, GPU2>::power_iteration(const faust_unsigned_int nbr_iter_max, const FPP2 threshold, int & flag, const bool rand_init/*=true*/) const
		{
			auto A = *this;

			const int nb_col = A.getNbCol();
			int i = 0;
			flag = 0;

			if (nbr_iter_max <= 0)
			{
				handleError("linear_algebra "," power_iteration :  nbr_iter_max <= 0");
			}
			if (nb_col != A.getNbRow())
			{
				handleError("linear_algebra "," power_iteration : Faust::Transform<FPP,GPU2> must be a square matrix");
			}
			Faust::Vect<FPP,GPU2> xk(nb_col);
			if(rand_init)
			{
				srand(0xF4+'u'+57); // always the same seed, not standard but allows reproducibility
				xk.setRand(); // the main objective is to make very unlikely to be orthogonal to the main eigenvector
			}
			else
				xk.setOnes();
			Faust::Vect<FPP,GPU2> xk_norm(nb_col);
			FPP lambda_old = 1.0;
			FPP lambda = 0.0;
			FPP alpha = 1.0;
			FPP beta = 0.0;
			while((Faust::fabs(lambda_old-lambda)> Faust::fabs(threshold) || Faust::fabs(lambda) <= Faust::fabs(threshold)) && i<nbr_iter_max)
			{
				i++;
				lambda_old = lambda;
				xk_norm = xk;
				xk_norm.normalize();
				xk = A.multiply(xk_norm);
				//				if(xk.isZero()) // A is most likely zero, lambda is zero //TODO
				//				{
				//					std::cerr << "WARNING: power_iteration product Ax leads to zero vector, A is most likely zero, lambda should be zero too." << std::endl;
				//					return FPP(0);
				//				}
				lambda = xk_norm.dot(xk);
				//std::cout << "i = " << i << " ; lambda=" << lambda << std::endl;
			}
			flag = (i<nbr_iter_max)?i:-1;

			return lambda;
		}


	template<typename FPP>
		Real<FPP> Transform<FPP,GPU2>::spectralNorm(const int nbr_iter_max, float threshold, int &flag) const // TODO: factorize with CPU code
		{
			if (size() == 0)
			{
				return 1; // TODO: why?
			}else
			{
				// if(this->is_zero) // The Faust is zero by at least one of its factors
				//return 0; //TODO
				// The Faust can still be zero (without any of its factor being)
				// this case will be detected in power_iteration
				Transform<FPP,GPU2> AtA((*this));
				AtA.adjoint();
				if (getNbCol() < getNbRow())
				{
					AtA.multiply(*this);
				}else
				{
					AtA.multiplyLeft(*this);
				}
				FPP maxAbsValue = std::sqrt(AtA.power_iteration(nbr_iter_max, threshold, flag));
				return absValue(maxAbsValue);

			}
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::Display(const bool transpose /* default to false */,const bool displaying_small_mat_elts /*false by default*/) const //TODO: factorize with CPU code
		{
			std::cout << to_string(transpose,displaying_small_mat_elts);
		}


	template<typename FPP>
		std::string Transform<FPP,GPU2>::to_string(const bool transpose /* default to false */, const bool displaying_small_mat_elts/* false by default */) const //TODO: factorize with CPU code
		{
			std::ostringstream str;

			if (size() == 0)
				str<<"empty Faust"<<std::endl;
			else
			{
				str<<"Faust size ";
				if(transpose)
					str << this->getNbCol() << "x" << this->getNbRow();
				else
					str << this->getNbRow()<<"x"<<this->getNbCol();
				str <<", density "<<1.0/getRCG()<< ", nnz_sum "<<this->get_total_nnz() << ", " << size() << " factor(s): "<< std::endl;
				int j;
				for (int i=0 ; i<size() ; i++)
				{
					if(transpose)
						j = size()-1-i;
					else
						j = i;
					str << "- FACTOR " << i;
					str << data[j]->to_string(transpose, displaying_small_mat_elts);
				}
			}
			return str.str();
		}

	template<typename FPP>
		void Transform<FPP, GPU2>::clear()
		{
			for (int i=0;i<data.size();i++)
			{
				if(dtor_delete_data)
					delete data[i];
				else
				{
					ref_man.release(data[i]);
				}
			}
			data.resize(0);
			total_nnz = 0;
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::insert(int32_t id, const MatGeneric<FPP,GPU2>* M, bool copying/*=true*/)
		{
			auto ins_M = const_cast<MatGeneric<FPP,GPU2>*>(M);
			if(copying)
				ins_M = M->clone();
			data.insert(data.begin()+id, ins_M);
			if(!dtor_delete_data) ref_man.acquire(ins_M);
			total_nnz += M->getNonZeros();
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::erase(int32_t id)
		{
			if(id < 0 || id >= size()) throw std::runtime_error("Transform<FPP, GPU2>: erase id is out of range");
			total_nnz -= data[id]->getNonZeros();
			if(!dtor_delete_data) ref_man.release(*(data.begin()+id));
			data.erase(data.begin()+id);
		}

	template<typename FPP>
		Transform<FPP,GPU2>::Transform(const std::vector<MatGeneric<FPP,GPU2>*> &factors, const FPP lambda_ /*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact/*=true*/) : Transform()
	{
		//TODO: take optional arguments into account
		GPUModHandler::get_singleton()->check_gpu_mod_loaded();
		for(auto m: factors)
			push_back(m, cloning_fact);
	}

	template<typename FPP>
		Transform<FPP, GPU2>::~Transform()
		{
			if(! this->dtor_disabled)
			{
				for (auto fac: data)
					if(this->dtor_delete_data)
						delete fac;
					else
						ref_man.release(fac);
			}
		}

}
