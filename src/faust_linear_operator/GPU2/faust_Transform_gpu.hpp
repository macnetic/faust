namespace Faust
{

	template<typename FPP>
		RefManager Transform<FPP,GPU2>::ref_man([](void *fact)
				{
#ifdef DEBUG
				std::cout << "Transform delete_fact" << std::endl;
#endif
				delete static_cast<MatGeneric<FPP,GPU2>*>(fact);
				});

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
			//			fact->set_id(M.is_id()); //TODO: in MatGeneric<FPP,GPU2> add is_id/set_id
		}

	template<typename FPP>
		void Transform<FPP,GPU2>::update_total_nnz() const
		{
			//nothing to do get_total_nnz doesn't use any cache by now
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
			if(transpose)
			{
				auto t_smat = smat->clone();
				t_smat->transpose();
				t_smat->tocpu(d_outer_count_ptr, d_inner_ptr, d_elts, (int32_t*) num_rows, (int32_t*) num_cols, (int32_t*) nnz);
				delete t_smat;
			}
			else
			{
				smat->tocpu(d_outer_count_ptr, d_inner_ptr, d_elts, (int32_t*) num_rows, (int32_t*) num_cols, (int32_t*) nnz);
			}
		}


	template<typename FPP>
		Vect<FPP,GPU2> Transform<FPP,GPU2>::multiply(const Vect<FPP,GPU2>& x, const char opThis/*='N'*/)
		{
			MatDense<FPP, GPU2> out = this->multiply((const MatDense<FPP,GPU2>)x, opThis);
			Vect<FPP,GPU2> v_out;
			v_out.gpu_mat = out.gpu_mat;
			out.gpu_mat = nullptr; // avoid freeing v_out.gpu_mat when out of scope
			return v_out;
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
		faust_unsigned_int Transform<FPP, GPU2>::get_fact_nnz(const faust_unsigned_int id) const
		{
			return this->get_fact(id, false)->getNonZeros();
		}

}
