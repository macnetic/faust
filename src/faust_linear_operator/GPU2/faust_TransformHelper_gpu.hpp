#include "faust_MatButterfly_gpu.h"
#include "faust_MatPerm_gpu.h"

namespace Faust
{
    template<typename FPP,FDevice DEVICE> class Transform;

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper() : TransformHelperGen<FPP,GPU2>()
    {
    }

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper(const TransformHelper<FPP,GPU2>& th, bool transpose, bool conjugate) : TransformHelperGen<FPP,GPU2>(&th, transpose, conjugate)
    {
    }

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_/*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/) : TransformHelper<FPP,GPU2>()
    {

        //if lambda is not 1.0 a factor will be multiplied and so it needs to be copied to preserve the original that could be used elsewhere
        // in an optimization purpose, the smallest factor is copied
        int min_size_id = 0;
        if(lambda_ != FPP(1.0))
        {
            std::vector<int> fact_ids(facts.size());
            int i = -1;
            std::generate(fact_ids.begin(), fact_ids.end(), [&i](){return ++i;});
            std::vector<int>::iterator result = std::min_element(fact_ids.begin(), fact_ids.end(), [&facts](const int &a, const int &b){return facts[a]->getNBytes() < facts[b]->getNBytes();});
            min_size_id = std::distance(fact_ids.begin(), result);
        }
        for(int i=0; i < facts.size(); i++)
        {
            if(i == min_size_id)
                this->push_back(facts[min_size_id], false, cloning_fact || lambda_ != (FPP) 1.0);
            else
                this->push_back(facts[i], false, cloning_fact);
        }
        this->transform->multiply(lambda_, min_size_id);
    }

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper(const TransformHelper<FPP,Cpu>& cpu_t, int32_t dev_id/*=-1*/, void* stream/*=nullptr*/) : TransformHelper<FPP,GPU2>()
    {
        for(auto cpu_fact: cpu_t)
            this->push_back(cpu_fact, false, dev_id, stream);
        this->is_transposed = cpu_t.is_transposed;
        this->is_conjugate = cpu_t.is_conjugate;
        //TODO: slice etc.
    }

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper(TransformHelper<FPP,GPU2>* th, Slice s[2]): TransformHelper<FPP,GPU2>()
    {
        this->init_sliced_transform(th, s);
    }

    template<typename FPP>
        TransformHelper<FPP,GPU2>::TransformHelper(TransformHelper<FPP,GPU2>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols): TransformHelper<FPP,GPU2>()
    {
        this->init_fancy_idx_transform(th, row_ids, num_rows, col_ids, num_cols);
    }

#ifndef IGNORE_TRANSFORM_HELPER_VARIADIC_TPL
    template<typename FPP>
        template<typename ... GList>
        TransformHelper<FPP,GPU2>::TransformHelper(GList& ... t): TransformHelper<FPP,GPU2>()
        {
            this->push_back_(t...);
        }
#endif



    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            //optimizedCopy is ignored because not handled yet by Transform<FPP,GPU2> // TODO ? (it's not used by wrappers anyway)
            this->transform->push_back(M, copying, transpose, conjugate);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy/*=false*/, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            MatGeneric<FPP,GPU2>* gpu_M = nullptr;
            const MatDense<FPP,Cpu>* cpu_dM = nullptr;
            const MatSparse<FPP,Cpu>* cpu_sM = nullptr;
            const MatBSR<FPP,Cpu>* cpu_bM = nullptr;
            const MatButterfly<FPP,Cpu>* cpu_bf = nullptr;
            const MatPerm<FPP,Cpu>* cpu_p = nullptr;
            if(nullptr != (cpu_dM = dynamic_cast<const MatDense<FPP,Cpu>*>(M)))
            {
                auto gpu_dM = new MatDense<FPP,GPU2>(M->getNbRow(), M->getNbCol(), cpu_dM->getData());
                gpu_M = gpu_dM;
            }
            else if(nullptr != (cpu_sM = dynamic_cast<const MatSparse<FPP,Cpu>*>(M)))
            {
                auto gpu_sM = new MatSparse<FPP,GPU2>(M->getNbRow(), M->getNbCol(), cpu_sM->getNonZeros(), cpu_sM->getValuePtr(), cpu_sM->getRowPtr(), cpu_sM->getColInd(), dev_id);
                gpu_M = gpu_sM;
            }
            else if(nullptr != (cpu_bM = dynamic_cast<const MatBSR<FPP,Cpu>*>(M)))
            {
                auto gpu_bM = new MatBSR<FPP,GPU2>(*cpu_bM);
                gpu_M = gpu_bM;
            }
            else if(nullptr != (cpu_bf = dynamic_cast<const MatButterfly<FPP,Cpu>*>(M)))
            {
                auto gpu_bf = new MatButterfly<FPP,GPU2>(*cpu_bf);
                gpu_M = gpu_bf;
            }
            else if(nullptr != (cpu_p = dynamic_cast<const MatPerm<FPP,Cpu>*>(M)))
            {
                auto gpu_p = new MatPerm<FPP,GPU2>(*cpu_p);
                gpu_M = gpu_p;
            }
	    else throw std::runtime_error("Unhandled CPU matrix type for conversion to GPU TransformHelper");
	    // won't work  anyway for MatPerm/Butterfly because of get_gpu_mat_ptr exception (but it will when Transform GPU2 will be more idenpendent from gpu_mod cuMatArray)
            this->transform->push_back(gpu_M, false);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back(const FPP* data, const int nrows, const int ncols, const bool optimizedCopy/*=false*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            auto dense_mat = new MatDense<FPP,GPU2>(nrows, ncols, data, /* no_alloc */ false);
            auto copying = transpose || conjugate || optimizedCopy;
            this->push_back(dense_mat, optimizedCopy, copying, transpose, conjugate); // optimizedCopy not supported on GPU2
            if(copying) delete dense_mat;
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back(const FPP* data, const int* row_ptr, const int* id_col, const int nnz, const int nrows, const int ncols, const bool optimizedCopy/*=false*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            auto sparse_mat = new MatSparse<FPP,GPU2>(nrows, ncols, nnz, data, row_ptr, id_col);
            auto copying = transpose || conjugate || optimizedCopy;
            this->push_back(sparse_mat, optimizedCopy, copying, transpose, conjugate); // optimizedCopy not supported on GPU2
            if(copying) delete sparse_mat;
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back(const FPP* bdata, const int* brow_ptr, const int* bcol_inds, const int nrows, const int ncols, const int bnnz, const int bnrows, const int bncols, const bool optimizedCopy/*=false*/, const bool transpose/*=false*/, const bool conjugate/*=false*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			auto bsr_mat = new MatBSR<FPP, GPU2>(nrows, ncols, bnrows, bncols, bnnz, bdata, brow_ptr, bcol_inds);
			auto copying = optimizedCopy || transpose || conjugate;
			this->push_back(bsr_mat, copying, transpose, conjugate);
			if(copying) delete bsr_mat;
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_first(const MatGeneric<FPP,GPU2>* M, const bool optimizedCopy/*=false*/, const bool copying/*=true*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            return this->transform->push_first(M, copying);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::display() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            this->transform->Display(this->is_transposed);
        }

    template<typename FPP>
        std::string TransformHelper<FPP,GPU2>::to_string() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->to_string(this->is_transposed);
        }

    template<typename FPP>
        template<typename Head, typename ... Tail>
        void TransformHelper<FPP,GPU2>::push_back_(Head& h, Tail&... t)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            for(auto it=h.begin(); it < h.end(); it++)
            {
                auto f = *it;
                this->push_back(f, false, false);
            }
            this->push_back_(t...);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::push_back_()
        {
            // do nothing, here just for empty tail of above function
        }

    template<typename FPP>
        MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::get_product(int prod_mod/*=-1*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			if(this->mul_order_opt_mode == DYNPROG)
			{
				std::vector<Faust::MatGeneric<FPP,GPU2>*> data = this->transform->data;
				if(this->is_transposed)
					std::reverse(data.begin(), data.end());
				auto P = std::move(dynprog_multiply(data, this->isTransposed2char()));
				return P;
			}
			else
				return this->transform->get_product(this->isTransposed2char(), this->is_conjugate);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_product(MatDense<FPP,GPU2>& M, int prod_mod/*=-1*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            this->transform->get_product(M, this->isTransposed2char(), this->is_conjugate);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_product(MatDense<FPP,Cpu>& M, int prod_mod/*=-1*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            MatDense<FPP,GPU2> gpuM;
            this->get_product(gpuM, prod_mod);
            M = gpuM.tocpu();
        }

    template<typename FPP>
        Real<FPP> TransformHelper<FPP,GPU2>::normFro(const bool full_array/*=true*/, const int batch_size/*=1*/) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->get_product().norm();
        }

    template<typename FPP>
        Real<FPP> TransformHelper<FPP,GPU2>::normL1(const bool full_array/*=true*/, const int batch_size/*=1*/) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->normL1(this->is_transposed, full_array, batch_size);
        }

    template<typename FPP>
        Real<FPP> TransformHelper<FPP,GPU2>::normInf(const bool full_array/*=true*/, const int batch_size/*=1*/) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->normL1(!this->is_transposed, full_array, batch_size);
        }

    template<typename FPP>
        faust_unsigned_int TransformHelper<FPP,GPU2>::size() const
        {
            return this->transform->size();
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::update_total_nnz() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            this->transform->update_total_nnz();
        }

    template<typename FPP>
        Real<FPP> TransformHelper<FPP,GPU2>::spectralNorm(int32_t nb_iter_max, float threshold, int& flag)
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->spectralNorm(nb_iter_max, threshold, flag);
        }
    template<typename FPP>
        FPP TransformHelper<FPP,GPU2>::power_iteration(const faust_unsigned_int nb_iter_max, const Real<FPP>& threshold, int& flag)
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->power_iteration(nb_iter_max, threshold, flag);
        }

    template<typename FPP>
        const MatGeneric<FPP,GPU2>* TransformHelper<FPP,GPU2>::get_gen_fact(const faust_unsigned_int id) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->get_fact(id, false);
        }

    template<typename FPP>
        MatGeneric<FPP,GPU2>* TransformHelper<FPP,GPU2>::get_gen_fact_nonconst(const faust_unsigned_int id) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->get_fact(id, false);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::update(const MatGeneric<FPP, GPU2>& M,const faust_unsigned_int id)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            return this->transform->update(M, id);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::replace(const MatGeneric<FPP, GPU2>* M,const faust_unsigned_int id)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            return this->transform->replace(M, id);
        }

    template<typename FPP>
        void TransformHelper<FPP, GPU2>::convertToSparse()
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            const MatDense<FPP,GPU2> * mat_dense;
            const MatSparse<FPP,GPU2> * mat_sparse;
            for(int i=0;i<this->size();i++)
            {
                if(mat_dense = dynamic_cast<const MatDense<FPP,GPU2>*>(this->get_gen_fact(i)))
                {
                    mat_sparse = new MatSparse<FPP,GPU2>(*mat_dense);
                    this->replace(mat_sparse, i);
                }
            }
        }

    template<typename FPP>
        void TransformHelper<FPP, GPU2>::convertToDense()
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            const MatDense<FPP,GPU2> * mat_dense;
            const MatSparse<FPP,GPU2> * mat_sparse;
            for(int i=0;i<this->size();i++)
            {
                if(mat_sparse = dynamic_cast<const MatSparse<FPP,GPU2>*>(this->get_gen_fact(i)))
                {
                    mat_dense = new MatDense<FPP,GPU2>(*mat_sparse);
                    this->replace(mat_dense, i);
                }
            }
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::multiply(const TransformHelper<FPP,GPU2>* right)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(right)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(right)->eval_fancy_idx_Transform();
            // The goal is to minimize the number of factors copied (but maybe the criterion should be the sum of the size of these factors rather than their number)
            //			std::cout << "===" << this->is_transposed << std::endl;
            //			this->display();
            //			std::cout << "first fact:" << std::endl;
            //			this->get_gen_fact(0)->Display();
            //			std::cout << "===" << right->is_transposed << std::endl;
            //			right->display();
            //			std::cout << "===" << std::endl;
            bool copying_this = false;
            bool copying_right = false;
            bool conj_diff = this->is_conjugate != right->is_conjugate;
            bool trans_diff = this->is_transposed != right->is_transposed;
            bool transconj_diff = conj_diff || trans_diff;
            bool out_transp = false, out_conj = false;
            bool transp_this = false, transp_right = false, conj_this = false, conj_right = false;
            if(transconj_diff)
                if(this->size() < right->size())
                {
                    copying_this = true;
                    out_transp = trans_diff && right->is_transposed;
                    out_conj = conj_diff && right->is_conjugate;
                    transp_this = trans_diff;
                    conj_this = conj_diff;
                }
                else
                {
                    copying_right = true;
                    out_transp = trans_diff && this->is_transposed;
                    out_conj = conj_diff && this->is_conjugate;
                    transp_right = trans_diff;
                    conj_right = conj_diff;
                }
            auto mul_faust = new TransformHelper<FPP,GPU2>();
            //			std::cout << "transp_this: " << transp_this << " conj_this: " << conj_this << std::endl;
            //			std::cout << "transp_right: " << transp_right << " conj_right: " << conj_right << std::endl;
            //			std::cout << "out_trans: " << out_transp << " out_conj: " << out_conj << std::endl;
            std::function<void()> push_right_facts = [&out_transp, &transp_right, &mul_faust, &right, &copying_right, &conj_right]()
            {
                if(transp_right)
                    for(int i=right->size()-1; i>= 0; i--)
                        mul_faust->push_back(right->get_gen_fact(i), /*OptimizedCopy*/ false, copying_right, /*transp_right*/ true, conj_right);
                else
                    for(auto f: *right)
                        mul_faust->push_back(f, /*OptimizedCopy*/ false, copying_right, /*transp_right*/ false, conj_right);
            };
            std::function<void()> push_this_facts = [&transp_this, &mul_faust, this, &copying_this, &conj_this]()
            {
                if(transp_this)
                    for(int i=size()-1; i>= 0; i--)
                        mul_faust->push_back(get_gen_fact(i), /*OptimizedCopy*/ false, copying_this, /*transp_this*/ true, conj_this);
                else
                    for(auto f: *this)
                        mul_faust->push_back(f, /*OptimizedCopy*/ false, copying_this, /*transp_this*/ false, conj_this);
            };
            if(out_transp)
            {
                push_right_facts();
                push_this_facts();
            }
            else
            {
                push_this_facts();
                push_right_facts();
            }
            mul_faust->is_transposed = out_transp;
            mul_faust->is_conjugate = out_conj;
            return mul_faust;
        }

    template<typename FPP>
        MatDense<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,GPU2> &A)
        {
            MatDense<FPP,GPU2> M = this->transform->multiply(A, this->isTransposed2char());
            return M;
        }

    template<typename FPP>
        MatDense<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::MatDense<FPP,Cpu> &A)
        {
            MatDense<FPP,GPU2> M = this->multiply(MatDense<FPP,GPU2>(A));
            return M.tocpu();
        }

    template<typename FPP>
        Vect<FPP,Cpu> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,Cpu> &A)
        {
			//TODO: delegate to function below
            Vect<FPP,GPU2> gpu_A(A.size(), A.getData());
            Vect<FPP,GPU2> v = this->multiply(gpu_A);
            return v.tocpu();
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::multiply(const FPP* cpu_in_buf, FPP* cpu_out_buf)
        {
			if(this->is_sliced)
			{
				this->sliceMultiply(this->slices, cpu_in_buf, cpu_out_buf, 1);
			}
			else if(this->is_fancy_indexed)
			{
				size_t id_lens[2] = {this->fancy_num_rows, this->fancy_num_cols};
				this->indexMultiply(this->fancy_indices, id_lens, cpu_in_buf, 1, cpu_out_buf);
			}
			else
			{
				int32_t in_vec_size = this->getNbCol();
				Vect<FPP,GPU2> gpu_A(in_vec_size, cpu_in_buf);
				Vect<FPP,GPU2> v = this->multiply(gpu_A); //TODO: handle transpose and conjugate
				v.tocpu(cpu_out_buf);
			}
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::multiply(const FPP* cpu_x_buf, int x_ncols, FPP* cpu_out_buf)
        {
			if(this->is_sliced)
			{
#if DEBUG
				std::cout << "calling sliceMultiply on GPU Faust" << std::endl;
#endif
				this->sliceMultiply(this->slices, cpu_x_buf, cpu_out_buf, x_ncols);
			}
			else if(this->is_fancy_indexed)
			{
				size_t id_lens[2] = {this->fancy_num_rows, this->fancy_num_cols};
				this->indexMultiply(this->fancy_indices, id_lens, cpu_x_buf, x_ncols, cpu_out_buf);
			}
			else
			{
				int32_t x_nrows;
				if(this->is_transposed)
					x_nrows = this->transform->getNbRow();
				else
					x_nrows = this->transform->getNbCol();
				MatDense<FPP,GPU2> gpu_x(x_nrows, x_ncols, cpu_x_buf, false);
				MatDense<FPP,GPU2> gpu_M = this->multiply(gpu_x); //TODO: handle transpose and conjugate
																  // TODO: fix this function, it works until here then it segfaults or gives a cuda error with tocpu (even if I use a cpu matdense set locally)
				gpu_M.tocpu(cpu_out_buf, nullptr);
			}
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::multiply(const FPP& a)
        {
            const vector<MatGeneric<FPP,GPU2>*>& vec = this->transform->data; //TransformHelper is a friend class of Transform // we can access private attribute data
            TransformHelper<FPP,GPU2>* th = new TransformHelper<FPP,GPU2>(vec, a, false, false, true);
            th->copy_transconj_state(*this);
            th->copy_slice_state(*this);
			th->copy_fancy_idx_state(*this);
            return th;
        }

    template<typename FPP>
        Vect<FPP,GPU2> TransformHelper<FPP,GPU2>::multiply(const Faust::Vect<FPP,GPU2>& a)
        {
            Vect<FPP,GPU2> v = this->transform->multiply(a, this->isTransposed2char());
            return v;
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::pop_front()
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            return this->transform->pop_front();
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::pop_back()
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            return this->transform->pop_back();
        }

    template<typename FPP>
        void Faust::TransformHelper<FPP,GPU2>::pack_factors(faust_unsigned_int start_id, faust_unsigned_int end_id,const int mul_order_opt_mode/*=DEFAULT*/)
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            if(start_id < 0 || start_id >= size())
                throw out_of_range("start_id is out of range.");
            if(end_id < 0 || end_id >= size())
                throw out_of_range("end_id is out of range.");
            Faust::MatDense<FPP,GPU2> * packed_fac = nullptr;
            if(end_id == start_id)
            {
                //nothing to do except converting to MatDense if start_id
                //factor is a MatSparse
                packed_fac = dynamic_cast<Faust::MatDense<FPP,GPU2>*>(*(begin()+start_id));
                if(packed_fac == nullptr)
                {// factor start_id is not at MatDense, convert it
                    packed_fac = new MatDense<FPP,GPU2>(*dynamic_cast<Faust::MatSparse<FPP,GPU2>*>(*(begin()+start_id)));
                }
                else
                {
                    return; //no change
                }
            }
            else
            {
                // we have to multiply factors from start_id to end_id into one matrix
                // simple way to do, 1) create a overhead-free TransformHelper with these factors
                // 2) call get_product() to override the start_id factors with the result on the end
                // 3) erase factors from start_id to end_id and insert packed factor too to replace them (that's Transform object responsibility).
                // 1)
                std::vector<Faust::MatGeneric<FPP,GPU2>*> topack_factors;
                for(int i=start_id;i <= end_id; i++)
                    topack_factors.push_back(get_gen_fact_nonconst(i));
                //				std::vector<Faust::MatGeneric<FPP,GPU2>*> topack_factors(begin()+start_id, begin()+end_id+1);
                Faust::TransformHelper<FPP,GPU2> t(topack_factors, 1.0, false, false, false);
                //TODO: not yet implemented for GPU2
                //				t.set_FM_mul_mode(mul_order_opt_mode);
                // 2)
                packed_fac = new MatDense<FPP,GPU2>(t.get_product());
            }
            // 3)
            faust_unsigned_int i = end_id;
            while(i>=start_id)
            {
                this->transform->erase(i);
                if(i == 0) break;
                i--;
            }
            this->transform->insert(start_id, packed_fac, false);
        }

    template<typename FPP>
        typename Transform<FPP,GPU2>::iterator TransformHelper<FPP,GPU2>::begin() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->begin();
        }

    template<typename FPP>
        typename Transform<FPP,GPU2>::iterator TransformHelper<FPP,GPU2>::end() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->end();
        }


    template<typename FPP>
        void TransformHelper<FPP,GPU2>::operator=(TransformHelper<FPP,GPU2>& th)
        {
            copy_state(th); // it copies the underlying Transform object too
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::save_mat_file(const char* filepath) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            this->transform->save_mat_file(filepath, this->is_transposed, this->is_conjugate);
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2> TransformHelper<FPP,GPU2>::read_from_mat_file(const char* filepath)
        {
            throw std::runtime_error("A .mat file is always read from CPU code not GPU's.");
        }

    template<typename FPP>
        int TransformHelper<FPP,GPU2>::get_mat_file_type(const char* filepath)
        {
            return TransformHelper<FPP,Cpu>::get_mat_file_type(filepath);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::tocpu(TransformHelper<FPP,Cpu>& cpu_transf) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            //TODO: tocpu support of arguments transpose and conjugate
            auto t = this->transform->tocpu();
            for(auto fac: t)
            {
                cpu_transf.push_back(fac, false, false);
            }
            cpu_transf.is_transposed = this->is_transposed;
            cpu_transf.is_conjugate = this->is_conjugate;
			// no need to handle slicing and indexing (evaluated above)
        }

    template<typename FPP>
        TransformHelper<FPP,Cpu>* TransformHelper<FPP,GPU2>::tocpu() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            auto cpu_t = new TransformHelper<FPP,Cpu>();
            tocpu(*cpu_t);
            return cpu_t;
        }


    template<typename FPP>
        faust_unsigned_int TransformHelper<FPP,GPU2>::get_total_nnz() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            return this->transform->get_total_nnz();
        }


    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::normalize(const int meth /* 1 for 1-norm, 2 for 2-norm (2-norm), -1 for inf-norm */) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            TransformHelper<FPP,Cpu> th;
            this->tocpu(th);
            auto thn = th.normalize(meth);
            auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
            delete thn;
            return gpu_thn;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::transpose()
        {
			//TODO: factor with CPU function in parent class
            auto t = new TransformHelper<FPP,GPU2>(*this, true, false);
            return t;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::conjugate()
        {
			//TODO: factor with CPU function in parent class
            auto t = new TransformHelper<FPP,GPU2>(*this, false, true);
            return t;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::adjoint()
        {
			//TODO: factor with CPU function in parent class
            auto t = new TransformHelper<FPP,GPU2>(*this, true, true);
            return t;
        }

    template<typename FPP>
        faust_unsigned_int TransformHelper<FPP,GPU2>::getNBytes() const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            faust_unsigned_int nbytes = 0;
            for(auto fac : *this)
            {
                if(dynamic_cast<Faust::MatDense<FPP, GPU2>*>(fac))
                    nbytes += fac->getNbCol() * fac->getNbRow() * sizeof(FPP);
                else if (dynamic_cast<Faust::MatSparse<FPP, GPU2>*>(fac))
                    nbytes += fac->getNonZeros() * (sizeof(FPP) + sizeof(int)) + (fac->getNbRow() + 1) * sizeof(int); // by default storage index is int
                else
                    throw runtime_error("Unknown matrix type.");
            }
            return nbytes;
        }
	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::optimize_multiply(std::function<void()> f, const bool transp /* deft to false */, const bool inplace, /* deft to 1 */ const int nsamples, const char* op_name)
        {
   			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
			std::vector<string> meth_names = {"DEFAULT_L2R", "DYNPROG"}; //TODO: it should be a function of faust_prod_opt module
			std::vector<int> meth_ids = {DEFAULT_L2R, DYNPROG};
			TransformHelper<FPP,GPU2>* t_opt = nullptr;
			int NMETS = 2;
			std::chrono::duration<double> * times = new std::chrono::duration<double>[NMETS]; //use heap because of VS14 (error C3863)
			int old_meth = this->get_mul_order_opt_mode();
			int nmuls = nsamples, opt_meth=0;
#if DEBUG_OPT_MUL
			cout << "nsamples used to measure time: " << nmuls << endl;
#endif
			for(int i=0; i < NMETS; i++)
			{
				this->set_FM_mul_mode(meth_ids[i]);
				auto start = std::chrono::system_clock::now();
				for(int j=0;j < nmuls; j++)
				{
					f();
				}
				auto end = std::chrono::system_clock::now();
				times[i] = end-start;
			}
			for(int i=0; i < NMETS-1; i++)
			{
				opt_meth = times[opt_meth]<times[i+1]?opt_meth:i+1;
			}
			if(inplace)
			{
				this->set_FM_mul_mode(meth_ids[opt_meth]);
				t_opt = this;
			}
			else
			{
				t_opt = new TransformHelper<FPP, GPU2>(this->transform->data, 1.0, false, false, true);
				cout << "best method measured in time on operation "<< op_name << " is: " << meth_names[opt_meth] << endl;
#if DEBUG_OPT_MUL
				cout << "all times: ";
				for(int i = 0; i < NMETS; i ++)
					cout << times[i].count() << " ";
				cout << endl;
#endif
				t_opt->set_FM_mul_mode(meth_ids[opt_meth]);
				// leave the current Faust unchanged
				this->set_FM_mul_mode(old_meth);
			}
			delete [] times;
			t_opt->copy_transconj_state(*this);
			return t_opt;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::fourierFaust(unsigned int n, const bool norma/*=true*/)
        {
            auto cpu_faust = TransformHelper<FPP,Cpu>::fourierFaust(n, norma);
            TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
            delete cpu_faust;
            return gpu_faust;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::pruneout(const int nnz_tres, const int npasses/*=-1*/, const bool only_forward/*=false*/)
        {
            TransformHelper<FPP,Cpu> th;
            this->tocpu(th);
            auto thn = th.pruneout(nnz_tres, npasses, only_forward);
            auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
            delete thn;
            return gpu_thn;
        }


    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::clone(int32_t dev_id/*=-1*/, void* stream/*=nullptr*/)
        {
            //TODO: take the requested device into account
            return TransformHelperGen<FPP,GPU2>::clone();
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::randFaust(int faust_nrows, int faust_ncols, RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row)
        {
            auto cpu_faust = TransformHelper<FPP,Cpu>::randFaust(faust_nrows, faust_ncols, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
            //	TransformHelper<FPP,GPU2>::TransformHelper(const std::vector<MatGeneric<FPP,GPU2> *>& facts, const FPP lambda_/*= (FPP)1.0*/, const bool optimizedCopy/*=false*/, const bool cloning_fact /*= true*/, const bool internal_call/*=false*/)
            TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust/*TODO: dev_id and stream ?*/);
            //		void TransformHelper<FPP,GPU2>::push_back(const MatGeneric<FPP,Cpu>* M, const bool optimizedCopy/*=false*/, const int32_t dev_id/*=-1*/, const void* stream/*=nullptr*/)
            //			for(auto cpu_fact: *cpu_faust)
            //				gpu_faust->push_back(cpu_fact, false, -1/*TODO: replace dev_id by an arg passed to the func */, nullptr);

            delete cpu_faust;
            return gpu_faust;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::randFaust(RandFaustType t, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int min_dim_size, unsigned int max_dim_size, float density, bool per_row)
        {
            return TransformHelper<FPP,GPU2>::randFaust(-1, -1, t, min_num_factors, max_num_factors, min_dim_size, max_dim_size, density, per_row);
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::hadamardFaust(unsigned int n, const bool norma/*=true*/)
        {
            //			std::cerr << "Warning: GPU2 hadamardFaust is implemented by copying the Faust on CPU RAM and copying them back." << std::endl;
            auto cpu_faust = TransformHelper<FPP,Cpu>::hadamardFaust(n, norma);
            TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
            delete cpu_faust;
            return gpu_faust;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::eyeFaust(unsigned int n, unsigned int m)
        {
            //			std::cerr << "Warning: GPU2 eyeFaust is implemented by copying the Faust on CPU RAM and copying them back." << std::endl;
            auto cpu_faust = TransformHelper<FPP,Cpu>::eyeFaust(n, m);
            TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
            delete cpu_faust;
            return gpu_faust;
        }


    template<typename FPP>
	TransformHelper<FPP, GPU2>* TransformHelper<FPP,GPU2>::randBSRFaust(unsigned int faust_nrows, unsigned int faust_ncols, unsigned int min_num_factors, unsigned int max_num_factors, unsigned int bnrows, unsigned int bncols, float density/*=.1f*/)
	{
		auto cpu_faust = TransformHelper<FPP,Cpu>::randBSRFaust(faust_nrows, faust_ncols, min_num_factors, max_num_factors, bnrows, bncols, density);
		TransformHelper<FPP,GPU2>* gpu_faust = new TransformHelper<FPP,GPU2>(*cpu_faust, -1, nullptr /*TODO: dev_id and stream ?*/);
		delete cpu_faust;
		return gpu_faust;
	}

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::swap_rows(const faust_unsigned_int id1,
                const faust_unsigned_int id2,
                const bool permutation/*=false*/,
                const bool inplace/*=false*/,
                const bool check_transpose/*=true*/)
        {
            TransformHelper<FPP,Cpu> th;
            this->tocpu(th);
            auto thn = th.swap_rows(id1, id2, permutation, inplace, check_transpose);
            auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
            delete thn;
            return gpu_thn;
        }

    template<typename FPP>
        TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::swap_cols(const faust_unsigned_int id1,
                const faust_unsigned_int id2,
                const bool permutation/*=false*/,
                const bool inplace/*=false*/,
                const bool check_transpose/*=true*/)
        {
            TransformHelper<FPP,Cpu> th;
            this->tocpu(th);
            auto thn = th.swap_cols(id1, id2, permutation, inplace, check_transpose);
            auto gpu_thn = new TransformHelper<FPP,GPU2>(*thn, -1, nullptr);
            delete thn;
            return gpu_thn;
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
                int* rowptr,
                int* col_ids,
                FPP* elts,
                faust_unsigned_int* nnz,
                faust_unsigned_int* num_rows,
                faust_unsigned_int* num_cols,
                const bool transpose /* = false*/) const
        {
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_sliced_Transform();
			const_cast<Faust::TransformHelper<FPP, GPU2>*>(this)->eval_fancy_idx_Transform();
            this->transform->get_fact(this->is_transposed?this->size()-id-1:id, rowptr, col_ids, elts, nnz, num_rows, num_cols, this->is_transposed ^ transpose);
            if(this->is_conjugate)
                Faust::conjugate(elts, *nnz);
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
                FPP* elts,
                faust_unsigned_int* num_rows,
                faust_unsigned_int* num_cols,
                const bool transpose /* = false*/) const
        {
            TransformHelperGen<FPP,GPU2>::get_fact(id, elts, num_rows, num_cols, transpose);
        }


    template<typename FPP>
        template<typename FPP2>
        TransformHelper<Real<FPP2>,GPU2>* TransformHelper<FPP,GPU2>::real()
        {
			this->eval_sliced_Transform();
			this->eval_fancy_idx_Transform();
            std::vector<MatGeneric<Real<FPP2>,GPU2>*> real_data;
            MatSparse<FPP, GPU2> *curfac_sp;
            MatDense<FPP, GPU2> *curfac_ds;
            for(auto curfac: this->transform->data)
            {
                if(curfac_ds = dynamic_cast<MatDense<FPP, GPU2>*>(curfac))
                {
                    auto real_fac = new MatDense<Real<FPP2>,GPU2>(curfac->getNbRow(), curfac->getNbCol(), /* cpu_data*/ nullptr, /* no_alloc*/ false, /* dev_id*/ -1, /* stream*/ nullptr);
                    *real_fac = curfac_ds->template to_real<Real<FPP2>>();
                    real_data.push_back(real_fac);
                }
                else if(curfac_sp = dynamic_cast<MatSparse<FPP, GPU2>*>(curfac))
                {
                    auto real_fac = new MatSparse<Real<FPP2>,GPU2>(curfac->getNbRow(), curfac->getNbCol(), /* nnz*/ 0, /* values */ nullptr, /* rowptr*/ nullptr, /* colinds*/ nullptr);
                    *real_fac = curfac_sp->template to_real<Real<FPP2>>();
                    real_data.push_back(real_fac);
                }
                else
                {
                    throw std::runtime_error("real() failed because a factor is neither a MatDense nor a MatSparse");
                }
            }
            return new TransformHelper<Real<FPP2>, GPU2>(real_data, 1.0, false, false, true);

        }

    template<typename FPP>
        FPP TransformHelper<FPP,GPU2>::get_item(faust_unsigned_int i, faust_unsigned_int j)
        {
            MatDense<FPP, GPU2> M;
            faust_unsigned_int out_id;
            TransformHelperGen<FPP,GPU2>::get_item(i, j, M, out_id);
            return M.tocpu().getData()[out_id];
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_fact_bsr_info(const faust_unsigned_int id,
                size_t& bdata_sz,
                size_t& browptr_sz,
                size_t& bcolinds_sz, 
                size_t& bnnz,
                size_t& bnrows,
                size_t& bncols) const
        {
            throw std::runtime_error("TransformHelper<FPP,GPU2>::get_fact_bsr_info error: GPU2 doesn't support the BSR matrix yet.");
        }

    template<typename FPP>
        void TransformHelper<FPP,GPU2>::get_fact(const faust_unsigned_int id,
                FPP* bdata,
                int* brow_ptr,
                int* bcol_inds) const
        {

            throw std::runtime_error("TransformHelper<FPP,GPU2>::get_fact error: GPU2 doesn't support the BSR matrix yet.");
        }

	template<typename FPP>
		void TransformHelper<FPP,GPU2>::init_fancy_idx_transform(TransformHelper<FPP,GPU2>* th, faust_unsigned_int* row_ids, faust_unsigned_int num_rows, faust_unsigned_int* col_ids, faust_unsigned_int num_cols)
		{
			TransformHelperGen<FPP, GPU2>::init_fancy_idx_transform(th, row_ids, num_rows, col_ids, num_cols);
		}

	template<typename FPP>
		FPP* Faust::TransformHelper<FPP,GPU2>::sliceMultiply(const Slice s[2], const FPP* cpu_X, FPP* cpu_out/*=nullptr*/, int X_ncols/*=1*/) const
		{
			//TODO: take care of eval_sliced_Transform indirect calls
            int32_t X_nrows;
			if(s[1].end_id != this->getNbCol() || s[1].start_id != 0)
				X_nrows = s[1].end_id-s[1].start_id;
			else
				X_nrows = this->getNbCol();
            MatDense<FPP,GPU2> gpu_X(X_nrows, X_ncols, cpu_X, false);
            MatDense<FPP,GPU2> gpu_M = this->transform->sliceMultiply(s, gpu_X, this->isTransposed2char());
			if(cpu_out == nullptr)
			{
				auto is_row_sliced = s[0].end_id != this->getNbRow() || s[0].start_id != 0;
				auto out_nrows = is_row_sliced?s[0].end_id-s[0].start_id:this->getNbRow();
				auto out_ncols = X_ncols;
				cpu_out = new FPP[out_nrows*out_ncols*sizeof(FPP)];
			}
            gpu_M.tocpu(cpu_out, nullptr);
			return cpu_out;
		}

	template<typename FPP>
		FPP* Faust::TransformHelper<FPP,GPU2>::indexMultiply(faust_unsigned_int* ids[2], size_t id_lens[2], const FPP* cpu_X, int X_ncols/*=1*/, FPP* cpu_out/*=nullptr*/) const
		{
			int32_t X_nrows;
			if(id_lens[1] > 0)
				X_nrows = id_lens[1];
			else
				X_nrows = this->getNbCol();
			MatDense<FPP,GPU2> gpu_X(X_nrows, X_ncols, cpu_X, false);
			MatDense<FPP,GPU2> gpu_M = this->transform->indexMultiply(ids, id_lens, gpu_X, this->isTransposed2char());
			if(cpu_out == nullptr)
			{
				auto out_nrows = id_lens[0]>0?id_lens[0]:this->getNbRow();
				auto out_ncols = X_ncols;
				cpu_out = new FPP[out_nrows*out_ncols*sizeof(FPP)];
			}
			gpu_M.tocpu(cpu_out, nullptr);
			return cpu_out;
		}


	template<typename FPP>
		TransformHelper<FPP, GPU2>* TransformHelper<FPP,GPU2>::optButterflyFaust(const TransformHelper<FPP, GPU2>* F)
		{
			return TransformHelperButterfly<FPP, GPU2>::optFaust(F);
		}

template<typename FPP>
	TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::fourierFaustOpt(unsigned int n, const bool norma/*=true*/)
	{
		// uncomment this code when TransformHelperButterfly will be replaced (as on CPU) and optButterflyFaust/optFaust implemented on GPU
//		auto F = TransformHelper<FPP,GPU2>::fourierFaust(n, norma);
//		auto Fo = TransformHelper<FPP,GPU2>::optButterflyFaust(F);
//		delete F;
//		return Fo;
		return TransformHelperButterfly<FPP, GPU2>::fourierFaust(n, norma);
	}

}

#include "faust_TransformHelper_cat_gpu.hpp"
