namespace Faust
{

	template<typename FPP, FDevice D>
		void faust_MatSparse_to_torch_Tensor(const Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev, const bool clone)
		{
			torch::Tensor values = torch::from_blob(const_cast<FPP*>(spm.getValuePtr()), {spm.getNonZeros()}, torch::TensorOptions().dtype(torch::kFloat64).device(dev));//.clone();
			//		cout << "tensor values:" << values << endl;
			torch::Tensor col = torch::from_blob(const_cast<int*>(spm.getInnerIndexPtr()), {spm.getNonZeros()}, torch::TensorOptions().dtype(torch::kI32).device(dev));//.clone();
			//		cout << "tensor col:" << col << endl;
			faust_unsigned_int* rows = new faust_unsigned_int[spm.getNonZeros()];
			for(int i = 0; i < spm.getNbRow(); i++)
			{
				for(int j=spm.getOuterIndexPtr()[i];j < spm.getOuterIndexPtr()[i+1]; j++)
				{
					rows[j] = i;
					//				cout << " " << rows[j];
				}
			}
			torch::Tensor row = torch::from_blob(rows, {spm.getNonZeros()}, torch::TensorOptions().dtype(torch::kI64).device(dev)).clone(); // clone is mandatory for rows
			if(clone)
			{
//				row = row.clone();
//				col = col.clone();
				values = values.clone();
			}
//			row = row.to(torch::kI64);
			col = col.to(torch::kI64); // mandatory conversion because torch forces to use same size types for indices and values (even if indices are integers and values floats)
			//		cout << "tensor row:" << row << endl;
			delete [] rows;
			torch::Tensor indices = at::stack({row, col}, /* dim */ 0);
			t = torch::sparse_coo_tensor(indices, values);
			//		cout << "tensor size: " << t.size(0) << " x " << t.size(1) << " t is sparse:" << t.is_sparse() << endl;
		}

	template<typename FPP, FDevice D>
		void faust_MatDense_to_torch_Tensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev, const bool clone)
		{
			t = torch::from_blob(const_cast<FPP*>(dm.getData()), {dm.getNbCol(), dm.getNbRow()},torch::TensorOptions().dtype(torch::kFloat64).device(dev));//.clone();
			t = t.t();
			if(clone)
				t = t.clone();
		}

	template<typename FPP, FDevice D>
		void torch_Tensor_to_faust_MatDense(const torch::Tensor & t, Faust::MatDense<FPP,D> & dm)
		{
			dm = Faust::MatDense<FPP,Cpu>(t.data_ptr<FPP>(), t.size(1), t.size(0));
			dm.transpose();
		}

	template<typename FPP, FDevice D>
		void faust_matvec_to_torch_TensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev, const bool clone)
		{
			const Faust::MatSparse<FPP,D> *spm;
			const Faust::MatDense<FPP,D> *dm;
//			torch::Tensor t;
			tl.resize(ml.size());
			int i = 0;
			for(auto m : ml)
			{
				if(spm = dynamic_cast<Faust::MatSparse<FPP,D>*>(m))
					faust_MatSparse_to_torch_Tensor(*spm, tl[i++], dev, clone);
//					faust_MatSparse_to_torch_Tensor(*spm, t, dev, clone);
				else if(dm = dynamic_cast<Faust::MatDense<FPP,D>*>(m))
					faust_MatDense_to_torch_Tensor(*dm, tl[i++], dev, clone);
//					faust_MatDense_to_torch_Tensor(*dm, t, dev, clone);
//				tl.push_back(t);
			}
		}

	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& ml, const torch::Tensor* op, at::DeviceType dev, const bool chain_opt, const bool contiguous_dense_to_torch)
	{
		bool all_dense = true;
		std::vector<torch::Tensor> mlc;
		for(auto t: ml)
		{
			all_dense &= !t.is_sparse();
			mlc.push_back(t);
		}
		if(op)
		{
			all_dense &= !op->is_sparse();
			mlc.push_back(*op);
		}
		if(all_dense)
			return torch::chain_matmul(mlc); //chain_opt is useless because I suppose torch does its own chain opt.
		if(chain_opt)
			return std::move(tensor_chain_mul_opt(mlc, nullptr, dev));
		auto it = mlc.end()-1;
		auto res = *(it);
		if(res.is_sparse())
			res = res.to_dense();
		std::vector<torch::Tensor> dense_contiguous_facts;
		while(it != mlc.begin())
		{
			auto f = *(--it);
			if(f.is_sparse())
			{
				if(dense_contiguous_facts.size() > 0)
				{
					//multiply chain of dense tensors before
					dense_contiguous_facts.push_back(res);
					res = torch::chain_matmul(dense_contiguous_facts);
					dense_contiguous_facts.erase(dense_contiguous_facts.begin(), dense_contiguous_facts.end());
				}
				res = at::_sparse_mm(f, res);
			}
			else
				if(contiguous_dense_to_torch)
					dense_contiguous_facts.insert(dense_contiguous_facts.begin(), f);
				else
					res = torch::matmul(f, res);
		}
		assert(res.size(0) == ml.size(0));
		assert(op == nullptr && res.size(1) == ml.size(1) || res.size(1) == op->size(1));
		return std::move(res); //explicit move but should work auto because Tensor class overrides move operator= and ctor
	}

	torch::Tensor tensor_chain_mul_opt(const std::vector<torch::Tensor>& ml, const torch::Tensor* op, at::DeviceType dev)
	{
		// cost to apply a on b
		auto cost = [](const torch::Tensor &a, const torch::Tensor &b)
		{
			uint64_t a_cost = a.size(0)*a.size(1);
			uint64_t b_cost;
			if(b.is_sparse())
				b_cost = b._nnz();	
			else
				b_cost = b.size(1);
			return a_cost*b_cost;
		};
		// argmin of costs v
		auto argmin = [](std::vector<uint64_t> v)
		{
			int i=0 /* argmin*/, j=0 /* element index*/;
			uint64_t min = numeric_limits<uint64_t>::max();
			for(auto e: v)
			{
				if(e < min)
				{
					min = e;
					i = j;
				}
				j++;
			}
			return i;
		};
		std::vector<const torch::Tensor*> mlc;
		for(int i=0;i<ml.size();i++) mlc.push_back(&ml[i]);
		if(op != nullptr)
			mlc.push_back(op);
		std::vector<uint64_t> costs(mlc.size()-1);
		for(int i=0;i<costs.size();i++)
			costs[i] = cost(*mlc[i], *mlc[i+1]);
		torch::Tensor* res = new torch::Tensor();
		vector<torch::Tensor*> res_list;
		res_list.push_back(res);
		while(mlc.size() > 1)
		{
//			for (int i = 0; i < mlc.size(); i++) cout << mlc[i] << "[" << mlc[i]->size(0) << "x" << mlc[i]->size(1) << "] "; cout << endl;
			auto i = argmin(costs);
			auto f1 = mlc[i];
			auto f2 = mlc[i+1];
			if(res != f1 && res != f2)
			{
				res = new torch::Tensor();
				res_list.push_back(res);
			}
//			cout << "argmin i:" << i << endl;
			if(f2->is_sparse())
				if(f1->is_sparse())
					*res = at::_sparse_mm(*f1, f2->to_dense());
				else
					*res = torch::t(torch::matmul(torch::t(*f2), torch::t(*f1)));
			else
				if(f1->is_sparse())
					*res = at::_sparse_mm(*f1, *f2);
				else
					*res = torch::matmul(*f1, *f2);
			mlc.insert(mlc.begin()+i, res);
			mlc.erase(mlc.begin()+i+1,mlc.begin()+i+3);
			if(mlc.size() > i+1)
			{
				costs.insert(costs.begin()+i, cost(*res, *mlc[i+1]));
				costs.erase(costs.begin()+i+1, costs.begin()+i+3);
			}
			else
				costs.erase(costs.end()-1);
		}
		assert(mlc.size() == 1);
		torch::Tensor t = *mlc[0];
		for(auto res: res_list)
			delete res;
		return std::move(t);
	}

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch)
		{
			std::vector<torch::Tensor> tl;
			faust_matvec_to_torch_TensorList(ml, tl, on_gpu?at::kCUDA:at::kCPU, clone);
			tensor_chain_mul(tl, out, op, on_gpu, clone, chain_opt, contiguous_dense_to_torch);
		}

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<torch::Tensor>& tl, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch)
	{
		torch::Tensor top, tres;
		const Faust::MatSparse<FPP,D> *spm;
		const Faust::MatDense<FPP,D> *dm;
		if(op)
		{
			if(spm = dynamic_cast<const Faust::MatSparse<FPP,D>*>(op))
				faust_MatSparse_to_torch_Tensor(*spm, top, on_gpu?at::kCUDA:at::kCPU, clone);
			else if(dm = dynamic_cast<const Faust::MatDense<FPP,D>*>(op))
				faust_MatDense_to_torch_Tensor(*dm, top, on_gpu?at::kCUDA:at::kCPU, clone);
			tres = tensor_chain_mul(tl, &top, on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch);
		}
		else
			tres = tensor_chain_mul(tl, static_cast<torch::Tensor*>(nullptr), on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch);
		out = Faust::MatDense<FPP,Cpu>(tres.data_ptr<FPP>(), tres.size(1), tres.size(0));
		out.transpose();
	}

}
