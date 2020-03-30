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
			t = torch::from_blob(const_cast<FPP*>(dm.getData()), {dm.getNbRow(), dm.getNbCol()},torch::TensorOptions().dtype(torch::kFloat64).device(dev));//.clone();
			if(clone)
				t = t.clone();
		}

	template<typename FPP, FDevice D>
		void faust_matvec_to_torch_TensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev, const bool clone)
		{
			const Faust::MatSparse<FPP,D> *spm;
			const Faust::MatDense<FPP,D> *dm;
			torch::Tensor t;
			for(auto m : ml)
			{
				if(spm = dynamic_cast<Faust::MatSparse<FPP,D>*>(m))
					faust_MatSparse_to_torch_Tensor(*spm, t, dev, clone);
				else if(dm = dynamic_cast<Faust::MatDense<FPP,D>*>(m))
					faust_MatDense_to_torch_Tensor(*dm, t, dev, clone);
				tl.push_back(t);
			}
		}

	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& ml, const torch::Tensor* op, at::DeviceType dev)
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
			return torch::chain_matmul(mlc);
		auto it = mlc.end()-1;
		auto res = *(it);
		if(res.is_sparse())
			res = res.to_dense();
		while(it != mlc.begin())
		{
			auto f = *(--it);
			torch::Tensor t;
			if(f.is_sparse())
				t = at::_sparse_mm(f, res);
			else
				t = torch::matmul(f, res);
			res = t;
		}
		return res;
	}

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone)
		{
			std::vector<torch::Tensor> tl;
			faust_matvec_to_torch_TensorList(ml, tl);
			torch::Tensor t = tensor_chain_mul(tl), top, tres;
			const Faust::MatSparse<FPP,D> *spm;
			const Faust::MatDense<FPP,D> *dm;
			if(op)
			{
				if(spm = dynamic_cast<const Faust::MatSparse<FPP,D>*>(op))
					faust_MatSparse_to_torch_Tensor(*spm, top, on_gpu?at::kCUDA:at::kCPU, clone);
				else if(dm = dynamic_cast<const Faust::MatDense<FPP,D>*>(op))
					faust_MatDense_to_torch_Tensor(*dm, top, on_gpu?at::kCUDA:at::kCPU, clone);
				tres = tensor_chain_mul(tl, &top, on_gpu?at::kCUDA:at::kCPU);
			}
			else
				tres = tensor_chain_mul(tl, static_cast<torch::Tensor*>(nullptr), on_gpu?at::kCUDA:at::kCPU);
			out = Faust::MatDense<FPP,Cpu>(tres.data_ptr<FPP>(), tres.size(0), tres.size(1));
		}
}
