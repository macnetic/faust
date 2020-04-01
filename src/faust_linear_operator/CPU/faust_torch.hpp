namespace Faust
{

	template<typename FPP, FDevice D>
		void convMatSparseToTensor(const Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev, const bool clone, const bool transpose /* = true*/)
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
			torch::Tensor indices;
			if(transpose)
			{
				//reverse row and col to take the matrix as a transpose mat
				indices = at::stack({col, row}, /* dim */ 0);
				t = torch::sparse_coo_tensor(indices, values);
				t.sparse_resize_({spm.getNbCol(), spm.getNbRow()}, t.sparse_dim(), t.dense_dim());
			}
			else
			{
				indices = at::stack({row, col}, /* dim */ 0);
				t = torch::sparse_coo_tensor(indices, values);
				t.sparse_resize_({spm.getNbRow(), spm.getNbCol()}, t.sparse_dim(), t.dense_dim());
			}
			//		cout << "tensor size: " << t.size(0) << " x " << t.size(1) << " t is sparse:" << t.is_sparse() << endl;
			assert(t._nnz() == spm.getNonZeros());
			assert(t.size(0) == spm.getNbCol() && t.size(1) == spm.getNbRow());
		}

	template<typename FPP, FDevice D>
		void convMatDenseToTensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev, const bool clone, const bool transpose /* = true*/)
		{
			uint64_t nrows, ncols;

			// number of nrows and ncols are inverted because the data is taken as a transpose matrix (Faust::MatDense is column-major order while torch is row-major order)
			// it saves the need/transpose
			nrows = dm.getNbCol();
			ncols = dm.getNbRow();
			t = torch::from_blob(const_cast<FPP*>(dm.getData()), {nrows, ncols},torch::TensorOptions().dtype(torch::kFloat64).device(dev));//.clone();
			if(! transpose)
				// need to transpose when transpose is false! conversion to torch row-major order
				// while Faust is in col-major order (the conversion is equivalent to a transpose)
				t = t.t();
			if(clone && transpose)
				t = t.clone();
			// if clone == true && transpose == false // the transpose above already cloned the data
		}

	template<typename FPP, FDevice D>
		void convTensorToMatDense(const torch::Tensor & t, Faust::MatDense<FPP,D> & dm, const bool transpose /* = true*/)
		{
			if(transpose)
			{
				dm = Faust::MatDense<FPP,Cpu>(t.data_ptr<FPP>(), t.size(1), t.size(0));
			}
			else
			{
				dm = Faust::MatDense<FPP,Cpu>(t.data_ptr<FPP>(), t.size(1), t.size(0));
				dm.transpose();
				// need to transpose when transpose is false! conversion from torch row-major order
				// while Faust is in col-major order (the conversion is equivalent to a transpose)
			}
		}

	template<typename FPP, FDevice D>
		void convMatGenListToTensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev, const bool clone /* = false*/, const bool transpose /* = true*/)
		{
			const Faust::MatSparse<FPP,D> *spm;
			const Faust::MatDense<FPP,D> *dm;
			tl.resize(ml.size());
			int i;
			if(transpose)
			{
				i = tl.size()-1; // transpose order == reverse order
				for(auto m : ml)
				{
					if(spm = dynamic_cast<Faust::MatSparse<FPP,D>*>(m))
						convMatSparseToTensor(*spm, tl[i--], dev, clone);
					else if(dm = dynamic_cast<Faust::MatDense<FPP,D>*>(m))
						convMatDenseToTensor(*dm, tl[i--], dev, clone);
				}
			}
			else
			{
				i = 0;
				for(auto m : ml)
				{
					if(spm = dynamic_cast<Faust::MatSparse<FPP,D>*>(m))
						convMatSparseToTensor(*spm, tl[i++], dev, clone);
					else if(dm = dynamic_cast<Faust::MatDense<FPP,D>*>(m))
						convMatDenseToTensor(*dm, tl[i++], dev, clone);
				}

			}
		}

	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& tl, const torch::Tensor* op, at::DeviceType dev, const bool chain_opt, const bool contiguous_dense_to_torch, const bool op_on_left /*=true*/)
	{
		bool all_dense = true;
		std::vector<torch::Tensor> tlc;
		for(auto t: tl)
		{
			all_dense &= !t.is_sparse();
			tlc.push_back(t);
		}
		if(op)
		{
			all_dense &= !op->is_sparse();
			if(op_on_left)
				tlc.insert(tlc.begin(), *op);
			else
				tlc.push_back(*op);
		}
		if(all_dense)
			return torch::chain_matmul(tlc); //chain_opt is useless because I suppose torch does its own chain opt.
		if(chain_opt)
			return std::move(tensor_chain_mul_opt(tlc, nullptr, dev));
		auto it = tlc.end()-1;
		auto res = *(it);
		if(res.is_sparse())
			res = res.to_dense();
		std::vector<torch::Tensor> dense_contiguous_facts;
		while(it != tlc.begin())
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
		if(contiguous_dense_to_torch && dense_contiguous_facts.size() > 0)
		{
			//multiply chain of dense tensors at the end/start of tlc
			dense_contiguous_facts.push_back(res);
			res = torch::chain_matmul(dense_contiguous_facts);
			dense_contiguous_facts.erase(dense_contiguous_facts.begin(), dense_contiguous_facts.end());
		}
		// don't worry assert is enabled only in debug mode (when DEBUG is defined)
		assert((op != nullptr && op_on_left && res.size(0) == op.size(0)) || ((! op_on_left || op == nullptr) && res.size(0) == tl[0].size(0)) || op == nullptr);
		assert(((op == nullptr || op_on_left) && res.size(1) == (*(tl.end()-1)).size(1)) || op != nullptr && res.size(1) == op->size(1));
		return std::move(res); //explicit move but should work auto because Tensor class overrides move operator= and ctor
	}

	torch::Tensor tensor_chain_mul_opt(const std::vector<torch::Tensor>& ml, const torch::Tensor* op, at::DeviceType dev, const bool op_on_left /* = true */)
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
			if(op_on_left)
				mlc.insert(mlc.begin(), op);
			else
				mlc.push_back(op);
		std::vector<uint64_t> costs(mlc.size()-1);
		for(int i=0;i<costs.size();i++)
			costs[i] = cost(*mlc[i], *mlc[i+1]);
		torch::Tensor* res = new torch::Tensor();
		vector<torch::Tensor*> res_list;
		res_list.push_back(res);
		while(mlc.size() > 1)
		{
			auto i = argmin(costs);
			auto f1 = mlc[i];
			auto f2 = mlc[i+1];
			if(res != f1 && res != f2)
			{
				res = new torch::Tensor();
				res_list.push_back(res);
			}
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
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch, const bool transpose /* = true */)
		{
			std::vector<torch::Tensor> tl;
			convMatGenListToTensorList(ml, tl, on_gpu?at::kCUDA:at::kCPU, clone, transpose);
			tensor_chain_mul(tl, out, op, on_gpu, clone, chain_opt, contiguous_dense_to_torch, transpose);
		}

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<torch::Tensor>& tl, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch, const bool transpose /* = true */)
	{
		torch::Tensor top, tres;
		const Faust::MatSparse<FPP,D> *spm;
		const Faust::MatDense<FPP,D> *dm;
		if(op)
		{
			if(spm = dynamic_cast<const Faust::MatSparse<FPP,D>*>(op))
				convMatSparseToTensor(*spm, top, on_gpu?at::kCUDA:at::kCPU, clone, transpose);
			else if(dm = dynamic_cast<const Faust::MatDense<FPP,D>*>(op))
				convMatDenseToTensor(*dm, top, on_gpu?at::kCUDA:at::kCPU, clone, transpose);
			tres = tensor_chain_mul(tl, &top, on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch, transpose /* op_on_left if transpose */);
		}
		else
			tres = tensor_chain_mul(tl, static_cast<torch::Tensor*>(nullptr), on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch, transpose /* op_on_left if transpose */);
		convTensorToMatDense(tres, out, transpose);
	}

	void display_TensorList(std::vector<torch::Tensor>& tl, const bool transpose /*= true*/)
	{
		if(transpose)
			for (int i = tl.size()-1; i >= 0; i--) {cout << "Tensor: " << tl.size()-1-i << " [" << tl[i].size(1) << "x" << tl[i].size(0) << "] " << (tl[i].is_sparse()?"SPARSE":"DENSE"); cout << endl;}
		else
			for (int i = 0; i < tl.size(); i++) {cout << "Tensor: " << i << " [" << tl[i].size(0) << "x" << tl[i].size(1) << "] " << (tl[i].is_sparse()?"SPARSE":"DENSE"); cout << endl;};
	}

}
