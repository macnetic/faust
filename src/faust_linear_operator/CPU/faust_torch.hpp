namespace Faust
{

	template<typename FPP, FDevice D>
		void convMatSparseToTensor(const Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev, const bool clone, const bool transpose /* = true*/)
		{
			long long spm_nnz = static_cast<long long>(spm.getNonZeros());
			long int spm_nrows = spm.getNbRow(), spm_ncols = spm.getNbCol();
			torch::Tensor values = torch::from_blob(const_cast<FPP*>(spm.getValuePtr()), {spm_nnz}, torch::TensorOptions().dtype(torch::kFloat64).device(dev));//.clone();
			//		cout << "tensor values:" << values << endl;
			torch::Tensor col = torch::from_blob(const_cast<int*>(spm.getInnerIndexPtr()), {spm_nnz}, torch::TensorOptions().dtype(torch::kI32).device(dev));//.clone();
			//		cout << "tensor col:" << col << endl;
			faust_unsigned_int* rows = new faust_unsigned_int[spm_nnz];
			for(int i = 0; i < spm_nrows; i++)
			{
				for(int j=spm.getOuterIndexPtr()[i];j < spm.getOuterIndexPtr()[i+1]; j++)
				{
					rows[j] = i;
					//				cout << " " << rows[j];
				}
			}
			torch::Tensor row = torch::from_blob(rows, {spm_nnz}, torch::TensorOptions().dtype(torch::kI64).device(dev)).clone(); // clone is mandatory for rows
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
				t.sparse_resize_({spm_ncols, spm_nrows}, t.sparse_dim(), t.dense_dim()); // this fixes low level bug (where the number of tensor columns is incorrect by one)
			}
			else
			{
				indices = at::stack({row, col}, /* dim */ 0);
				t = torch::sparse_coo_tensor(indices, values);
				t.sparse_resize_({spm_nrows, spm_ncols}, t.sparse_dim(), t.dense_dim());
			}
			//		cout << "tensor size: " << t.size(0) << " x " << t.size(1) << " t is sparse:" << t.is_sparse() << endl;
			assert(t._nnz() == spm.getNonZeros());
			assert(t.size(0) == spm_ncols && t.size(1) == spm_nrows);
		}

	template<typename FPP, FDevice D>
		void convMatDenseToTensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev, const bool clone, const bool transpose /* = true*/)
		{
			long int nrows, ncols;

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
				dm = Faust::MatDense<FPP,D>(t.data_ptr<FPP>(), t.size(1), t.size(0));
			}
			else
			{
				dm = Faust::MatDense<FPP,D>(t.data_ptr<FPP>(), t.size(1), t.size(0));
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

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,D> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch, const bool transpose /* = true */)
		{
			std::vector<torch::Tensor> tl;
			convMatGenListToTensorList(ml, tl, on_gpu?at::kCUDA:at::kCPU, clone, transpose);
			tensor_chain_mul(tl, out, op, on_gpu, clone, chain_opt, contiguous_dense_to_torch, transpose);
		}

	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<torch::Tensor>& tl, Faust::MatDense<FPP,D> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch, const bool transpose /* = true */)
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
			else
				throw std::runtime_error("Only MatSparse and MatDense conversions to tensor are handled.");
			tres = tensor_chain_mul(tl, &top, on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch, transpose /* op_on_left if transpose */);
		}
		else
			tres = tensor_chain_mul(tl, static_cast<torch::Tensor*>(nullptr), on_gpu?at::kCUDA:at::kCPU, chain_opt, contiguous_dense_to_torch, transpose /* op_on_left if transpose */);
		convTensorToMatDense(tres, out, transpose);
	}


}
