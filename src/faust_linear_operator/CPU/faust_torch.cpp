#include "faust_torch.h"

void display_TensorList(std::vector<torch::Tensor>& tl, const bool transpose /*= true*/)
{
	if(transpose)
		for (int i = tl.size()-1; i >= 0; i--) {cout << "Tensor: " << tl.size()-1-i << " [" << tl[i].size(1) << "x" << tl[i].size(0) << "] " << (tl[i].is_sparse()?"SPARSE":"DENSE"); cout << endl;}
	else
		for (int i = 0; i < tl.size(); i++) {cout << "Tensor: " << i << " [" << tl[i].size(0) << "x" << tl[i].size(1) << "] " << (tl[i].is_sparse()?"SPARSE":"DENSE"); cout << endl;};
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
		assert((op != nullptr && op_on_left && res.size(0) == op->size(0)) || ((! op_on_left || op == nullptr) && res.size(0) == tl[0].size(0)) || op == nullptr);
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
