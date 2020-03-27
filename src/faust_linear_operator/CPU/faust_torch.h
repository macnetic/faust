#ifndef __FAUST_TORCH__
#define __FAUST_TORCH__
#include "faust_MatSparse.h"
#include <vector>
#include <torch/torch.h>

namespace Faust
{
	/**
	 * Converts a Faust::MatDense to a torch::Tensor.
	 *
	 * Unfortunately this is not a view but a copy which is not a problem when dev == at::kCUDA.
	 */
	template<typename FPP, FDevice D>
		void faust_MatDense_to_torch_Tensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev = at::kCPU);


	/**
	 * Converts a Faust::MatSparse to a torch::Tensor.
	 *
	 * Unfortunately this is not a view but a copy which is not a problem when dev == at::kCUDA.
	 */
	template<typename FPP, FDevice D>
		void faust_MatSparse_to_torch_Tensor(Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev = at::kCPU);

	/**
	 * Converts a Faust::MatGeneric vector to a torch::TensorList (vector alias).
	 *
	 * Unfortunately this is not a view but a copy which is not a problem when dev == at::kCUDA.
	 */
	template<typename FPP, FDevice D>
		void faust_matvec_to_torch_TensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev = at::kCPU);

	/**
	 * Computes the tensor chain product of ml and applies it optionally to the tensor op.
	 *
	 * Returns the result as a Tensor.
	 */
	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& ml, const torch::Tensor* op= nullptr, at::DeviceType dev = at::kCPU);

	/**
	 * Computes the matrix chain product of ml and applies it optionally to the matrix op.
	 *
	 * This function converts all the matrices to Tensors before and then computes the tensor product.
	 *
	 * Returns the result as a Faust::MatDense.
	 */
	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op = nullptr, const bool on_gpu = false);
}
#include "faust_torch.hpp"
#endif
