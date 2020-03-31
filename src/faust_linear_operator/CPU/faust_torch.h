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
	 *\note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 * \param dev the device at::kCPU or at::kCUDA.
	 * \param clone true to copy the Faust matrix data to create the tensor, false to use the same data without copying (false by default).
	 */
	template<typename FPP, FDevice D>
		void faust_MatDense_to_torch_Tensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev = at::kCPU, const bool clone = false);

	template<typename FPP, FDevice D>
		void torch_Tensor_to_faust_MatDense(const torch::Tensor & t, Faust::MatDense<FPP,D> & dm);


	/**
	 * Converts a Faust::MatSparse to a torch::Tensor.
	 *
	 * \note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 *\param dev the device at::kCPU or at::kCUDA.
	 *\param clone true to copy the Faust matrix data to create the tensor, false to use the same data without copying (false by default).
	 */
	template<typename FPP, FDevice D>
		void faust_MatSparse_to_torch_Tensor(Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev = at::kCPU,  const bool clone = false);

	/**
	 * Converts a Faust::MatGeneric vector to a torch::TensorList (vector alias).
	 *
	 *\param dev the device at::kCPU or at::kCUDA.
	 * \param clone true to copy the Faust matrices data to create the tensors, false to use the same data without copying (false by default).
	 */
	template<typename FPP, FDevice D>
		void faust_matvec_to_torch_TensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev = at::kCPU, const bool clone = false);

	/**
	 * Computes the tensor chain product of ml and applies it optionally to the tensor op.
	 *
	 *\param dev the device at::kCPU or at::kCUDA.
	 * Returns the result as a Tensor.
	 */
	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& ml, const torch::Tensor* op= nullptr, at::DeviceType dev = at::kCPU, const bool chain_opt = false,  const bool contiguous_dense_to_torch = false);

	torch::Tensor tensor_chain_mul_opt(const std::vector<torch::Tensor>& ml, const torch::Tensor* op, at::DeviceType dev = at::kCPU);

	/**
	 * Computes the matrix chain product of ml and applies it optionally to the matrix op if provided.
	 *
	 * This function converts all the matrices to Tensors before and then computes the tensor product.
	 *
	 *\note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 * \param on_gpu true ot use the GPU backend, false for the CPU backend (false by default).
	 * \param clone true to copy the Faust matrices data to create the tensors, false to use the same data without copying (false by default).
	 * \param contiguous_dense_to_torch if true then torch::chain_matmul is used to computed intermediary product of dense contiguous factors. Note that if chain_opt is true, this option can't be true and is internally set to false silently.
	 * Returns the result as a Faust::MatDense.
	 */
	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op = nullptr, const bool on_gpu = false, const bool clone = false, const bool chain_opt = false,  const bool contiguous_dense_to_torch = false);

	/**
	 * Computes the matrix chain product of tl and applies it optionally to the matrix op if provided.
	 *
	 * op is converted to Tensor before computing the tensor product.
	 *
	 *\note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 * \param on_gpu true ot use the GPU backend, false for the CPU backend (false by default).
	 * \param clone true to copy the Faust matrices data to create the tensors, false to use the same data without copying (false by default).
	 * \param contiguous_dense_to_torch if true then torch::chain_matmul is used to computed intermediary product of dense contiguous factors. Note that if chain_opt is true, this option can't be true and is internally set to false silently.
	 * Returns the result as a Faust::MatDense.
	 */
	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<torch::Tensor>& tl, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch);

}
#include "faust_torch.hpp"
#endif
