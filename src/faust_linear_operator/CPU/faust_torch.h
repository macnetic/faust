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
	 * \param transpose is true by default because it's more efficient to handle the difference of data storage between torch and faust ; row-major order for torch, column-major order for faust.
	 */
	template<typename FPP, FDevice D>
		void convMatDenseToTensor(const Faust::MatDense<FPP,D> & dm, torch::Tensor & t, at::DeviceType dev = at::kCPU, const bool clone = false, const bool transpose = true);

	/**
	 * Converts a torch::Tensor (t) to Faust::MatDense (dm).
	 *
	 * \param transpose see convMatDenseToTensor.
	 *
	 */
	template<typename FPP, FDevice D>
		void convTensorToMatDense(const torch::Tensor & t, Faust::MatDense<FPP,D> & dm, const bool transpose = true);


	/**
	 * Converts a Faust::MatSparse to a torch::Tensor.
	 *
	 * \note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 *\param dev the device at::kCPU or at::kCUDA.
	 *\param clone true to copy the Faust matrix data to create the tensor, false to use the same data without copying (false by default).
	 *\param transpose is true by default because it's more efficient to handle the difference of data storage between torch and faust ; row-major order for torch, column-major order for faust. For MatSparse it does not really matter (because data copying is still necessary) but it's preferable to be consistent with what is done for MatDense.
	 */
	template<typename FPP, FDevice D>
		void convMatSparseToTensor(const Faust::MatSparse<FPP,D> & spm, torch::Tensor & t, at::DeviceType dev, const bool clone, const bool transpose = true);
	/**
	 * Converts a Faust::MatGeneric vector to a torch::TensorList (vector alias).
	 *
	 *\param dev the device at::kCPU or at::kCUDA.
	 * \param clone true to copy the Faust matrices data to create the tensors, false to use the same data without copying (false by default).
	 * \param transpose to true implies that the ml factors will be converted and stored into tl in reverse order (the goal is to compute efficiently the ml product -- cf. tensor_chain_mul, through tl defined as the transpose product. It's more efficient).
	 */
	template<typename FPP, FDevice D>
		void convMatGenListToTensorList(const std::vector<Faust::MatGeneric<FPP,D>*> & ml, std::vector<torch::Tensor> &tl, at::DeviceType dev = at::kCPU, const bool clone = false, const bool transpose = true);

	/**
	 * Computes the tensor chain product of tl and applies it optionally to the tensor op.
	 *
	 * \param tl the sequence of tensors to compute the product.
	 * \param if op is not nullptr the functions returns the product of tl and op.
	 *\param dev the device at::kCPU or at::kCUDA.
	 * \param chain_opt if true then the function pass the hand to tensor_chain_mul_opt.
	 * \param contiguous_dense_to_torch If true consecutive/contiguous dense factors will be computed (as intermediary product) through torch::chain_matmul() as it is always done when tl is full of only dense factors (e.g.: if tl = {S1, S2, D3, D4, D5, S6} and the letter D represents a dense factor, S a sparse factor, D3*D4*D4 will be calculated in one call of torch::chain_matmul while the remaining products will be computed one by one. This option and chain_opt are exclusive. If chain_opt is true, this boolean is forced to false. Nota: as far as I've tested torch::chain_matmul can't work with sparse Tensor-s (hence these option and function).
	 * \param op_on_left: if true op*tl is computed, otherwise tl*op is computed. It is set to true by default because this is the optimal scenario (the transpose product scenario) regarding the column-major order of Faust::Matdense and the row-major order of torch. This default value is consistent with the default value of transpose in other functions.
	 *
	 * Returns the result as a Tensor.
	 */
	torch::Tensor tensor_chain_mul(const std::vector<torch::Tensor>& tl, const torch::Tensor* op= nullptr, at::DeviceType dev = at::kCPU, const bool chain_opt = false,  const bool contiguous_dense_to_torch = false, const bool op_on_left = true);

	/**
	 * Computes the tensor chain product of tl and applies it optionally to the tensor op.
	 *
	 * This function does the same as tensor_chain_mul except that it optimizes the matrix chain product choosing an order of computation that minimizes the cost.
	 *
	 * \param tl the sequence of tensors to compute the product.
	 * \param if op is not nullptr the functions returns the product of tl and op.
	 *\param dev the device at::kCPU or at::kCUDA.
	 * \param op_on_left: if true op*tl is computed, otherwise tl*op is computed. It is set to true by default because this is the optimal scenario (the transpose product scenario) regarding the column-major order of Faust::Matdense and the row-major order of torch. This default value is consistent with the default value of transpose in other functions.
	 *
	 * Returns the result as a Tensor.
	 */
	torch::Tensor tensor_chain_mul_opt(const std::vector<torch::Tensor>& tl, const torch::Tensor* op, at::DeviceType dev = at::kCPU, const bool op_on_left = true);

	/**
	 * Computes the matrix chain product of ml and applies it optionally to the matrix op if provided.
	 *
	 * This function converts all the matrices to Tensors before and then computes the tensor product (using tensor_chain_mul(TensorList, Tensor) just above).
	 *
	 * \note Complex tensors are not available in libtorch, an exception is thrown when FPP is complex.
	 *
	 * \param on_gpu true ot use the GPU backend, false for the CPU backend (false by default).
	 * \param clone true to copy the Faust matrices data to create the tensors, false to use the same data without copying (false by default).
	 * \param contiguous_dense_to_torch if true then torch::chain_matmul is used to computed intermediary product of dense contiguous factors. Note that if chain_opt is true, this option can't be true and is internally set to false silently.
	 * Returns the result as a Faust::MatDense.
	 */
	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<Faust::MatGeneric<FPP,D>*>& ml, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op = nullptr, const bool on_gpu = false, const bool clone = false, const bool chain_opt = false,  const bool contiguous_dense_to_torch = false, const bool transpose = true);

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
	 *
	 * Returns the result as a Faust::MatDense.
	 */
	template<typename FPP, FDevice D>
		void tensor_chain_mul(const std::vector<torch::Tensor>& tl, Faust::MatDense<FPP,Cpu> & out, const Faust::MatGeneric<FPP,D>* op, const bool on_gpu, const bool clone, const bool chain_opt, const bool contiguous_dense_to_torch, const bool transpose = true);

	/**
	 * This function display a Tensor list (size and storage format of Tensor-s, DENSE or SPARSE).
	 *
	 * \param transpose if true the transpose TensorList of tl is displayed. It's true by default to be consistent with other functions.
	 */
	void display_TensorList(std::vector<torch::Tensor>& tl, const bool transpose = true);
}
#include "faust_torch.hpp"
#endif
