
%> \brief Enumeration class of all matrix chain multiplication methods available to multiply a Faust to a matrix or to compute Faust.full().
%> These methods are used by Faust.optimize_time().
% =========================================================
classdef FaustMulMode
	properties(SetAccess = public, Constant)
        %> \brief The default method, it computes the product from the right to the left.
		DEFAULT_L2R=0
        %% \brief This method follows a greedy principle: it chooses to multiply the less costly product of two matrices at each step of the whole product computation.
        %%
        %% The computational cost depends on the matrix dimensions and the number of nonzeros (when a matrix is in sparse format).
		GREEDY=4
        %> \brief This method implements the classic dynamic programming
        %> solution of the chain matrix problem (see
        %> https://en.wikipedia.org/wiki/Matrix_chain_multiplication#A_dynamic_programming_algorithm).
        DYNPROG=5
		%> \brief This method computes the product of the matrix chain from the left to the right using the Torch C++ library (CPU backend).
		%>
		%> This method is only available for the specific packages pyfaust_torch.
		TORCH_CPU=8
		%> \brief This method computes the product following the minimal cost order using the Torch C++ library (CPU backend).
		%>
        %> The method is basically the same as DYNPROG but it is implemented by the Torch library (chain_matmul function).
        %>
		%> References:
        %> https://pytorch.org/cppdocs/api/function_namespaceat_1aee491a9ff453b6033b4106516bc61a9d.html?highlight=chain_matmul
        %> https://pytorch.org/docs/stable/generated/torch.chain_matmul.html?highlight=chain_matmul#torch.chain_matmul
        %>
		%> This method is only available for the specific packages pyfaust_torch.
		TORCH_CPU_BEST_ORDER=9
		%> \brief The same as TORCH_CPU except that torch::chain_matmul is used to
		%> compute in one call the intermediary product of dense contiguous
		%> factors, then the result is multiplied by sparse factors if any remains.
		%>
		%> This method is only available for the specific packages pyfaust_torch.
		TORCH_CPU_DENSE_ROW_TORCH=10
	end
end
