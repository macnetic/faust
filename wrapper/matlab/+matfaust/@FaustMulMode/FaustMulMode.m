
%> \brief Enumeration class of all matrix chain multiplication methods available to multiply a Faust to a matrix or to compute Faust.full().
%> These methods are used by Faust.optimize_time().
% =========================================================
classdef FaustMulMode
	properties(SetAccess = public, Constant)

		%> \brief The default method. Multiplying from the right to the left.
		DEFAULT_L2R=0
		%> \brief This method computes the product by its ends.
		%>
		%> For each iteration/multiplication it chooses to multiply the most right
		%> or the most left pair of matrices (in order to decrease the computation cost).
		%> The computational cost depends on the matrix dimensions.
		GREEDY_ALL_ENDS=1
		%> \brief This method computes the product starting by the pair of matrices whose the computation cost is the smallest.
		%>
		%> After this first multiplication the rest of the factors are multiplied
		%> on the left side (from the right to the left) and then the right factors are
		%> multiplied (from the left to the right).
		%> The computational cost depends on the matrix dimensions.
		GREEDY_1ST_BEST=2
		%> \brief This method computes a chain of matrices by ordering the product according to the minimal computation cost order.
		%>
		%> The computational cost depends on the matrix dimensions.
		GREEDY_ALL_BEST_CONVDENSE=3
		%> \brief This method follows the same principle as GREEDY_ALL_BEST_CONVDENSE method but is capable to multiply dense matrices as well as sparse matrices.
		%>
		%> The computational cost depends on the matrix dimensions and the number
		%> of nonzeros (when a matrix is in sparse format).
		GREEDY_ALL_BEST_GENMAT=4
        %> \brief This method implements the classic dynamic programming
        %> solution of the chain matrix problem (see
        %> https://en.wikipedia.org/wiki/Matrix_chain_multiplication#A_dynamic_programming_algorithm).
        DYNPROG=5
		%> \brief This method computes the product performing a parallel reduction of the product.
		%>
		%> It uses as many threads as C++ STL advises (std::thread::hardware_concurrency() -- https://en.cppreference.com/w/cpp/thread/thread/hardware_concurrency).
		%>
		%> Reference: https://en.wikipedia.org/wiki/Reduce_%29parallel_pattern%29
		CPP_PROD_PAR_REDUC=6
		%> \brief This method is equivalent to CPP_PROD_PAR_REDUC but is implemented using OpenMP.
		OMP_PROD_PAR_REDUC=7
		%> \brief This method computes the product of the matrix chain from the left to the right using Torch C++ library (CPU backend).
		%>
		%> This method is only available for the specific packages pyfaust_torch.
		TORCH_CPU=8
		%> \brief This method computes the product following the minimal cost order using Torch C++ library (CPU backend).
		%>
		%> The method is basically the same as DYNPROG but it is implemented with Torch library.
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
