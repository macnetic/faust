% ==================================================
%> @brief Functor for the SPSYMM projector.
%>
%> A, the image matrix of M, is such that A is symmetric and \f$ k \le \| A \|_0 \le k + 1,  \| A\|_F = 1 \f$ (if normalized == True), assuming that \f$\| M \|_0 >= k\f$.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b spsymm(@b shape, @b k): returns a SPSYMM projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), k defines the sparsity of the output matrix (k nnz coefficients).<br/>
%> &nbsp;&nbsp;&nbsp; @b spsymm(@b shape,@b k,@b 'normalized', bool,@b 'pos', bool): the optional parameters are set. By default both normalized and pos are false.
%>
%> @param shape: vector of size 2, to define the size of the input matrix.
%> @param k: the sparsity parameter (the number of nonzeros of the projection image. The result might
%> be k+1 nonzeros in case of an odd number of nonzeros on the diagonal.
%> @param 'normalized', true: normalizes the projection image according to its Frobenius norm.
%> @param 'normalized', false: (the default) no normalization.
%> @param 'pos', true: skips the negative values (replaced by zero) of the input matrix.
%> @param 'pos', false: (the default) negative values are not skipped.
%>
%> @retval proj the spsymm projector.
%>
%> <br/>
%>
%> @b Example
%> @code
%> >> import matfaust.proj.spsymm
%> >> rng(42) % for reproducibility
%> >> M = rand(5, 5);
%> >> p = spsymm(size(M), 2);
%> >> p(M)
%>
%> ans =
%>
%>          0         0         0         0         0
%>          0         0    0.7071         0         0
%>          0    0.7071         0         0         0
%>          0         0         0         0         0
%>          0         0         0         0         0
%>
%> >>
%> @endcode
% ==================================================
classdef spsymm < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = spsymm(shape, k, varargin)
			import matfaust.factparams.ConstraintInt
			% default values
			proj.constraint = ConstraintInt('spsymm', shape(1), shape(2), k, varargin{:});
		end
	end
end
