%==================================================
%> @brief Functor for the NORMCOL projector.
%>
%> A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \f$ the j-th column \f$  A_{*,j} \f$ is such that \f$ \| A_{*,j} \|_2 = s  \f$.
%==================================================
classdef normcol < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b normcol(shape, normval): returns a NORMCOL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), normval defines the sparsity of the output matrix (the 2-norm of each column).<br/>
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param normval: the sparsity parameter.
		%===============================================
		function proj = normcol(shape, normval)%, varargin)
			import matfaust.factparams.ConstraintReal
			proj.constraint = ConstraintReal('normcol', shape(1), shape(2), normval)%, varargin{:});
		end
	end
end
