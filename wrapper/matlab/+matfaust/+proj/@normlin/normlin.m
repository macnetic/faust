%==================================================
%> @brief Functor for the NORMLIN projector.
%>
%> A, the image matrix, is defined by \f$ \forall i \in \{1,\ldots,shape(1)\} \f$ the i-th row \f$  A_{i,*} \f$ is such that \f$ \| A_{i, *} \|_2 = s \f$.
%==================================================
classdef normlin < matfaust.proj.proj_gen
	properties
	end
	methods
		%===============================================
		%> @b Usage
		%>
		%> &nbsp;&nbsp;&nbsp; @b normlin(shape, normval): returns a NORMLIN projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), normval defines the sparsity of the output matrix (the 2-norm of each row).<br/>
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param normval: the sparsity parameter.
		%===============================================
		function proj = normlin(shape, normval, varargin)
			import matfaust.factparams.ConstraintReal
			if(normval < 0)
				error('A norm can''t be negative')
			end
			proj.constraint = ConstraintReal('normlin', shape(1), shape(2), normval, varargin{:});
		end
	end
end
