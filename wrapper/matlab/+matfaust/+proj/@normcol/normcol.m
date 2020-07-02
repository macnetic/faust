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
		%> &nbsp;&nbsp;&nbsp; @b normcol(shape, s): returns a NORMCOL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), s defines the sparsity of the output matrix (the 2-norm of each column).<br/>
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param 's', value: (optional) the sparsity parameter, the 2-norm of the column (1 by default).
		%===============================================
		function proj = normcol(shape, varargin)
			import matfaust.factparams.ConstraintReal
			% TODO: factorize this sparsing with normlin
			if(length(varargin) > 0)
				if(strcmp(varargin{1}, 's'))
					disp('s in varargin')

					if(length(varargin) < 2)
						error('s value is missing')
					end
					s = varargin{2};
					varargin(1:2) = []; % deleting the used cell
					if(~ isreal(s) || ~ isscalar(s))
						error('s must be a real scalar')
					end
					if(s < 0)
						error('A norm can''t be negative')
					end
				% else % do nothing (ConstraintReal role)
				end
			end
			proj.constraint = ConstraintReal('normcol', shape(1), shape(2), s)%, varargin{:});
		end
	end
end
