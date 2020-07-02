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
		%> &nbsp;&nbsp;&nbsp; @b normlin(shape, s): returns a NORMLIN projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), s defines the sparsity of the output matrix (the 2-norm of each row).<br/>
		%>
		%> @param shape: vector of size 2, to define the size of the input matrix.
		%> @param 's', value: (optional) the sparsity parameter, the 2-norm of the row (1 by default).
		%===============================================
		function proj = normlin(shape, varargin)
			import matfaust.factparams.ConstraintReal
			s = 1;
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

			proj.constraint = ConstraintReal('normlin', shape(1), shape(2), s, varargin{:});
		end
	end
end
