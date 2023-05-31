% ==================================================
%> @brief Functor for the NORMCOL projector.
%>
%> A, the image matrix, is defined by \f$ \forall j \in \{1,\ldots,shape(2)\} \f$ the j-th column \f$  A_{*,j} \f$ is such that \f$ \| A_{*,j} \|_2 = s  \f$.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b normcol(shape, s): returns a NORMCOL projector (functor), shape defines the size of the input matrix (e.g. [1, 10]), s defines the sparsity of the output matrix (the 2-norm of each column).<br/>
%>
%> @param shape: vector of size 2, to define the size of the input matrix.
%> @param 's', value: (optional) the sparsity parameter, the 2-norm of the column (1 by default).
%> @retval normcol projector
%>
%> @b Example
%> @code
%> >> import matfaust.proj.normcol
%> >> rng(42)
%> >> M = rand(5, 5)
%>
%> M =
%>
%>     0.3745    0.1560    0.0206    0.1834    0.6119
%>     0.9507    0.0581    0.9699    0.3042    0.1395
%>     0.7320    0.8662    0.8324    0.5248    0.2921
%>     0.5987    0.6011    0.2123    0.4319    0.3664
%>     0.1560    0.7081    0.1818    0.2912    0.4561
%>
%> >> p = normcol(size(M), .01);
%> >> pM = p(M);
%> >> norm(pM(:,1), 2) % doctest: +ELLIPSIS
%>
%> ans =
%>
%>     1...
%>
%> >>
%> @endcode
% ==================================================
classdef normcol < matfaust.proj.proj_gen
	properties
	end
	methods

		function proj = normcol(shape, varargin)
			import matfaust.factparams.ConstraintReal
			s = 1;
			% TODO: factorize this parsing with normlin
			if(length(varargin) > 0)
				if(strcmp(varargin{1}, 's'))
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
			proj.constraint = ConstraintReal('normcol', shape(1), shape(2), s); %, varargin{:});
		end
	end
end
