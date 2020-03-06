% =========================================================
%> @brief This class defines a StoppingCriterion for the FAuST's algorithms.
%>
%> A stopping criterion can be of two kinds:
%>
%>     - number of iterations,
%>     - error treshold for the approximation of the matrix.
%>
% =========================================================
classdef StoppingCriterion
	properties (SetAccess = public)
		is_criterion_error
		num_its
		tol
		maxiter
		relerr
	end
	methods
		function stop_crit = StoppingCriterion(varargin)
			%set default values
			is_criterion_error = false;
			tol = 0.3;
			num_its = 500;
			maxiter = 10000;
			relerr = false;
			matrix = []
			if(nargin < 1)
				error('matfaust.factparams.StoppingCriterion needs at least one argument.')
			else
				if(nargin == 1)
					% only one arg
					if(isscalar(varargin{1}))
						if(floor(varargin{1}) == varargin{1})
							% arg is int: num_its
							num_its = real(varargin{1});
						else
							% real
							tol = real(varargin{1});
							is_criterion_error = true;
						end
					else
						error('when only argument is passed to StoppingCriterion it must be a scalar')
					end
				else
					for i=1:nargin
						switch(varargin{i})
							case 'maxiter'
								if(nargin == i || ~ isscalar(varargin{i+1}))
									error('maxiter keyword arg. is not followed by a number')
								else
									maxiter = real(floor((varargin{i+1}))); % real in case of cplx num
								end
							case 'tol'
								if(nargin == i || ~ isscalar(varargin{i+1}))
									error('tol keyword arg. is not followed by a number')
								else
									tol = real(varargin{i+1}); % real in case of cplx num
									is_criterion_error = true;
								end
							case 'num_its'
								if(nargin == i || ~ isscalar(varargin{i+1}))
									error('num_its keyword arg. is not followed by a number')
								else
									num_its = real(floor((varargin{i+1}))) % real in case of cplx num
									is_criterion_error = false;
								end
							case 'relerr'
								if(nargin == i || ~ islogical(varargin{i+1}))
									error('relerr keyword argument is not followed by a logical')
								else
									relerr = varargin{i+1}
								end
							case 'relmat'
								if(nargin == i || ~ ismatrix(varargin{i+1}))
									error('mat keyword argument is not followed by a matrix')
								else
									matrix = varargin{i+1}
									varargin{i+1} = 0 % erase matrix in args because switch can't evaluate a matrix
									% so otherwise it would fail the next ite
								end
							end
						end
					end
					if(relerr)
						if(size(matrix, 1) == 0 && size(matrix, 2) == 0)
							error('when error is relative (relerr == true) the reference matrix ''relmat'' must be specified')
						end
						% adjust tol to relerr
						tol = tol * norm(matrix, 'fro')
					end
					stop_crit.is_criterion_error = is_criterion_error;
					stop_crit.num_its = num_its;
					stop_crit.maxiter = maxiter;
					stop_crit.tol = tol;
					stop_crit.relerr = relerr;
				end
			end
		end
	end
