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
		%> see matfaust.factparams.StoppingCriterion
		num_its
		%> see matfaust.factparams.StoppingCriterion
		tol
		%> see matfaust.factparams.StoppingCriterion
		maxiter
		%> see matfaust.factparams.StoppingCriterion
		relerr
	end
	properties(SetAccess = private)
		is_criterion_error
	end
	properties(Constant)
		DEFAULT_MAXITER = 10000
		DEFAULT_NUMITS = 500
		DEFAULT_TOL = .3
	end
	methods
		% =========================================================
		%>	@brief Displays a StoppingCriterion instance.
		%>
		% =========================================================
		function display(self)
			if(self.is_criterion_error)
				fprintf('tol: %f\n', self.tol)
				fprintf('relerr: %d\n', self.relerr)
				fprintf('maxiter: %d\n',self.maxiter)
			else
				fprintf('num_its: %d\n', self.num_its)
				fprintf('maxiter: %d\n', self.maxiter)
			end
		end
		% =========================================================
		%>	@brief Constructor.
		%>
		%>	@param 'numits', int the fixed number of iterations of the algorithm. By default the value is DEFAULT_NUMITS. If arguments num_its and tol are used together only the last one in the argument list is taken into account.
		%>	@param 'tol', real error target according to the algorithm is stopped. If arguments num_its and tol are used together only the last one in the argument list is taken into account.
		%> @param 'maxiter', int the maximum number of iterations to run the algorithm, whatever is the criterion used (tol or num_its).
		%> @param 'relerr', bool false to define a absolute error with tol, true for a relative error
		%> (in this case the 'relmat' matrix will be used to convert internally the given 'tol' to the corresponding absolute error).
		%> @param 'relmat', matrix the matrix against which is defined the relative error. if relerr is True, this argument is mandatory.
		%>
		%> @b Example:
		%> @code
		%>>> import matfaust.factparams.StoppingCriterion
		%>>> s = StoppingCriterion(5)
		%>num_its 5:
		%>maxiter 10000:
		%>>> s = StoppingCriterion(.5)
		%>tol: 0.500000
		%>relerr: 0
		%>maxiter: 10000
		%>>> s = StoppingCriterion('tol', .5)
		%>tol: 0.500000
		%>relerr: 0
		%>maxiter: 10000
		%>>> s = StoppingCriterion('numits', 5)
		%>num_its 500:
		%>maxiter 10000:
		%>>> s = StoppingCriterion('tol', .2, 'relerr', true, 'relmat', rand(10,10))
		%>tol: 1.210149
		%>relerr: 1
		%>maxiter: 10000
		%> @endcode
		% =========================================================
		function stop_crit = StoppingCriterion(varargin)
			%set default values
			is_criterion_error = false;
			tol = matfaust.factparams.StoppingCriterion.DEFAULT_TOL;
			num_its = matfaust.factparams.StoppingCriterion.DEFAULT_NUMITS;
			maxiter = matfaust.factparams.StoppingCriterion.DEFAULT_MAXITER;
			relerr = false;
			matrix = [];
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
									num_its = real(floor((varargin{i+1}))); % real in case of cplx num
									is_criterion_error = false;
								end
							case 'relerr'
								if(nargin == i || ~ islogical(varargin{i+1}))
									error('relerr keyword argument is not followed by a logical')
								else
									relerr = varargin{i+1};
								end
							case 'relmat'
								if(nargin == i || ~ ismatrix(varargin{i+1}))
									error('mat keyword argument is not followed by a matrix')
								else
									matrix = varargin{i+1};
									varargin{i+1} = 0; % erase matrix in args because switch can't evaluate a matrix
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
						tol = tol * norm(matrix, 'fro');
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
