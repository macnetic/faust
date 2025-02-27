%===============================================================================
%>  Runs the greedy OMP algorithm optimized by Cholesky decomposition.
%===
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b x=omp(y, D) runs the algorithm with the default configuration which is equivalent to <b> x=omp(y, D, 'maxiter', length(y), 'verbose', false)</b>.<br/>
%> &nbsp;&nbsp;&nbsp; <b> x=omp(y, D, 'maxiter', N)</b> stops the algorithm after N iterations.<br/>
%> &nbsp;&nbsp;&nbsp; <b> x=omp(y, D, 'tol', 10^-16)</b> runs the algoritm until the relative error is lower or equal to 10**-16. This is equivalent to <b>x=omp(y, D, 'tol', 10^-16, 'relerr', true)</b><br/>
%> &nbsp;&nbsp;&nbsp; <b> x=omp(y, D, 'maxiter', N, 'tol', 10^-16)</b> runs at most N iterations until the <code>tol</code> precision is reached.<br/>
%> &nbsp;&nbsp;&nbsp; <b>x=omp(y, D, 'tol', 10^-16, 'relerr', false)</b> runs the algorithm until the absolute error is lower or equal to 10**-16.<br/><br/>
%>
%>
%> @param y The vector to approximate by D*x.
%> @param D The dictionary as a matrix or a Faust.
%> @param 'maxiter', number (optional) To define a maximum number of iterations (by default this length(y)).
%> @param 'tol', number (optional) To define the error value used for the relative or absolute error (by default this is 0, for not stopping on any error criterion).
%> @param 'relerr', true (optional) To define a stopping criterion based on the relative error (this is the default error).
%> @param 'relerr', false (optional) To define a stopping criterion based on the absolute error.
%> @param 'verbose', true (optional) To enable the verbosity.
%> @param 'verbose', false (optional) To disable the verbosity (this is the default option).
%>
%>
%> @return x the solution of y = D*x (according to the error).
%>
%> @b Example:
%> @code
%> >> rng(42)
%> >> F = matfaust.rand(15, 10);
%> >> x = rand(10, 1)
%>
%> x =
%>
%>     0.3745
%>     0.9507
%>     0.7320
%>     0.5987
%>     0.1560
%>     0.1560
%>     0.0581
%>     0.8662
%>     0.6011
%>     0.7081
%>
%> >> y = F * x;
%> >> x_ = matfaust.tools.omp(y, F)
%> Stopping. Exact signal representation found!
%>
%> x_ =
%>
%>     0.3745
%>     0.9507
%>     0.7320
%>     0.5987
%>     0.1560
%>     0.1560
%>     0.0581
%>     0.8662
%>     0.6011
%>     0.7081
%>
%> @endcode
%>
%===============================================================================
function x = omp(y, D, varargin)
	argc = length(varargin);
	% set parameter default values
	maxiter = length(y);
	tol = 0;
	relerr = true;
	verbose = false;
	if(argc > 0)
		for i=1:argc
			switch(varargin{i})
				case 'maxiter'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('maxiter keyword arg. is not followed by a number');
					else
						maxiter = real(floor((varargin{i+1}))); % real in case of cplx num
					end
				case 'tol'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('tol keyword arg. is not followed by a number');
					else
						tol = real(varargin{i+1}); % real in case of cplx num
					end
				case 'relerr'
					if(argc == i || ~ islogical(varargin{i+1}))
						error('relerr keyword argument is not followed by a logical');
					else
						relerr = varargin{i+1};
					end
				case 'verbose'
					if(argc == i || ~ islogical(varargin{i+1}))
						error('verbose keyword argument is not followed by a logical');
					else
						verbose = varargin{i+1};
					end
				otherwise
					if(isstr(varargin{i}))
						error([ varargin{i} ' unrecognized argument']);
					end
			end
		end
	end
	if(tol > 0)
		tol
		if(relerr)
			if(size(y,1) == 1)
				y_sqr_norm = y*y';
			else
				y_sqr_norm = y'*y;
			end
			x = greed_omp_chol(y, D, size(D,2), 'stopCrit', 'mse',  'stopTol', tol*y_sqr_norm/length(y), 'verbose', verbose);
		else % absolute error
			x = greed_omp_chol(y, D, size(D,2), 'stopCrit', 'mse',  'stopTol', tol^2/length(y), 'verbose', verbose);
		end
	else
		% maxiter
		x = greed_omp_chol(y, D, size(D,2), 'stopCrit', 'M', 'stopTol', maxiter, 'verbose', verbose);
	end
end
%> @package matfaust.tools @brief The matfaust tools namespace


