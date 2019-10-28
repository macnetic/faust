%==========================================================================================
%> @brief Computes the FGFT of the Laplacian matrix Lap (using fact.eigtj).
%>
%>
%> @b Usage
%>  &nbsp;&nbsp;&nbsp; @b See fact.eigtj
%>
%> @param Lap the Laplacian matrix (which is symmetric or hermitian).
%> @param maxiter see fact.eigtj
%> @param 'nGivens_per_fac', integer see fact.eigtj
%> @param 'tol', number see fact.eigtj
%> @param 'relerr', bool see fact.eigtj
%> @param 'verbosity', integer see fact.eigtj
%> @param 'order', integer see fact.eigtj
%>
%> @retval [FGFT,D]:
%> - with FGFT being the Faust object representing the Fourier transform and,
%> -  D as a sparse diagonal matrix of the eigenvalues by default in ascendant order along the rows/columns.
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420 <https://hal.inria.fr/hal-01416110>
%>
%>
%> <p> @b See @b also fact.eigtj, fact.fgft_palm
%>
%==========================================================================================
function [FGFT,D] = fgft_givens(Lap, maxiter, varargin)
	import matfaust.Faust
	nGivens_per_fac = 1; % default value
	verbosity = 0; % default value
%	if(~ ismatrix(Lap) || ~ isreal(Lap))
%		error('Lap must be a real matrix.')
%	end
	if(size(Lap,1) ~= size(Lap,2))
		error('Lap must be square')
	end
	if(~ isnumeric(maxiter) || maxiter-floor(maxiter) > 0 || maxiter <= 0)
		error('maxiter must be a positive integer.')
	end
	bad_arg_err = 'bad number of arguments.';
	tol = 0;
	relerr = true;
	verbosity = 0;
	argc = length(varargin);
	order = 1; % ascending order
	if(argc > 0)
		for i=1:argc
			switch(varargin{i})
				case 'tol'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('tol keyword arg. is not followed by a number')
					else
						tol = real(varargin{i+1}); % real in case of cplx num
					end
				case 'relerr'
					if(argc == i || ~ islogical(varargin{i+1}))
						error('relerr keyword argument is not followed by a logical')
					else
						relerr = varargin{i+1};
					end
				case 'verbosity'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('verbose keyword argument is not followed by a number')
					else
						verbosity = floor(real(varargin{i+1}));
					end
				case 'nGivens_per_fac'
					if(argc == i || ~ isscalar(varargin{i+1}) || ~ isnumeric(varargin{i+1}))
						error('nGivens_per_fac must be followed by a positive integer.')
					else
						nGivens_per_fac = floor(abs(real(varargin{i+1})));
						nGivens_per_fac = min(nGivens_per_fac, maxiter);
						nGivens_per_fac = max(1, nGivens_per_fac);
					end
				case 'order'
					if(argc == i || (~ strcmp(varargin{i+1}, 'ascend') && ~ strcmp(varargin{i+1}, 'descend') && ~ strcmp(varargin{i+1}, 'undef')))
						error('order must be followed by a char array among ''ascend'', ''descend'' or ''undef''.')
					else
						order = varargin{i+1};
						if(order(1) == 'a')
							order = 1
						elseif(order(1) == 'd')
							order = -1
						else
							order = 0
						end
					end
				otherwise
					if(isstr(varargin{i}) && (~ strcmp(varargin{i}, 'ascend') && ~ strcmp(varargin{i}, 'descend') && ~ strcmp(varargin{i}, 'undef')) )
						error([ varargin{i} ' unrecognized argument'])
					end
			end
		end
	end
	[core_obj, D] = mexfgftgivensReal(Lap, maxiter, nGivens_per_fac, verbosity, tol, relerr, order);
	D = sparse(diag(real(D)));
	FGFT = Faust(core_obj, isreal(Lap));
end
