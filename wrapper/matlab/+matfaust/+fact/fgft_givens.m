%==========================================================================================
%> @brief Computes the FGFT for the Laplacian matrix Lap.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b fgft_givens(Lap, J) calls the non-parallel Givens algorithm.<br/>
%> &nbsp;&nbsp;&nbsp; @b fgft_givens(Lap, J, 0) or fgft_givens(Lap, J, 1) do the same as in previous line.<br/>
%> &nbsp;&nbsp;&nbsp; @b fgft_givens(Lap, J, t) calls the parallel Givens algorithm (if t > 1, otherwise it calls basic Givens algorithm), see eigtj. <br/>
%> &nbsp;&nbsp;&nbsp; @b fgft_givens(Lap, J, t, 'verbosity', 2) same as above with a level of verbosity of 2 in output. <br/>
%>
%> @param Lap the Laplacian matrix as a numpy array. Must be real and symmetric.
%> @param J see eigtj
%> @param t see eigtj
%> @param verbosity see eigtj
%>
%> @retval [FGFT,D]:
%> - with FGFT being the Faust object representing the Fourier transform and,
%> -  D as a sparse diagonal matrix of the eigenvalues in ascendant order along the rows/columns.
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420
%>
%>
%> <p> @b See @b also fact.eigtj, fact.fgft_palm
%>
%==========================================================================================
function [FGFT,D] = fgft_givens(Lap, J, varargin)
	import matfaust.Faust
	t = 1; % default value
	verbosity = 0; % default value
	if(~ ismatrix(Lap) || ~ isreal(Lap))
		error('Lap must be a real matrix.')
	end
	if(size(Lap,1) ~= size(Lap,2))
		error('Lap must be square')
	end
	if(~ isnumeric(J) || J-floor(J) > 0 || J <= 0)
		error('J must be a positive integer.')
	end
	bad_arg_err = 'bad number of arguments.';
	if(length(varargin) >= 1)
		t = varargin{1};
		if(~ isnumeric(t))
			error('t must be a positive or nul integer.')
		end
		t = floor(abs(t));
		t = min(t, J);
		if(length(varargin) >= 2)
			if(~ strcmp(varargin{2}, 'verbosity'))
				error('arg. 4, if used, must be the str `verbosity''.')
			end
			if(length(varargin) == 3)
				if(isnumeric(varargin{3}))
					verbosity = floor(real(varargin{3}));
				else
					error('verbosity must be numeric')
				end
			else
				error(bad_arg_err)
			end
		end
	end
	[core_obj, D] = mexfgftgivensReal(Lap, J, t, verbosity);
	D = sparse(diag(D));
	FGFT = Faust(core_obj, true);
end
