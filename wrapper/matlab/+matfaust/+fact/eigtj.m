%==========================================================================================
%> @brief Performs an approximate eigendecomposition of M and returns the eigenvalues in D along with the corresponding left eigenvectors (as the columns of the Faust object V).
%>
%> The output is such that <code>V*D*V'</code> approximates M. V is a product of Givens rotations obtained by truncating the Jacobi algorithm.
%>
%> The trade-off between accuracy and sparsity of V can be set through the parameters nGivens and nGivens_per_fac or concurrently with the arguments tol and relerr that define the targeted error.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; <b>[V,D] = eigtj(M, 'nGivens', n)</b> computes V with n Givens rotations grouped into factors containing ideally <code>m = floor(size(M,1)/2)</code> Givens rotation each (if m doesn't divide n, the first factor at least contains less than m rotations).<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, 'nGivens', n, 'nGivens_per_fac', 1) does the same as in previous line but limiting to 1 the number of Givens rotation per factor. Setting nGivens_per_fac to 0 produces the same results. <br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, 'nGivens', n, 'nGivens_per_fac', t) as above but with t Givens rotatons per factor.<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, 'nGivens', n, 'nGivens_per_fac', t, 'verbosity', 2) same as above with a level of output verbosity of 2. <br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, 'nGivens', n, 'nGivens_per_fac', t, 'tol', 0.01, 'relerr', true) Uses a stopping criterion based on relative error <code>norm(V*D*V'-M, 'fro')/norm(M, 'fro')</code>. This criterion is concurrent to nGivens (here n).<br/>
%> &nbsp;&nbsp;&nbsp; @b eigtj(M, 'nGivens', n, 'nGivens_per_fac', t, 'tol', 0.01, 'relerr', true) Uses a stopping criterion based on absolute error <code>norm(V*D*V'-M, 'fro')</code>. This criterion is concurrent to nGivens (here n).<br/>
%>
%> @param M the matrix to diagonalize. Must be real and symmetric, or complex hermitian. Can be in dense or sparse format.
%> @param nGivens, integer (optional if tol is set) the maximum number of iterations which is defined by the
%> number of Givens rotations that can be computed in eigenvector transform V.
%> The number of rotations per factor of V is defined by nGivens_per_fac.
%> @param 'tol', number (optional if nGivens is set) the tolerance error at which the algorithm stops. By default, it's zero for not stopping on an error criterion. Note that the error reaching is not guaranteed (in particular, if the error starts to increase from one iteration to another then the algorithm is stopped).
%> @param 'order', char (optional) 'descend' for a descending order of eigenvalues, 'ascend' for an ascending order (default value) or 'undef' for no sort.
%> @param 'nGivens_per_fac', integer (optional) the number of Givens rotations per factor of V, must be an integer between 1 to <code>floor(size(M, 1)/2)</code> which is the default value.
%> @param 'relerr', true (optional) For a stopping criterion based on the relative error (this is the default error).
%> @param 'relerr', false (optional) For a stopping criterion based on the absolute error.
%> @param 'verbosity', integer (optional) the level of verbosity, the greater the value the more info is displayed. It can be helpful to understand for example why the  algorithm stopped before reaching the tol error or the number of Givens (nGivens).
%>
%> @retval [V,D]
%> - V the Faust object representing the approximate eigenvector transform. The column <code>V(:, i)</code> is the eigenvector corresponding to the eigenvalue <code>D(i,i)</code>.
%> - D the sparse diagonal matrix of the approximate eigenvalues (by default in ascendant order along the rows/columns).
%>
%> @b Example
%> @code
%> import matfaust.fact.eigtj
%>
%> % get a Laplacian to diagonalize
%> load('Laplacian_256_community.mat')
%> % do it
%> [Uhat, Dhat] = eigtj(Lap, 'nGivens', size(Lap,1)*100, 'nGivens_per_fac', size(Lap, 1)/2, 'verbosity', 2)
%> % Uhat is the Fourier matrix/eigenvectors approximation as a Faust (200 factors)
%> % Dhat the eigenvalues diagonal matrix approx.
%> % Computing the decomposition of the same matrix but targeting a precise relative error
%> [Uhat2, Dhat2] = eigtj(Lap, 'tol', 0.01)
%> assert(norm(Lap-Uhat2*Dhat2*Uhat2', 'fro')/norm(Lap, 'fro') < .011)
%> % and then asking for an absolute error
%> [Uhat3, Dhat3] = eigtj(Lap, 'tol', 0.1, 'relerr', false)
%> assert(norm(Lap-Uhat3*Dhat3*Uhat3', 'fro') < .11)
%> % now recompute Uhat2, Dhat2 but asking a descending order of eigenvalues
%> [Uhat4, Dhat4] = eigtj(Lap, 'tol', 0.01, 'order', 'descend')
%> assert(all(all(Dhat4(end:-1:1,end:-1:1) == Dhat2(1:end,1:end))))
%> % and now with no sort
%> [Uhat5, Dhat5] = eigtj(Lap, 'tol', 0.01, 'order', 'undef');
%> assert(all(sort(diag(Dhat5)) == diag(Dhat2)))
%> @endcode
%>
%>
%>
%> @b References:
%> - [1]   Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast
%> graph Fourier transforms via multi-layer sparse approximations",
%> IEEE Transactions on Signal and Information Processing
%> over Networks 2018, 4(2), pp 407-420 <https://hal.inria.fr/hal-01416110>
%>
%> <p> @b See @b also fact.svdtj
%>
%==========================================================================================
function [V,D] = eigtj(M, varargin)
	import matfaust.Faust
	nGivens_per_fac = 1; % default value
	verbosity = 0; % default value
%	if(~ ismatrix(M) || ~ isreal(M))
%		error('M must be a real matrix.')
%	end
	if(size(M,1) ~= size(M,2))
		error('M must be square')
	end
	bad_arg_err = 'bad number of arguments.';
	nGivens = 0;
	tol = 0;
	relerr = true;
	verbosity = 0;
	argc = length(varargin);
	order = 1; % ascending order
	if(argc > 0)
		for i=1:argc
			switch(varargin{i})
				case 'nGivens'
					if(argc == i || ~ isnumeric(varargin{i+1}) || varargin{i+1}-floor(varargin{i+1}) > 0 || varargin{i+1} <= 0 )
						error('nGivens keyword arg. is not followed by an integer')
					end
					nGivens = varargin{i+1};
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
						nGivens_per_fac = max(1, nGivens_per_fac);
					end
				case 'order'
					if(argc == i || (~ strcmp(varargin{i+1}, 'ascend') && ~ strcmp(varargin{i+1}, 'descend') && ~ strcmp(varargin{i+1}, 'undef')))
						error('order must be followed by a char array among ''ascend'', ''descend'' or ''undef''.')
					else
						order = varargin{i+1};
						if(order(1) == 'a')
							order = 1;
						elseif(order(1) == 'd')
							order = -1;
						else
							order = 0;
						end
					end
				otherwise
					if(isstr(varargin{i}) && (~ strcmp(varargin{i}, 'ascend') && ~ strcmp(varargin{i}, 'descend') && ~ strcmp(varargin{i}, 'undef')) )
						error([ varargin{i} ' unrecognized argument'])
					end
			end
		end
	end
	if(nGivens == 0 && tol == 0)
		error('Either nGivens or tol must be greater than zero.')
	end
	if(nGivens > 0)
		nGivens_per_fac = min(nGivens_per_fac, nGivens);
	end
	[core_obj, D] = mexfgftgivensReal(M, nGivens, nGivens_per_fac, verbosity, tol, relerr, order);
	D = sparse(diag(real(D)));
	V = Faust(core_obj, isreal(M));
end
