%==========================================================================================
%> @brief Performs an approximate eigendecomposition of M and returns the eigenvalues in D along with the corresponding normalized right eigenvectors (as the columns of the Faust object V).
%>
%> The output is such that <code>V*D*V'</code> approximates M. V is a product of Givens rotations obtained by truncating the Jacobi algorithm.
%>
%> The trade-off between accuracy and complexity of V can be set through the parameters nGivens and tol that define the targeted number of Givens rotations and targeted error.
%>
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; Primary examples of calls include:<br/>
%> &nbsp;&nbsp;&nbsp; <b>[V,D] = eigtj(M, 'nGivens', n)</b> outputs V as a Faust object made of n elementary Givens rotations grouped into factors containing ideally <code>m = floor(size(M,1)/2)</code> such Givens rotations, and the diagonal entries of D are the approximate eigenvalues in ascending order.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D]  =  eigtj(M,’tol’,0.01)</b> same as above with n determined adaptively by a relative approximation error 0.01.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D] = eigtj(M,’tol’,0.01)</b> same as above with n determined adaptively by a relative approximation error 0.01.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D] = eigtj(M,’tol’,0.01,’relerr’,False)</b> same as above with an absolute approximation error.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D] = eigtj(M,’nGivens’,n,’tol’,0.01)</b> same as above with a number of elementary Givens bounded by n even if the targeted approximation error is not achieved.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D]= eigtj(M,’nGivens’,n,’tol’,0.01,’nGivens_per_fac’,t)</b> same as above with (up to) t Givens rotations per factor.<br/>
%>&nbsp;&nbsp;&nbsp; <b>[V,D]= eigtj(M,’nGivens’,n,’order’,’descend’) </b> same as above where the diagonal entries of D are the approximate eigenvalues in descending order (and with columns of V permuted accordingly).<br/>
%>&nbsp;&nbsp;&nbsp; <b>eigtj(M, 'nGivens', n, 'nGivens_per_fac', t, 'tol', 0.01, 'relerr', true)</b> uses a stopping criterion based on absolute error norm(V*D*V'-M, 'fro'). This criterion is concurrent to nGivens (here n).<br/>
%>
%> @param M the matrix to diagonalize. Must be real and symmetric, or complex hermitian. Can be in dense or sparse format. The class(M) value can be double or single (only if isreal(M) is true).
%> @param nGivens, integer [optional if tol is set] targeted number of Givens rotations.
%> The number of rotations per factor of V is defined by nGivens_per_fac.
%> @param 'tol', number [optional if nGivens is set] the tolerance error at which the algorithm stops. The default value is zero so that stopping is based on reaching the targeted nGivens.
%> @param 'err_period', int:  it defines the period, in number of factors of V,
%> the error is compared to tol (reducing the period spares some factors but increases slightly the computational cost because the error
%> is computed more often).
%> @param 'order', char [optional, default is ‘ascend’] order of eigenvalues, possible choices are ‘ascend, 'descend' or 'undef' (to avoid a sorting operation and save some time).
%> @param 'nGivens_per_fac', integer [optional, default is <code>floor(size(M, 1)/2)</code>] targeted number of Givens rotations per factor of V. Must be an integer between 1 to <code>floor(size(M, 1)/2)</code>.
%> @param 'relerr', bool [optional, default is true] the type of error used as stopping criterion. (true) for the relative error norm(V*D*V'-M, 'fro')/norm(M, 'fro'), (false) for the absolute error norm(V*D*V'-M, 'fro').
%> @param 'enable_large_Faust', bool [optional, default is false] if true, it allows to compute a transform that doesn't worth it regarding its complexity relatively to the matrix M. Otherwise, by default, an exception is raised before the algorithm starts.
%> @param 'verbosity', integer [optional] verbosity level. The greater the value the more info is displayed. It can be helpful to understand for example why the algorithm stopped before reaching the tol error or the number of Givens (nGivens).
%>
%>
%> @retval [V,D]
%> - V the Faust object representing the approximate eigenvector transform. The column <code>V(:, i)</code> is the eigenvector corresponding to the eigenvalue <code>D(i,i)</code>.
%> - D the sparse (real) diagonal matrix of the approximate eigenvalues (by default in ascending order along the diagonal).
%>
%> @note
%>	- When  ‘nGivens’ and ‘tol’ are used simultaneously, the number of Givens rotations in V may be smaller than specified by ‘nGivens’ if the error criterion is met first, and the achieved error may be larger than specified if ‘nGivens’ is reached first during the iterations of the truncated Jacobi algorithm.
%> @note
%> - When nGivens_per_fac > 1, all factors have exactly nGivens_per_fac except the leftmost one which may have fewer if the total number of Givens rotations is not a multiple of nGivens_per_fact
%>
%>
%> @b Example
%> @code
%> import matfaust.fact.eigtj
%>
%> % get a Laplacian to diagonalize
%> load('Laplacian_256_community.mat')
%> % do it
%> [Uhat, Dhat] = eigtj(Lap, 'nGivens', size(Lap,1)*100, 'nGivens_per_fac', size(Lap, 1)/2, 'verbosity', 2, 'enable_large_Faust', true)
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
	nGivens_per_fac = floor(size(M,1)/2); % default value
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
	enable_large_Faust = false;
	argc = length(varargin);
	order = 1; % ascending order
	err_period = 100;
	if(argc > 0)
		for i=1:2:argc
			switch(varargin{i})
				case 'enable_large_Faust'
					if(argc == i || ~ islogical(varargin{i+1}))
						error('enable_large_Faust keyword argument is not followed by a logical')
					else
						enable_large_Faust = varargin{i+1};
					end
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
				case 'err_period'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('err_period keyword argument is not followed by a number')
					else
						err_period = floor(real(varargin{i+1}));
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
	if(strcmp(class(M), 'single'))
		[core_obj, D] = mexfgftgivensRealFloat(M, nGivens, nGivens_per_fac, verbosity, tol, relerr, order, enable_large_Faust, err_period);
		D = sparse(diag(real(double(D))));
		V = Faust(core_obj, isreal(M), 'cpu', 'float');
	else
		[core_obj, D] = mexfgftgivensReal(M, nGivens, nGivens_per_fac, verbosity, tol, relerr, order, enable_large_Faust, err_period);
		D = sparse(diag(real(D)));
		V = Faust(core_obj, isreal(M));
	end
end
