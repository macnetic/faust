%====================================================================
%> @brief Performs a singular value decomposition and returns the left and
%> right singular vectors as Faust transforms.
%>
%> @note this function is based on fact.eigtj.
%>
%> @param M: a real matrix.
%> @param maxiter: see fact.eigtj
%> @param 'nGivens_per_fac',integer see fact.eigtj
%> @param nGivens_per_fac: see fact.eigtj
%> @param 'tol', number see fact.eigtj
%> @param 'relerr', bool see fact.eigtj
%> @param 'verbosity', integer see fact.eigtj
%>
%> @retval [U,S,V]: U*full(S)*V' being the approximation of M.
%>      - S: (sparse diagonal matrix) S the singular values in
%>		descendant order.
%>      - U: (Faust object) U the left-singular transform.
%>      - V: (Faust object) V the right-singular transform.
%>
%> @Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.fact.svdtj
%> >> M = rand(128,128)
%> >> [U,S,V] = svdtj(M,1024,'nGivens_per_fac', 64)
%> @endcode
%>
%====================================================================
function [U,S,V] = svdtj(M, maxiter, varargin)
%	[W1,D1] = matfaust.fact.eigtj(M*M', maxiter, varargin{:});
%	[W2,D2] = matfaust.fact.eigtj(M'*M, maxiter, varargin{:});
%	S = diag(W1'*M*W2);
%	[~,I] = sort(abs(S), 'descend');
%	S = sparse(diag(S(I)));
%	sign_S = sign(S);
%	S = S*sign_S;
%	Id = eye(size(S));
%	U = W1(:,1:size(Id,1))*matfaust.Faust({Id(:,I),sign_S});
%	V = W2(:,1:size(Id,1))*matfaust.Faust(Id(:,I));
%	TODO: factorize argument parsing code with fgft_givens
	import matfaust.Faust
	nGivens_per_fac = 1; % default value
	verbosity = 0; % default value
%	if(~ ismatrix(M) || ~ isreal(M))
%		error('M must be a real matrix.')
%	end
	if(size(M,1) ~= size(M,2))
		error('M must be square')
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
	[core_obj1, S, core_obj2] = mexsvdtjReal(M, maxiter, nGivens_per_fac, verbosity, tol, relerr, order);
	S = sparse(diag(real(S)));
	U = Faust(core_obj1, isreal(M));
	V = Faust(core_obj2, isreal(M));
end
