%====================================================================
%> @brief Performs a singular value decomposition and returns the left and
%> right singular vectors as Faust transforms.
%>
%> @note this function is based on fact.eigtj.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b [U,S,V] = svdtj(M, varargin). See below for further details on how svdtj is defined using eigtj.<br/>
%%>
%> @param M: a real matrix.
%> @param maxiter: see fact.eigtj
%> @param 'nGivens_per_fac',integer see fact.eigtj
%> @param nGivens_per_fac: see fact.eigtj
%> @param 'tol', number see fact.eigtj (the error tolerance is not exactly for the svdtj but for the subsequent eigtj calls).
%> @param 'relerr', bool see fact.eigtj
%> @param 'verbosity', integer see fact.eigtj
%>
%> @retval [U,S,V]: U*S*V' being the approximation of M.
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
%> >> [U,S,V] = svdtj(M, 'maxiter', 1024,'nGivens_per_fac', 64)
%> @endcode
%>If we call svdtj on the matrix M, it makes two internal calls to eigtj.
%>
%>        In Matlab it would be:
%>        1.  [D1, W1] = eigtj(M*M, varargin)
%>        2.  [D2, W2] = eigtj(M'*M, varargin)
%>
%>        It gives the following equalities (ignoring the fact that eigtj computes approximations):
%>        \f[
%>            W_1 D_1 W_1^* = M M^*
%>        \f]
%>        \f[
%>                W_2 D_2 W_2^* = M^* M
%>        \f]
%>        But because of the SVD \f$ M = USV^* \f$ we also have:
%>        \f[MM^* = U S V^* V S U^* = U S^2 U^* = W_1 D_1 W_1^*\f]
%>        \f[M^* M = V S U^* U S V^* = V S^2 V^* = W_2 D_2  W_2^*\f]
%>        It allows to identify the left singular vectors of M to W1,
%>        and likewise the right singular vectors to W2.
%>
%>        To compute a consistent approximation of S we observe that U and V are orthogonal/unitary hence \f$ S  = U^* M V \f$ so we ignore the off-diagonal coefficients of the approximation and take \f$ S = diag(U^* M V)  \approx diag(W_1^* M W_2)\f$
%>
%>        The last step performed by svdtj() is to sort the singular values of S in descending order and build a signed permutation matrix to order the left singular vectors of W1 accordingly. The -1 elements of the signed permutation matrix allow to change the sign of each negative values of S by reporting it on the corresponding left singular vector (\f$ \sigma v_i = (-\sigma_i) (-v_i )\f$).<br/>
%>        To sum up W1 is replaced by W1 P and W2 by W2 abs(P) (because W2 also needs to be ordered), with P the signed permutation resulting of the descending sort of S. That new transforms/Fausts W1 and W2 are returned by svdtj along with the ordered S. Note that the permutation factor is not append to the transform W1 (or W2) but multiplied directly to the last factor of W1 (or W2).
%>
%> <p> @b See @b also fact.eigtj
%>
%>
%====================================================================
function [U,S,V] = svdtj(M, varargin)
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
	bad_arg_err = 'bad number of arguments.';
	maxiter = 0;
	tol = 0;
	relerr = true;
	verbosity = 0;
	argc = length(varargin);
	order = 1; % ascending order
	if(argc > 0)
		for i=1:argc
			switch(varargin{i})
				case 'maxiter'
					maxiter = varargin{i+1};
					if(argc == i || ~ isnumeric(maxiter) || maxiter-floor(maxiter) > 0 || maxiter <= 0 )
						error('maxiter keyword arg. is not followed by an integer')
					end
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
	if(maxiter == 0 && tol == 0)
		error('Either maxiter or tol must be greater than zero.')
	end
	if(maxiter > 0)
		nGivens_per_fac = min(nGivens_per_fac, maxiter);
	end
	[core_obj1, S, core_obj2] = mexsvdtjReal(M, maxiter, nGivens_per_fac, verbosity, tol, relerr, order);
	S = sparse(diag(real(S)));
	U = Faust(core_obj1, isreal(M));
	V = Faust(core_obj2, isreal(M));
end
