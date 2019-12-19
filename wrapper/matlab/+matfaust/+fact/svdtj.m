%====================================================================
%> @brief Performs an approximate singular value decomposition and returns the left and
%> right singular vectors as Faust transforms.
%>
%> @note this function is based on fact.eigtj which relies on the truncated Jacobi algorithm, hence the 'tj' in the name. See below the example for further details on how svdtj is defined using eigtj.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b Primary examples of calls include:<br/>
%> &nbsp;&nbsp;&nbsp; @b <b>[U,S,V] = svdtj(M, ‘nGivens’,n) </b> outputs U,V as Faust objects made of n elementary Givens rotations grouped into factors. By default each factor gathers (up to) m = floor(size(M,1)/2) such rotations. By default vector S contains the approximate singular values in descending order.<br/>
%> &nbsp;&nbsp;&nbsp; @b    <b>[U,S,V] = svdtj(M,’tol’,0.01)</b> same as above with n determined adaptively by a relative approximation error 0.01<br/>
%> &nbsp;&nbsp;&nbsp; @b       <b>[U,S,V] = svdtj(M,’tol’,0.01,’relerr’,false)</b> same as above with an absolute approximation error<br/>
%> &nbsp;&nbsp;&nbsp; @b       <b>[U,S,V] = svdtj(M,’nGivens’,n,’tol’,0.01)</b> same as above with a number of elementary Givens bounded by n even if the targeted approximation error is not achieved<br/>
%> &nbsp;&nbsp;&nbsp; @b       <b>[U,S,V] = svdtj(M,’nGivens’,n,’tol’,0.01,’nGivens_per_fac’,t) </b>same as above with (up to) t Givens rotations per factor<br/>
%>
%>
%> @param M: a real or complex, dense or sparse matrix.
%> @param nGivens, integer: see fact.eigtj
%> @param 'tol', number see fact.eigtj (NB: as described below, the error tolerance is not exactly for the approximate SVD but for the subsequent eigtj calls).
%> @param 'relerr', bool see fact.eigtj
%> @param 'nGivens_per_fac',integer see fact.eigtj
%>
%> @retval [U,S,V]: such that U*S*V' is the approximate of M with:
%>      - S: (sparse real diagonal matrix) the singular values in descendant order.
%>      - U, V: (Faust objects) unitary transforms.
%>
%> @Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.fact.svdtj
%> >> M = rand(128,128)
%> >> [U,S,V] = svdtj(M, 'nGivens', 1024, 'nGivens_per_fac', 64)
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
%>        To sum up W1 is replaced by W1 P and W2 by W2 abs(P) (because W2 also needs to be ordered), with P the signed permutation resulting of the descending sort of S. The resulting transforms/Fausts W1 and W2 are returned by svdtj along with the ordered S. Note that the permutation factor P (resp. abs(P)) is fused with the rightmost factor of the Faust object W1 (resp. W2).
%>
%> <p> @b See @b also fact.eigtj
%>
%>
%====================================================================
function [U,S,V] = svdtj(M, varargin)
%	[W1,D1] = matfaust.fact.eigtj(M*M', nGivens, varargin{:});
%	[W2,D2] = matfaust.fact.eigtj(M'*M, nGivens, varargin{:});
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
	verbosity = 0; % default value
%	if(~ ismatrix(M) || ~ isreal(M))
%		error('M must be a real matrix.')
%	end
	if(size(M,1) ~= size(M,2))
		error('M must be square')
	end
	bad_arg_err = 'bad number of arguments.';
	nGivens_per_fac = floor(size(M,1)/2); % default value
	nGivens = 0;
	tol = 0;
	relerr = true;
	verbosity = 0;
	argc = length(varargin);
	order = -1; % descending order
	if(argc > 0)
		for i=1:argc
			switch(varargin{i})
				case 'nGivens'
					nGivens = varargin{i+1};
					if(argc == i || ~ isnumeric(nGivens) || nGivens-floor(nGivens) > 0 || nGivens <= 0 )
						error('nGivens keyword arg. is not followed by an integer')
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
	if(nGivens == 0 && tol == 0)
		error('Either nGivens or tol must be greater than zero.')
	end
	if(nGivens > 0)
		nGivens_per_fac = min(nGivens_per_fac, nGivens);
	end
	[core_obj1, S, core_obj2] = mexsvdtjReal(M, nGivens, nGivens_per_fac, verbosity, tol, relerr, order);
	S = sparse(diag(real(S)));
	U = Faust(core_obj1, isreal(M));
	V = Faust(core_obj2, isreal(M));
end
