%====================================================================
%> @brief Performs an approximate singular value decomposition and returns the left and
%> right singular vectors as Faust transforms.
%>
%> @note this function is based on fact.eigtj which relies on the truncated Jacobi algorithm, hence the 'tj' in the name. See below the examples for further details on how svdtj is defined using eigtj.
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
%> @param M: a real or complex, dense or sparse matrix. The class(M) value can be double or single (only if isreal(M) is true).
%> @param nGivens, int|matrix: defines the number of Givens rotations
%> that will be used at most to compute U and V.
%> If it is an integer, it will apply both to U and V.
%> If it is a tuple of two integers as nGivens = [JU, JV], JU
%> will be the limit number of rotations for U and JV the same for V.
%> nGivens argument is optional if tol is set but becomes mandatory otherwise.
%> @param 'tol', number: this is the error target on U' S  V against M.
%> if error <= tol, the algorithm stops. See relerr below for the error formula.
%> 	This argument is optional if nGivens is set, otherwise it becomes mandatory.
%> @param 'relerr', bool: defines the type of error used to evaluate the stopping criterion defined by tol. true for a relative error (\f$\| U' S V - M\|_F \over \| M \|_F\f$) otherwise this is an absolute error (\f$\| U' S V - M\|_F \f$).
%> @param 'nGivens_per_fac',int: this argument is the number of Givens
%> rotations to set at most by factor of U and V.
%> If this is an integer it will be the same number of U and V.
%> Otherwise, if it is a tuple of integers [tU, tV], tU will be the number
%> of Givens rotations per factor for U and tV the same for V.
%> By default, this parameter is maximized for U and V,
%> i.e. tU = size(M, 1) / 2, tV = size(M, 2) / 2.
%> @param 'enable_large_Faust',bool see fact.eigtj
%> @param 'err_period', int:  it defines the period, in number of factors of U
%> or V, the error is compared to tol (reducing the period spares
%> some factors but increases slightly the computational cost because the error
%> is computed more often).
%>
%>
%> \note In order to speed up the error computation of the approximate Faust, the
%> algorithm follows a two-times strategy:
%> 1. In the first iterations a rough error but less costly error is computed \f$\| M \|^2_F - \| S \|^2_F\f$
%> 2. In the latest iterations the precise error (relative or absolute) indicated above
%> (see relerr argument) is computed to reach the targeted accuracy. Don't forget that
%> the err_period argument determines the frequency at which the error is
%> calculated, you might decrease its value in order to obtain an error closer to
%> what you asked. In last resort you can disable the two-times strategy by
%> setting the environment related variable like this (defaulty the value is '0', which means
%> that the two-times strategy is used).
%> @code
%> setenv('SVDTJ_ALL_TRUE_ERR', '1')
%> @endcode
%>
%> @retval [U,S,V]: such that U*S*V' is the approximate of M with:
%>      - S: (sparse real diagonal matrix) the singular values in descending order.
%>      - U, V: (Faust objects) orthonormal transforms.
%>
%> @b Example
%> @code
%> >> import matfaust.fact.svdtj
%> >> rng(42) % for reproducibility
%> >> M = rand(16, 32);
%> >> % Factoring by specifying the numebr of Givens rotations
%> >> [U1, S1, V1] = svdtj(M, 'nGivens', 4096, 'enable_large_Faust', true);
%> >> % verify the approximate is accurate
%> >> norm(U1 * S1 * V1' - M) / norm(M)
%>
%> ans =
%>
%>    3.1328e-15
%> >> % Specifying a different number of rotations for U and V
%> >> % Because U is smaller it should need less rotations
%> >> [U2, S2, V2] = svdtj(M, 'nGivens', [2400, 3200], 'enable_large_Faust', true);
%> >> norm(U2 * S2 * V2' - M) / norm(M)
%>
%> ans =
%>
%>    3.1141e-15
%> >> % Factoring according to an approximate accuracy target
%> >> [U3, S3, V3] = svdtj(M, 'tol', 1e-12, 'enable_large_Faust', false);
%> Warning: eigtj stopped on U approximate which is too long to be worth it, set enable_large_Faust to true if you want to continue anyway...
%> Warning: eigtj stopped on V approximate which is too long to be worth it, set enable_large_Faust to true if you want to continue anyway...
%> >> norm(U3 * S3 * V3' - M) / norm(M)
%>
%> ans =
%>
%>      1.5950e-13
%>
%> >> % try with an absolute toelrance (the previous one was relative to M norm)
%> >> [U4, S4, V4] = svdtj(M, 'tol', 1e-6, 'relerr', false, 'enable_large_Faust', true);
%> >> % verify the absolute error is lower than 1e-6
%> >> norm(U4 * S4 * V4' - M)
%>
%> ans =
%>
%>    3.6203e-14
%> >> % try a less accurate approximate to get less factors
%> >> [U5, S5, V5] = svdtj(M, 'nGivens', [256, 512], 'tol', 1e-1, 'enable_large_Faust', false);
%> >> norm(U5 * S5 * V5' - M) / norm(M)
%>
%> ans =
%>
%>     0.0500
%> >> %%% Now let's see the lengths of the different U, V Fausts
%> >> length(V1) % it should be 4096 / nGivens_per_fac, which is (size(M, 2) / 2) = 256
%>
%> ans =
%>
%>    256
%>
%> >> length(U1) % it should be 4096 / nGivens_per_fac, which is (size(M, 1) / 2) = 512
%>
%> ans =
%>
%>    100
%>
%> >> % but it is not, svdtj stopped automatically on U1 because its error stopped enhancing
%> >> % (it can be verified with: 'verbosity', 1)
%> >> [length(U2), length(V2)]
%>
%> ans =
%>
%>    100   200
%>

%> >> [length(U3), length(V3)]
%>
%> ans =
%>
%>     64   256
%>
%> >> [length(U4), length(V4)]
%>
%> ans =
%>
%>    100   200
%>
%> >> % not surprisingly U5 and V5 use the smallest number of factors (nGivens and tol were the smallest)
%> >> [length(U5), length(V5)]
%>
%> ans =
%>
%>     32    32
%>
%> @endcode
%>
%> Explanations:
%>
%>If we call svdtj on the matrix M, it makes two internal calls to eigtj.
%>
%> In Matlab it would be:
%> 1.  [D1, W1] = eigtj(M*M, varargin)
%> 2.  [D2, W2] = eigtj(M'*M, varargin)
%>
%> It gives the following equalities (ignoring the fact that eigtj computes approximations):
%> \f[
%>     W_1 D_1 W_1^* = M M^*
%> \f]
%> \f[
%>    W_2 D_2 W_2^* = M^* M
%> \f]
%> But because of the SVD \f$ M = USV^* \f$ we also have:
%> \f[MM^* = U S V^* V S U^* = U S^2 U^* = W_1 D_1 W_1^*\f]
%> \f[M^* M = V S U^* U S V^* = V S^2 V^* = W_2 D_2  W_2^*\f]
%> It allows to identify the left singular vectors of M to W1,
%> and likewise the right singular vectors to W2.
%>
%> To compute a consistent approximation of S we observe that U and V are orthogonal/unitary hence \f$ S  = U^* M V \f$ so we ignore the off-diagonal coefficients of the approximation and take \f$ S = diag(U^* M V)  \approx diag(W_1^* M W_2)\f$
%>
%> The last step performed by svdtj() is to sort the singular values of S in descending order and build a signed permutation matrix to order the left singular vectors of W1 accordingly. The -1 elements of the signed permutation matrix allow to change the sign of each negative values of S by reporting it on the corresponding left singular vector (\f$ \sigma v_i = (-\sigma_i) (-v_i )\f$).<br/>
%> To sum up W1 is replaced by W1 P and W2 by W2 abs(P) (because W2 also needs to be ordered), with P the signed permutation resulting of the descending sort of S. The resulting transforms/Fausts W1 and W2 are returned by svdtj along with the ordered S. Note that the permutation factor P (resp. abs(P)) is fused with the rightmost factor of the Faust object W1 (resp. W2).
%>
%> <p> @b See @b also fact.eigtj, fact.pinvtj
%>
%>
%====================================================================
function [U,S,V] = svdtj(M, varargin)
%	TODO: factorize argument parsing code with eigetj
	import matfaust.Faust
	% default values
	verbosity = 0;
	bad_arg_err = 'bad number of arguments.';
	nGivens_per_fac = [floor(size(M, 1) / 2), floor(size(M, 2) / 2)];
	nGivens = [0, 0];
	tol = 0;
	relerr = true;
	argc = length(varargin);
	enable_large_Faust = false;
	order = -1; % descending order
	err_period = 100;
	if(argc > 0)
		for i=1:2:argc
			switch(varargin{i})
				case 'nGivens'
					nGivens = varargin{i+1};
					if size(nGivens) == [1, 1]
						nGivens = [nGivens, nGivens];
					end
					if(argc == i || ~ isreal(nGivens) || any(nGivens-floor(nGivens) ~= [0, 0]))
						error('nGivens keyword arg. is not followed by an integer or a pair of integers')
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
					nGivens_per_fac = varargin{i+1};
					if size(nGivens_per_fac) == [1, 1]
						nGivens_per_fac = [nGivens_per_fac, nGivens_per_fac];
					end
					if(argc == i || ~ isreal(nGivens_per_fac) || any(nGivens_per_fac-floor(nGivens_per_fac) ~= [0, 0]))
						error('nGivens_per_fac must be followed by a positive integer or a pair of positive integers.')
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
				case 'enable_large_Faust'
					if(argc == i || ~ islogical(varargin{i+1}))
						error('enable_large_Faust keyword argument is not followed by a logical')
					else
						enable_large_Faust = varargin{i+1};
					end
				case 'err_period'
					if(argc == i || ~ isscalar(varargin{i+1}))
						error('err_period keyword argument is not followed by a number')
					else
						err_period = floor(real(varargin{i+1}));
					end
				otherwise
					varargin{i}
					error([ 'unrecognized argument'])
			end
		end
	end
	if(all(nGivens == 0) && tol == 0)
		error('Either nGivens or tol must be greater than zero.')
	end
	if(nGivens ~= [0, 0])
		nGivens_per_fac = min(nGivens_per_fac, nGivens);
	end
	mex_args = {M, nGivens(1), nGivens(2), nGivens_per_fac(1), nGivens_per_fac(2), verbosity, tol, relerr, order, enable_large_Faust, err_period};
	if isreal(M)
		if(strcmp(class(M), 'single'))
			[core_obj1, S, core_obj2] = mexsvdtjfloat(mex_args{:});
			S = spdiags(double(S), 0, size(M, 1), size(M, 2)); % matlab doesn't support single precision sparse matrix
			U = Faust(core_obj1, isreal(M), 'cpu', 'float');
			V = Faust(core_obj2, isreal(M), 'cpu', 'float');
		else
			[core_obj1, S, core_obj2] = mexsvdtjdouble(mex_args{:});
			S = spdiags(double(S), 0, size(M, 1), size(M, 2));
			U = Faust(core_obj1, isreal(M));
			V = Faust(core_obj2, isreal(M));
		end
	else % M is complex
			[core_obj1, S, core_obj2] = mexsvdtjdouble(mex_args{:});
			S = spdiags(S, 0, size(M, 1), size(M, 2));
			U = Faust(core_obj1, isreal(M), 'cpu', 'double', true);
			V = Faust(core_obj2, isreal(M), 'cpu', 'double', true);
	end
end
