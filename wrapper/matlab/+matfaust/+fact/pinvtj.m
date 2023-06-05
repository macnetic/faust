%====================================================================
%> @brief Computes the pseudoinverse of M using svdtj.
%>
%> @param M: the matrix to compute the pseudoinverse.
%> @param nGivens, int|matrix: see fact.svdtj
%> @param 'tol', number: see fact.svdtj
%> @param 'relerr', bool: see fact.svdtj
%> @param 'nGivens_per_fac',int: see fact.svdtj
%> @param 'enable_large_Faust',bool: see fact.svdtj
%> @param 'err_period', int: see fact.svdtj
%>
%>
%> @retval [V,Sp,Uh]: such that V*Sp*Uh is the approximate of M^+ with:
%>      - Sp: (sparse real diagonal matrix) the pseudoinverse of the matrix of the singular values.
%>      - V, Uh: (Faust objects) orthonormal transforms.
%>
%> @b Example
%> @code
%> >> import matfaust.fact.pinvtj
%> >> rng(42) % for reproducibility
%> >> M = rand(128, 64);
%> >> [V, Sp, Uh] = pinvtj(M, 'tol', 1.5e-2);
%> >> norm(V * Sp * Uh * M - eye(64)) / norm(eye(64))
%>
%> ans =
%>
%>      0.0036
%>
%> @endcode
%====================================================================
function [V,Sp,Uh] = pinvtj(M, varargin)
    %%
    import matfaust.fact.svdtj
    setenv('PINVTJ_ERR', '1');
    [U, S, V] = svdtj(M, varargin{:});
    setenv('PINVTJ_ERR', '0');
    Sp = 1 ./ nonzeros(S);
    Sp(Sp == Inf) = 0;
    Sp = spdiags(Sp, 0, size(V, 1), size(U, 1));
    Uh = U';
end
