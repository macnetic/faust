% experimental block start
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
	[W1,D1] = matfaust.fact.eigtj(M*M', maxiter, varargin{:});
	[W2,D2] = matfaust.fact.eigtj(M'*M, maxiter, varargin{:});
	S = diag(W1'*M*W2);
	[~,I] = sort(abs(S), 'descend');
	S = sparse(diag(S(I)));
	sign_S = sign(S);
	S = S*sign_S;
	Id = eye(size(S));
	U = W1(:,1:size(Id,1))*matfaust.Faust({Id(:,I),sign_S});
	V = W2(:,1:size(Id,1))*matfaust.Faust(Id(:,I));
end
% experimental block end
