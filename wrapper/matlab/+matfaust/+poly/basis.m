%======================================================================
%> @brief  TODO
%>
%=======================================================================
function F = basis(L, K, basis_name, varargin)
	% TODO: T0 in varargin
	if(~ ismatrix(L) || ~ issparse(L) || size(L, 1) ~= size(L, 2))
		error("L must be a square matrix.")
	end
	
	is_real = isreal(L)
	if(is_real)
		core_obj = mexPolyReal('chebyshev', L, K);
	else
		core_obj = mexPolyCplx('chebyshev', L, K);
	end

	P = matfaust.Faust(core_obj, is_real)
	disp(P)
end
