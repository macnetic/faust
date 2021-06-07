% experimental block start
%==========================================================================
%> @brief Factorizes the matrix M according to the butterfly support.
%==========================================================================
function F = butterfly(M)
	import matfaust.Faust
	if(isreal(M))
		core_obj = mexButterflyReal(M);
	else
		core_obj = mexButterflyCplx(M);
	end
	F = Faust(core_obj, isreal(M));
end
% experimental block end
