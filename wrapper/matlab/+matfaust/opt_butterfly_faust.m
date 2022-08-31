%=================================
%> See matfaust.Faust.opt_butterfly alias.
%>
%=================================
function G = opt_butterfly_faust(F)
	import matfaust.Faust;
	G = Faust.opt_butterfly(F)
%	if ~ matfaust.isFaust(F)
%		error('F must be a Faust.')
%	end
%	 G = matfaust.Faust(F, call_mex(F, 'opt_butterfly', F.matrix.objectHandle));
%	func_name = 'opt_butterfly';
%	if (strcmp(F.dev, 'cpu'))
%		if(F.is_real)
%			if(strcmp(F.dtype, 'double'))
%				G = matfaust.Faust(F, mexFaustReal(func_name, F.matrix.objectHandle));
%			else % float
%				G = matfaust.Faust(F, mexFaustRealFloat(func_name, F.matrix.objectHandle));
%			end
%		else
%			G = matfaust.Faust(F, mexFaustCplx(func_name, F.matrix.objectHandle));
%		end
%	elseif(startsWith(F.dev, 'gpu'))
%		if(F.is_real)
%			if(strcmp(F.dtype, 'double'))
%				G = matfaust.Faust(F, mexFaustGPUReal(func_name, F.matrix.objectHandle));
%			else % float
%				G = matfaust.Faust(F, mexFaustGPURealFloat(func_name, F.matrix.objectHandle));
%			end
%		else
%			G = matfaust.Faust(F, mexFaustGPUCplx(func_name, F.matrix.objectHandle));
%		end
%	end
end
