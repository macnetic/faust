function check_fact_mat(funcname, M)
	if(~ ismatrix(M) || isscalar(M))
		error([funcname,'() 1st argument (M) must be a matrix.'])
	end
	if(~ isnumeric(M))
		error([funcname, '() 1st argument (M) must be real or complex.'])
	end
%			if(~ isreal(M))
%				error([funcname, '() doesn''t yet support complex matrix factorization.'])
%			end
end

