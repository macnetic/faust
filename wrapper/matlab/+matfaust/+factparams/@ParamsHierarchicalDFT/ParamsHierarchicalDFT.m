% =========================================================
%> @brief The simplified parameterization class for factorizing a DFT matrix using the hierarchical factorization algorithm.
%>
%> <p>@b See @b also matfaust.fact.hierarchical</p>
% =========================================================
classdef ParamsHierarchicalDFT < matfaust.factparams.ParamsHierarchical
	properties (SetAccess = public)
		supports
	end

	methods
		function p = ParamsHierarchicalDFT(n)
			import matfaust.factparams.*
			n = floor(n);
			d = 2^n;
			stop_crit = StoppingCriterion(8);
			supports = support_DFT(n);
			fac_cons = cell(1,n);
			for i=1:n
				fac_cons{1, i} = ConstraintMat('supp', complex(supports{i}));
			end
			res_cons = cell(1,n);
			for j=1:n-1
				supp = complex(eye(size(supports{1})));
				for i=j+1:n+1
					supp = supp*supports{i};
				end
				supp(nonzeros(supp)) = 1;
				res_cons{1,j} = ConstraintMat('supp', complex(supp));
			end
			res_cons{1,n} = ConstraintMat('supp', complex(supports{n+1}));
			p = p@matfaust.factparams.ParamsHierarchical(fac_cons, res_cons, stop_crit,...
			stop_crit, 'is_update_way_R2L', true);
			p.supports = supports;
		end

	end

	methods(Static)

		function sp = createParams(M, p)
			import matfaust.factparams.ParamsHierarchicalDFT
			pot = log2(size(M,1));
			if(size(M,1) ~= size(M,2) || pot-floor(pot) > 0)
				error('M must be a square matrix of order a power of two.')
			end
			sp = ParamsHierarchicalDFT(pot);
		end


	end
end

function v=bit_rev_permu(n)
	n = floor(n);
	if n == 0
		v = [1];
	end
	size = bitshift(1, n);
	v = 0:size-1;
	lower_mask = 1;
	upper_mask = bitshift(1, n-1);
	shift = 0;
	while lower_mask < upper_mask
		for i=1:size
			lobit = bitshift(bitand(v(i), lower_mask), -shift);
			hibit = bitshift(bitand(v(i), upper_mask), -n+shift+1);
			if lobit > hibit
				v(i) = bitxor(v(i), lower_mask);
				v(i) = bitor(v(i), upper_mask);
			elseif lobit < hibit
				v(i) = bitor(v(i), lower_mask);
				v(i) = bitxor(v(i), upper_mask);
			end
		end
		lower_mask = bitshift(lower_mask, 1);
		upper_mask = bitshift(upper_mask, -1);
		shift = shift + 1;
	end
	v = v + 1;
end

function supports=support_DFT(n)
	size = 2^n;
	supports = cell(1,n+1);
	for i=0:n-1
		supp_bf = kron(ones(2,2), eye(2 ^ ((n-i)-1)));
		supports{i+1} = complex(kron(eye(2^i), supp_bf));
	end
	% bit-reversal permutation
	row_ids = 1:size;
	col_ids = bit_rev_permu(n);
	br_permu = complex(zeros(size,size));
	br_permu(row_ids, col_ids) = 1;
	supports{n+1} = br_permu;
end
