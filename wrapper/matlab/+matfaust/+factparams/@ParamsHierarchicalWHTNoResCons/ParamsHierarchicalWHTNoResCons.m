% =================================================
%> The simplified parameterization class for factorizing a Hadamard matrix with the hierarchical factorization algorithm.

%> This type of parameters is typically used for a Hadamard matrix factorization.
%> This is a variant of ParamsHierarchicalWHT. Here the intermediate residual
%> factors are not constrained at all, the other factors are constrained with
%> matfaust.proj.skperm.
%>
%>
%> <b/> See also matfaust.fact.hierarchical, matfaust.demo.hadamard
% =================================================
classdef ParamsHierarchicalWHTNoResCons < matfaust.factparams.ParamsHierarchicalNoResCons
    methods
        %=================================================
        %> @param n: the number of output factors (the input matrix to factorize must
        %> be of size [2^n, 2^n]) .
        %>
        %> @b Example
        %> @code
        %> >> import matfaust.wht
        %> >> import matfaust.factparams.*
        %> >> import matfaust.fact.hierarchical
        %> >> H = full(wht(32));
        %> >> p = ParamsHierarchicalWHTNoResCons(5);
        %> >> F = hierarchical(H, p);
        %> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 1/4
        %> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 2/4
        %> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 3/4
        %> Faust::HierarchicalFact<FPP,DEVICE,FPP2>::compute_facts : factorization 4/4
        %>
        %> >>
        %> @endcode
        %=================================================
        function p = ParamsHierarchicalWHTNoResCons(n, varargin)
            import matfaust.factparams.*
            proj_name = 'skperm'; % because splincol doesn't work well
            cons_name = ConstraintName.str2name_int(proj_name);
            d = 2^ceil(n);
            stop_crit = StoppingCriterion(30);
            cons = {};
            for i=1:n-1
                cons = {cons{:} ConstraintInt(ConstraintName(cons_name), d, d, 2)};
            end
            cons = {cons{:} ConstraintInt(ConstraintName(cons_name), d, d, ceil(d/2^(n-1)))};
            p@matfaust.factparams.ParamsHierarchicalNoResCons(cons, stop_crit, ...
                stop_crit, 'is_update_way_R2L', true, 'packing_RL', false, ...
                varargin{:});
        end
    end
    methods(Static)
        function sp = createParams(M, p)
            import matfaust.factparams.ParamsHierarchicalWHTNoResCons
            pot = log2(size(M,1));
            if(size(M,1) ~= size(M,2) || pot-floor(pot) > 0)
                error('M must be a square matrix of order a power of two.')
            end
            sp = ParamsHierarchicalWHTNoResCons(pot);
        end
    end
end
