% ===================================
%> @brief This class defines the same parameterization as ParamsHierarchicalRectMat except that there is no constraint on the residual factors (cf.  matfaust.proj.proj_id).
%>
%> <p>@b See @b  matfaust.fact.hierarchical</p>
% ===================================
classdef ParamsHierarchicalRectMatNoResCons < matfaust.factparams.ParamsHierarchicalRectMat
    methods

        %===================================
        %> @brief Class constructor.
        %>
        %> @param m cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param n cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param j cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param k cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param s cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param rho (varargin{1}) cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %> @param P (varargin{2}) cf. ParamsHierarchicalRectMat.ParamsHierarchicalRectMat
        %>
        %===================================
        function p = ParamsHierarchicalRectMatNoResCons(m, n, j, k, s, varargin)
            %%
            import matfaust.proj.proj_id
            p = p@matfaust.factparams.ParamsHierarchicalRectMat(m, n, j, ...
                k, s, varargin{:});
            % Remove  all constraints on residual factors except the last one
            n_cons = length(p.constraints);
            for i=1:n_cons/2-1
                p.constraints{i} = proj_id([m, m]).constraint;
            end
        end
    end
    methods(Static)

        %===================================
        %> @brief Static member function to create a ParamsHierarchicalRectMatNoResCons instance by a simplified parameterization expression.
        %>
        %>
        %> @b Example
        %> @code
        %> >> import matfaust.factparams.ParamsHierarchicalRectMatNoResCons
        %> >> num_facts = 9;
        %> >> k = 10;
        %> >> s = 8;
        %> >> p = ParamsHierarchicalRectMatNoResCons.createParams(rand(256, 1024), {'rectmat_simple', num_facts, k, s});
        %> @endcode
        %===================================
        function sp = createParams(M, p)
            %%
            import matfaust.factparams.*
            ParamsHierarchicalRectMat.createParams(M, p); % just to verify
            % input arguments
            sp = ParamsHierarchicalRectMatNoResCons(size(M,1), size(M,2), p{2:end});
        end
    end
end
