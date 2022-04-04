classdef ParamsFactFactory
    properties(Constant)
        SIMPLIFIED_PARAMS = {
        {'squaremat', 'hadamard'},
        {'rectmat', 'meg'},
        {'dft'},
        {'hadamard_simple', 'hadamard_no_rescons'},
        {'rectmat_simple', 'meg_simple'}
        };
        SQRTMAT_ID = 1;
        RECTMAT_ID = 2;
        DFTMAT_ID = 3;
        WHT_SIMPLE_ID = 4
        RECTMAT_SIMPLE_ID = 5;
    end
    methods(Static)
        function sp = createParams(M, p)
            import matfaust.factparams.*
            import matfaust.factparams.ParamsFactFactory.*
            param_id = -1;
            if(~ is_a_valid_simplification(p))
                error('Invalid p to represent a simplified parametrization')
            end
            if(ischar(p))
                param_id = lower(p);
            else % p{1} is a cell with a char array at first pos
                param_id = lower(p{1});
            end
            if(any(strcmp(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.SQRTMAT_ID})))
                sp = ParamsHierarchicalWHT.createParams(M, p);
            elseif(any(strcmp(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.WHT_SIMPLE_ID})))
                sp = ParamsHierarchicalWHTNoResCons.createParams(M, p);
            elseif(any(strcmp(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.RECTMAT_ID})))
                sp = ParamsHierarchicalRectMat.createParams(M, p);
            elseif(any(strcmp(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.RECTMAT_SIMPLE_ID})))
                sp = ParamsHierarchicalRectMatNoResCons.createParams(M, p);
            elseif(any(strcmp(param_id, ParamsFactFactory.SIMPLIFIED_PARAMS{ParamsFactFactory.DFTMAT_ID})))
                sp = ParamsHierarchicalDFT.createParams(M, p);
            else
                error('p is not a known simplified parametrization')
            end
        end
        function bool = is_a_valid_simplification(p)
            bool = ischar(p) || iscell(p) && length(p) > 0 && ischar(p{1});
        end
    end
end
