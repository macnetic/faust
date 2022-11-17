function L = LazyLinearOperator(shape, varargin)
    import matfaust.lazylinop.LazyLinearOp
    p = inputParser;
    validDtype = @(dtype) any(strcmp(dtype, {'complex', 'double', 'single', 'undefined'}));
    addOptional(p, 'dtype', 'undefined', validDtype)
    addOptional(p, 'matvec', 'none', @(l) isa(l, 'function_handle'))
    addOptional(p, 'rmatvec', 'none', @(l) isa(l, 'function_handle'))
    addOptional(p, 'matmat', 'none', @(l) isa(l, 'function_handle'))
    addOptional(p, 'rmatmat', 'none', @(l) isa(l, 'function_handle'))

    parse(p, varargin{:})
    dtype = p.Results.dtype;
    matvec = p.Results.matvec;
    rmatvec = p.Results.rmatvec;
    matmat = p.Results.matmat;
    rmatmat = p.Results.rmatmat;

    if strcmp(matvec, 'none') && strcmp(matmat, 'none')
        error(['At least a matvec or a matmat function must be passed in'
        ' kwargs'])
    end

    % define a matmat function from matvec in case matmat is not provided
    function out = matmat_(M, matvec_)
        out = zeros(shape(1), size(M, 2));
        for i=1:size(M, 2)
            out(:, i) = matvec_(M(:, i));
        end
    end

    if strcmp(matmat, 'none')
        matmat = @(M) matmat_(M, matvec);
    end

    if strcmp(rmatmat, 'none') && ~ strcmp(rmatvec, 'none')
        rmatmat = @(M) matmat_(M, rmatvec);
    end

    L = LazyLinearOp.create_from_funcs(matmat, rmatmat, shape);
end
