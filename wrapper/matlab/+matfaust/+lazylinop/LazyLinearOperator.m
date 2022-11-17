%=============================================================
%> Returns a LazyLinearOp L defined by shape and at least a matvec or matmat function.
%>
%> @param shape: dimensions of the operator (M, N),
%> @param 'matvec', function_handle: returns L * v (v a vector of size N).
%> @param 'rmatvec', function_handle: (optional) returns L' * v (v a vector of size M).
%> @param 'matmat', function_handle: returns A * V (V a dense matrix of dimensions (N, K)).
%> @param 'rmatmat', function_handle: (optional) returns A' * V (V a dense matrix of dimensions (M, K)).
%> @param dtype, str: (optional) complex, double, single or undefined (by default).
%>
%>
%> @b Example:
%> In this example we create a LazyLinearOp for the DFT using the fft matlab
%> built-in.
%> @code
%> >> import matfaust.lazylinop.LazyLinearOperator
%> >> lfft = LazyLinearOperator([8, 8], 'matvec', @(x) fft(x), 'rmatvec', @(x) 8 * ifft(x))
%> >> x = rand(8, 1);
%> >> norm(lfft * x - fft(x)) < 1e-12
%>
%> ans =
%>
%>  logical
%>
%>  1
%> @endcode
%>
%=============================================================
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
