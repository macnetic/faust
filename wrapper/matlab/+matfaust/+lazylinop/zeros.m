%=============================================================
%> @brief Returns a zero LazyLinearOp.
%>
%> @param shape: (1x2 real matrix) the shape of the zero LazyLinearOp.
%>
%> @b Example:
%> @code
%> >> import matfaust.lazylinop.zeros
%> >> Lz = zeros([10, 12]);
%> >> x = rand(12, 1);
%> >> Lz * x
%>
%> ans =
%>
%>      0
%>      0
%>      0
%>      0
%>      0
%>      0
%>      0
%>      0
%>      0
%>      0
%> @endcode
%>
%=============================================================
function lz = zeros(shape)
    import matfaust.lazylinop.LazyLinearOp
    import matfaust.lazylinop.LazyLinearOperator
        if ~ ismatrix(shape) || ~ isnumeric(shape) || any(size(shape) ~= [1, 2]) || ~ isreal(shape)
            error('shape is not valid 1x2 double matrix')
        end
    function mm = matmat_(op, shape)
        ops = size(op);
        if ops(1) ~= shape(2)
            error('Dimensions must agree')
        end
        if LazyLinearOp.isLazyLinearOp(op)
            mm = matfaust.lazylinop.zeros([shape(1), size(op, 2)]);
        elseif numel(size(op)) == 2
            mm = zeros([shape(1), size(op, 2)]);
            if isreal(op)
                if strcmp(class(op), 'single')
                    mm = single(mm);
                end
            else
                % op is complex
                mm = complex(mm);
            end
        else % op ndim > 2
            zargs = num2cell(ops(3:end));
            mm = zeros(shape(1), ops(2), zargs{:});
        end
    end
    lz = LazyLinearOperator(shape, 'matmat', @(x) matmat_(x, shape), 'rmatmat', ...
        @(x) matmat_(x, [shape(2), shape(1)]));
end
