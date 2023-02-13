%=============================================================
%> @brief Returns a LazyLinearOp of op padded with zeros on one or two dimensions.
%>
%> @param op: the array/operator to pad with zeros.
%> @param z_sizes: a cell array of vectors of two integers. The first vector defined the number of zeros to prepend and append to op on row dimension. The second vector is the same for the column dimension.
%>
%> @b Example:
%> @code
%> >> import matfaust.lazylinop.zpad
%> >> A = reshape(1:18*2, 18, 2)
%>
%> A =
%>
%>      1    19
%>      2    20
%>      3    21
%>      4    22
%>      5    23
%>      6    24
%>      7    25
%>      8    26
%>      9    27
%>     10    28
%>     11    29
%>     12    30
%>     13    31
%>     14    32
%>     15    33
%>     16    34
%>     17    35
%>     18    36
%>
%> >> lz = zpad(A, { [2, 3]; [4, 1]});
%> >> lz
%>
%> lz =
%>
%>   23x7 LazyLinearOp array with no properties.
%>
%> >> full(lz)
%>
%> ans =
%>
%>      0     0     0     0     0     0     0
%>      0     0     0     0     0     0     0
%>      0     0     0     0     1    19     0
%>      0     0     0     0     2    20     0
%>      0     0     0     0     3    21     0
%>      0     0     0     0     4    22     0
%>      0     0     0     0     5    23     0
%>      0     0     0     0     6    24     0
%>      0     0     0     0     7    25     0
%>      0     0     0     0     8    26     0
%>      0     0     0     0     9    27     0
%>      0     0     0     0    10    28     0
%>      0     0     0     0    11    29     0
%>      0     0     0     0    12    30     0
%>      0     0     0     0    13    31     0
%>      0     0     0     0    14    32     0
%>      0     0     0     0    15    33     0
%>      0     0     0     0    16    34     0
%>      0     0     0     0    17    35     0
%>      0     0     0     0    18    36     0
%>      0     0     0     0     0     0     0
%>      0     0     0     0     0     0     0
%>      0     0     0     0     0     0     0
%>
%> @endcode
%>
%=============================================================
function out = zpad(op, z_sizes)

    import matfaust.lazylinop.zeros
    import matfaust.lazylinop.aslazylinearoperator
    import matfaust.lazylinop.LazyLinearOp
    % TODO: sanitize op
    if ~ iscell(z_sizes)
        error('z_sizes must be a cell array')
    end
    z_sizes = cell2mat(z_sizes);
    if numel(size(z_sizes)) > 2
        error('Cannot pad zeros on more than two dimensions')
    end
    if ~ LazyLinearOp.isLazyLinearOp(op)
        op = aslazylinearoperator(op);
    end
    out = op;
    for i=1:size(z_sizes, 1)
        bz = z_sizes(i, 1);
        az = z_sizes(i, 2);
        if bz > 0
            if i == 1
                out = vertcat(matfaust.lazylinop.zeros([bz, size(out, 2)]), out);
            else
                out = horzcat(matfaust.lazylinop.zeros([size(out, 1), bz]), out);
            end
        if az > 0
            if i == 1
                out = vertcat(out, matfaust.lazylinop.zeros([az, size(out, 2)]));
            else
                out = horzcat(out, matfaust.lazylinop.zeros([size(out, 1), az]));
            end
        end
    end
end
