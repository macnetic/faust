%=========================================
%> Returns  a toeplitz Faust whose the first column is c and the first row r.
%>
%> @b Usage
%>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c), T is a symmetric Toeplitz Faust whose the first column is c. <br/>
%> &nbsp;&nbsp;&nbsp; @b T = toeplitz(c, r), T is a Toeplitz Faust whose the first column is c and the first row is [c(1), r(2:)] <br/>

%> @param c: the first column of the toeplitz Faust.
%> @param r: (2nd argument) the first row of the toeplitz Faust. Defaulty r =
%>         c. r(1) is ignored, the first row is always [c(1),
%>         r(2:)].
%>
%> @b See also matfaust.circ, matfaust.anticirc
%=========================================
function T = toeplitz(c, varargin)
    if (length(varargin) > 0)
        r = varargin{1};
        if(~ ismatrix(r) || size(r, 1) ~= 1 && size(r, 2) ~= 1)
            error('The second argument must be a vector')
        end
    else
        r = c; % default r
    end
    if size(c, 2) == 1
        c = c.';
    elseif size(c, 1) ~= 1
        error('c must be a vector')
    end
    if size(r, 2) == 1
        r = r.';
    elseif size(r, 1) ~= 1
        error('r must be a vector')
    end
    m = numel(c);
    n = numel(r);
    N = 2 ^ ceil(log2(max(m, n)));
    c_ = [c, zeros(1, N-m+1+N-n), r(end:-1:2)];
    C = matfaust.circ(c_);
    T = C(1:m, 1:n);
end
