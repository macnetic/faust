%=========================================
%> @brief Returns the Direct Sine Transform (Type II) Faust of order n.
%>
%> The analytical formula of DST II used here is:
%> \f$2 \sum_{n=0}^{N-1} x_n sin \left( {\pi (k+1) (2n + 1)} \over {2N} \right)\f$
%>
%> @param n: the order of the DST (it must be a power of two).
%> @param 'dev', str: 'gpu' or 'cpu' to create the Faust on CPU or GPU ('cpu' is the default).
%>
%> @b Example
%> @code
%> % in a matlab terminal
%> >> import matfaust.*
%> >> F = dct(8);
%> >> x = ones(8, 1);
%> >> % apply the DCT to x
%> >> real(F*x)
%>
%> ans =
%>
%>   10.2517
%>    0.0000
%>    3.5999
%>    0.0000
%>    2.4054
%>    0.0000
%>    2.0392
%>    0.0000
%> @endcode
%>
%>@b See also matfaust.dft, matfaust.dct
%=========================================
function D = dst(n, varargin)
    import matfaust.Faust
    % bit reversal permutation
    MDFT = mod_fft(n, varargin{:});
    d1 = zeros(1, n);
    for k=1:n
        d1(k) = -1j * exp(-1j * pi / 2 / n * k);
    end
    % D1 = -2 * sparse(diag(d1));
    D1 = -2 * sparse(1:n, 1:n, d1, n, n, n);
    d2 = zeros(1, n);
    for k=1:n
        d2(k) = exp(-1j * pi / n * k);
    end
    % D2 = sparse(diag(d2));
    D2 = sparse(1:n, 1:n, d2, n, n, n);
    %    P_ = zeros(n*2, n);
    %    for i=1:n/2
    %        P_(i, 2*(i-1)+1) = 1;
    %        P_(i+n, 2*i) = 1;
    %    end
    % P_ as sparse matrix
    P_row_inds = [ 1:n/2, 1+n:n/2+n ];
    P_col_inds = [ 1:2:n-1, 2:2:n ];
    P_ = sparse(P_row_inds, P_col_inds, ones(n, 1), n*2, n, n);
    F_even = Faust(D1, varargin{:}) * MDFT;
    F_odd = Faust(D1, varargin{:}) * Faust(D2, varargin{:}) * MDFT;
    F = [F_even, F_odd];
    F = F * Faust(P_);
    D = real(F);
end

function P = bitrev_perm(N)
    index = 1:N;
    new_index = BitReversalPermutation(index);
    P = sparse(index, new_index, ones(1, N), N, N);
end

function O = omega(N)
    %% List of n-th root of unit raised to the power of -(k+1) (instead of k as
    %in FFT, the purpose is to write the DST).
    o = exp(pi*1j/N);
    vo = zeros(1, N);
    for k=1:N
        vo(k) = o^-k;
    end
    % O = sparse(diag(vo));
    O = sparse(1:N, 1:N, vo, N, N, N);
end

function  B = butterfly(N)
    %% Butterfly factor of order N.
    I_N2 = speye(N/2);
    O_N2 = omega(N/2);
    B = [[ I_N2, O_N2]; [I_N2, - O_N2]];
end

function F = mod_fft(N, varargin)
    %% Modified FFT for DST computation (it would have been the standard FFT
    %% if omega was raised to the power of -k instead of -(k+1)).
    N_ = N;
    Bs = cell(1, log2(N));
    i = 1;
    while N_ ~= 1
        B = butterfly(N_);
        Bs{i} = kron(speye(N/N_), B);
        i = i + 1;
        N_ = N_ / 2;
    end
    F = matfaust.Faust([Bs, {complex(bitrev_perm(N))}], varargin{:});
end
