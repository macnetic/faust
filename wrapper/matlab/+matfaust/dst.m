%=========================================
%> Returns the Direct Sine Transform (Type II) Faust of order n.
%=========================================
function D = dst(n, varargin)
    import matfaust.Faust
    % bit reversal permutation
    MDFT = mod_fft(n, varargin{:});
    d1 = zeros(1, n);
    for k=1:n
        d1(k) = -1j * exp(-1j * pi / 2 / n * k);
    end
    D1 = -2 * sparse(diag(d1));
    d2 = zeros(1, n);
    for k=1:n
        d2(k) = exp(-1j * pi / n * k);
    end
    D2 = sparse(diag(d2));
    P_ = zeros(n*2, n);
    for i=1:n/2
        P_(i, 2*(i-1)+1) = 1;
        P_(i+n, 2*i) = 1;
    end
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
    O = sparse(diag(vo));
end

function  B = butterfly(N)
    %% Butterfly factor of order N.
    I_N2 = eye(N/2);
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
        Bs{i} = kron(eye(N/N_), B);
        i = i + 1;
        N_ = N_ / 2;
    end
    F = matfaust.Faust([Bs, {complex(bitrev_perm(N))}], varargin{:});
end
