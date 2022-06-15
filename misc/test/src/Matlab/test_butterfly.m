function test_butterfly()
    import matfaust.rand_butterfly
    import matfaust.wht
    import matfaust.fact.butterfly
    N = 32;
    P_funcs = {@(N) bitrev_perm(N), @(N) rand_perm(N)};
    P_types = {'bitrev_perm', 'rand_perm'};
    M_funcs = {@(N) rand_butterfly(N), @(N) rand_butterfly(N, 'field', 'complex'), @(N) wht(N)};
    M_types = {'rand_butterfly_real', 'rand_butterfly_complex', 'hadamard'};
    em = {'expectation: should be small', 'expectation: should be large'};
    expectations = {em{1}, em{2}, em{2}, em{2}, em{1}, em{2}, em{2}, em{2}, em{1}};
    for ind_P=1:length(P_funcs)
        P_func = P_funcs{ind_P};
        P = P_func(N);
        P_type = P_types{ind_P};
        for ind_M=1:length(M_funcs)
            M_func = M_funcs{ind_M};
            M_type = M_types{ind_M};
            M = M_func(N);
            disp(['======================= matrix: ', M_type, ' P: ', P_type])
            assert(all(all(full(indices2perm(perm2indices(P))) == full(P))))
            Ma = full(M);
            if ~ isreal(Ma)
                P = complex(P);
            else
                P = real(P);
            end
            G1 = M*P;
            G2 = M*P.';
            F = cell(1, 9);
            for i=0:(length(F)-1)
                if i < 3
                    disp('mat=M')
                    mat = Ma;
                elseif i < 6
                    disp('mat=M*P')
                    mat = G1;
                else
                    disp('mat=M*P.''')
                    mat = G2;
                end
                if mod(i, 3)  == 0
                    disp('perm=None')
                    p = [];
                elseif mod(i, 3) == 1
                    disp('perm=P')
                    p = P;
                else
                    disp('perm=P.''')
                    p = P.';
                end
                butterfly_args = {mat, 'type', 'bbtree'};
                if numel(p) > 0
                    butterfly_args = [ butterfly_args{:}, 'perm', {perm2indices(p)}];
                end
                F{i+1} = butterfly(butterfly_args{1:end});
                disp(['F', int2str(i+1), ' err: ', num2str(err(mat, F{i+1})), ' ', expectations{i+1}])
                if ~ is_as_expected(err(mat, F{i+1}), expectations{i+1})
                    disp('(NOT OK)')
                    save([M_type, 'mat-', P_type, '-', int2str(i+1), '.mat'], 'mat', 'p', 'P')
                    save(F{i+1}, ['F', int2str(i+1), '.mat'])
                end
                disp('-----')
            end
        end
    end
end

function e = err(a, F)
    e = norm(F-a) / norm(a);
end

function b = is_as_expected(err, expected)
    b =  err < 1e-12 && numel(strfind(expected, 'small')) || ...
            err > 1e-12 && numel(strfind(expected, 'large'));
end

function P = rand_perm(n)
    J = randperm(n);
    P = indices2perm(J);
end

function J = perm2indices(P)
    [J, ~, ~] = find(P);
    J = J.';
end

function P = indices2perm(J)
    n = numel(J);
    P = zeros(n);
    for j=1:n
        P(J(j), j) = 1;
    end
    P = sparse(P);
end
