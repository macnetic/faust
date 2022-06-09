function P = bitrev_perm(N)
    index = 1:N;
    new_index = BitReversalPermutation(index);
    P = sparse(index, new_index, ones(1, N), N, N);
end
