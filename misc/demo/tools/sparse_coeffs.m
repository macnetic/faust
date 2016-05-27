function Gamma = sparse_coeffs(D, Ntraining, Sparsity)
%sparse_coeffs. Generation of sparse coefficients
%  Gamma = sparse_coeffs(D, Ntraining, Sparsity) generates Ntraining sparse
%  vectors stacked in a matrix Gamma. Each sparse vector is of size the 
%  number of atoms in the dictionary D, its support is drawn uniformly at
%  random and each non-zero entry is iid Gaussian. 
%
%  References:
%  [1] Le Magoarou L. and Gribonval R., "Learning computationally efficient
%  dictionaries and their implementation as fast transforms", submitted to
%  NIPS 2014

%  Luc Le Magoarou
%  PANAMA team
%  Inria, Rennes - Bretagne Atlantique
%  luc.le-magoarou@inria.fr
%
%  June 2014

Natoms = size(D,2);
Gamma = zeros(Natoms,Ntraining);
for i = 1:Ntraining
    r = randn(Sparsity ,1);
    pos_temp = randperm(Natoms);
    pos = pos_temp(1:Sparsity);
    Gamma(pos,i) = r;
end
end

