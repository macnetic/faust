function [Fourier,facts] = CooleyTukeyFact(dim)
%  Cooley-Tukey Factorization of the Fourier Matrix.
%
%  [Fourier,facts] = CooleyTukeyFact(dim) factorize the Fourier matrix of dimension dim
%                    Fourier is the full-storage Fourier matrix
%                    facts is a cell-array of sparse matrix
%                    
%		     Fourier = facts{1}*facts{2}*...*facts{log2(dim)+1}.
%                    
%                    warning : dim must be a power of 2
%                               

%
%% source http://www.cs.cornell.edu/~bindel/class/cs5220-s10/slides/FFT.pdf
% 
%
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot,  Luc Le Magoarou, Remi Gribonval
%			INRIA Rennes, FRANCE
%			http://www.inria.fr/
%
% The FAuST Toolbox is distributed under the terms of the GNU Affero
% General Public License.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% Contacts:
%   Nicolas Bellot	: nicolas.bellot@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
%	approximations of matrices and applications", Journal of Selected
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>

nb_fact = log2(dim)+1;

if (floor(nb_fact) ~= nb_fact)
	error('dim must be a power of 2');
end


facts=cell(1,nb_fact);
%% init permutation matrix
index=1:dim;
new_index = BitReversalPermutation(index);
P=sparse(index,new_index,ones(1,dim),dim,dim);
facts{nb_fact}=P;



for q=1:nb_fact-1
    L=2^q;
    r=dim/L;
    disp([ 'L*r' : num2str(L*r)]);
    wL=exp(-2*pi*i/L);
    index=1:L/2;
    OmegaDiag=wL.^(index-1);
    Omega=sparse(index,index,OmegaDiag,L/2,L/2);
    %%Omega=diag(OmegaDiag);
    Rhs=[speye(L/2) Omega;speye(L/2) -Omega];
    Aq = kron(speye(r),Rhs);
    facts{nb_fact-q}=Aq;
end


Fourier = dvp(facts);



end

