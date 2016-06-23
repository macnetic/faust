%% Description mexHierarchical_fact 
%  Hierarchical matrix factorization.
%  [lambda, facts, errors] = mexHierarchical_fact(params) runs the hierarchical
%  matrix factorization algorithm (Algorithm 3 of [1])on the specified
%  input matrix, returning the factors in "facts" (cell array of sparse matrices), 
%  the multiplicative scalar in "lambda" and the errors in "errors".
%
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'data' - Training data.
%      A matrix to hierarchically factorize.
%
%    'nfacts' - Number of factors.
%      Specifies the desired number of factors.
%
%    'cons' - Constraint sets.
%      Specifies the constraint sets in which eafaust_params_palm.ch factor should lie. It
%      should be a cell-array of size 2*(nfacts-1), where the jth columns
%      sets the constraints for the jth factorization in two factors:
%      cons(1,j) specifies the constraints for the left factor and
%      cons(2,j) for the right factor. cons(i,j) should be itself a
%      cell-array of size 1*4 taking this form for a factor of size m*n:
%      {'constraint name', 'constraint parameter', m, n}
%
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'niter1' - Number of iterations for the 2-factorisation.
%      Specifies the desired number of iteration for the factorisations in
%      2 factors. The default value is 500.
%
%    'niter2' - Number of iterations for the global optimisation.
%      Specifies the desired number of iteration for the global
%      optimisation. The default value is 500.
%
%    'fact_side' - Side to be factorized iteratively: 1-the left or
%      0-the right. The default value is 0;
%
%    'verbose' - Verbosity of the function. if verbose=1, the function
%      outputs the error at each iteration. if verbose=0, the function runs
%      in silent mode. The default value is 0.
%
%    'update_way' - Way in which the factors are updated. If update_way = 1
%      ,the factors are updated from right to left, and if update_way = 0,
%      the factors are updated from left to right. The default value is 0.
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Luc Le Magoarou, Remi Gribonval
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
%   Nicolas Bellot : nicolas.bellot@inria.fr
%   Adrien Leman   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%
