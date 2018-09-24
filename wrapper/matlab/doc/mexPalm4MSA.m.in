%% Description: mexPalm4MSA. 
%  Factorization of a data matrix into multiple factors using PALM.
%  [lambda, facts] = mexPalm4MSA(params) runs the PALM algorithm on the
%  specified set of signals (Algorithm 2 of [1]), returning the factors in
%  "facts" and the multiplicative scalar in "lambda".
%
%  Required fields in PARAMS:
%  --------------------------
%
%    'data' - Training data.
%      A matrix containing the training signals as its columns.
%
%    'nfacts' - Number of factors.
%      Specifies the desired number of factors.
%
%    'cons' - Constraint sets.
%      Specifies the constraint sets in which each factor should lie. It
%      should be a cell-array of size 1*nfacts, where the jth column sets
%      the constraints for the jth factor (starting from the left). cons(j)
%      should be itself a cell-array of size 1*4 taking this form for a
%      factor of size m*n:
%      {'constraint name', 'constraint parameter', m, n}
%
%    'niter' - Number of iterations.
%      Specifies the number of iterations to run.
%
%    'init_facts' - Initialization of "facts".
%      Specifies a starting point for the algorithm.
%
%  Optional fields in PARAMS:
%  --------------------------
%
%    'init_lambda' - Initialization of "lambda".
%      Specifies a starting point for the algorithm. The default value is 1
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
% the project :  <http://faust.inria.fr>
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
