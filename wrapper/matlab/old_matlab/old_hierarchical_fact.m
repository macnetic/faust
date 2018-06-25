function [lambda, facts, errors] = old_hierarchical_fact(matrix,params)
%% Description hierarchical_fact 
%  Hierarchical matrix factorization.
%  [lambda, facts, errors] = hierarchical_fact(matrix,params) runs the hierarchical
%  matrix factorization algorithm (Algorithm 3 of [1])on the specified
%  input matrix, returning the factors in "facts" (cell array of sparse matrices), 
%  the multiplicative scalar in "lambda" and the errors in "errors".
%
%
%  
%
%
%  Required fields in PARAMS:
%  --------------------------
%
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
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Luc Le Magoarou, Remi Gribonval
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
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%



%%%% Setting parameters values %%%%
if (isfield(params,'niter1'))
    niter1 = params.niter1;
else
    niter1 = 500;
end

if (isfield(params,'niter2'))
    niter2 = params.niter2;
else
    niter2 = 500;
end

if (isfield(params,'verbose'))
    verbose = params.verbose;
else
    verbose = 0;
end

if (isfield(params,'update_way'))
    update_way = params.update_way;
else
    update_way = 0;
end

if (isfield(params,'fact_side'))
    fact_side = params.fact_side ;
else
    fact_side = 0;
end



%%%% Verifying the validity of the constraints %%%%
verifSize = params.nrow == params.cons{1,1}{3} && params.cons{1,1}{4}...
    == params.cons{2,1}{3} && params.ncol == params.cons{2,1}{4};
for i = 2:params.nfacts-1
    if fact_side
        verifSize = verifSize && params.cons{2,i-1}{3} == params.cons{2,i}{4}...
            && params.cons{1,i}{4} == params.cons{2,i}{3} && params.nrow == params.cons{1,i}{3};
    else
        verifSize = verifSize && params.cons{1,i-1}{4} == params.cons{1,i}{3}...
            && params.cons{1,i}{4} == params.cons{2,i}{3} && params.ncol == params.cons{2,i}{4};
    end
end

if ~verifSize
    error('Size incompatibility in the constraints')
end

if params.nfacts-1 ~= size(params.cons,2)
    error('The number of constraints is in conflict with the number of factors')
end


if ((size(matrix,1) ~= params.nrow) || (size(matrix,2) ~= params.ncol))
	error('The config doesn''t match the size of matrix to be factorized')
end



% Initialization
lambda = 1;
facts = cell(1,params.nfacts);
Res = matrix;
errors = zeros(params.nfacts-1,2);

%%%%%% Main loop %%%%%%
for k=1:params.nfacts-1
    cons = params.cons(:,k);
    
    %%%% Factorization in 2 %%%%
    params2.niter = niter1;
    params2.nfacts = 2;
    params2.data = Res;
    params2.verbose = verbose;
    params2.update_way = update_way;
    params2.cons = [cons(1), cons(2)];
    params2.init_facts = {zeros(cons{1}{3},cons{1}{4}), eye(cons{2}{3},cons{2}{4})};
    if update_way
        params2.init_facts = {eye(cons{1}{3},cons{1}{4}), zeros(cons{2}{3},cons{2}{4})};
    end
    params2.init_lambda = 1;
    [lambda2, facts2] = old_palm4MSA(params2);
    %[lambda2, facts2] = palm4MSA_backprop(params2);
    if fact_side
        facts(3:end) = facts(2:end-1);
        facts(1:2) = facts2;
    else
        facts(k:k+1) = facts2;
    end
    lambda = lambda*lambda2;
    
    %%%% Global optimization %%%%
    params3.niter = niter2;
    params3.nfacts = k+1;
    params3.data = matrix;
    params3.verbose = verbose;
    params3.update_way = update_way;
    if fact_side
        params3.cons = [cons(1), params.cons(2,k:-1:1)];
    else
        params3.cons = [params.cons(1,1:k), cons(2)];
    end   
    params3.init_facts = facts(1:k+1);
    params3.init_lambda = lambda;
    [lambda, facts3] = old_palm4MSA(params3);
    %[lambda, facts3] = palm4MSA_backprop(params3);
    facts(1:k+1) = facts3;
    if fact_side
        Res = facts3{1};
    else
        Res = facts3{k+1};
    end
    errors(k,1) = norm(matrix - lambda*dvp(facts3))/norm(matrix);
    errors(k,2) = nnzero_count(facts3)/numel(matrix);   
end
end
