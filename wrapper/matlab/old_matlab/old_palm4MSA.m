function [lambda, facts] = old_palm4MSA(params)
%% Description: palm4MSA. 
%  Factorization of a data matrix into multiple factors using PALM.
%  [lambda, facts] = palm4MSA(params) runs the PALM algorithm on the
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


%%%% Setting optional parameters values %%%%
if (isfield(params,'init_lambda'))
    init_lambda = params.init_lambda ;
else
    init_lambda = 1;
end

if (isfield(params,'verbose'))
    verbose = params.verbose ;
else
    verbose = 0;
end

if (isfield(params,'update_way'))
    update_way = params.update_way;
else
    update_way = 0;
end

if params.nfacts ~= size(params.init_facts,2)
    error('Wrong initialization: params.nfacts and params.init_facts are in conflict')
end

%%%%%% Setting of the constraints %%%%%%%
handles_cell = cell(1,params.nfacts);
for ii=1:params.nfacts
    cons = params.cons{ii};
    if strcmp(cons{1},'sp')
        handles_cell{ii} = @(x) prox_sp(x,cons{2});
    elseif strcmp(cons{1},'spcol')
        handles_cell{ii} = @(x) prox_spcol(x,cons{2});
    elseif strcmp(cons{1},'splin')
        handles_cell{ii} = @(x) prox_splin(x,cons{2});
    elseif strcmp(cons{1},'normcol')
        handles_cell{ii} = @(x) prox_normcol(x,cons{2});
    elseif strcmp(cons{1},'splincol')
        handles_cell{ii} = @(x) prox_splincol(x,cons{2});
    elseif strcmp(cons{1},'l0pen')
        handles_cell{ii} = @(x) prox_l0pen(x,cons{2});
    elseif strcmp(cons{1},'l1pen')
        handles_cell{ii} = @(x) prox_l1pen(x,cons{2});
    elseif strcmp(cons{1},'const')
        handles_cell{ii} = @(x) cons{2};
    elseif strcmp(cons{1},'wav')
        handles_cell{ii} = @(x) prox_wav(x,cons{2});
    elseif strcmp(cons{1},'sppos')
        handles_cell{ii} = @(x) prox_sp_pos(x,cons{2});
    elseif strcmp(cons{1},'blkdiag')
        handles_cell{ii} = @(x) prox_blkdiag(x,cons{2});
    elseif strcmp(cons{1},'splin_test')
        handles_cell{ii} = @(x) prox_splin_test(x,cons{2});
    elseif strcmp(cons{1},'supp')
        handles_cell{ii} = @(x) prox_supp(x,cons{2});
    elseif strcmp(cons{1},'normlin')
        handles_cell{ii} = @(x) prox_normlin(x,cons{2});
    else
        error('The expressed type of constraint is not known')
    end
end

% Initialization
facts = params.init_facts;
lambda = init_lambda;

%%%%%% Main Loop %%%%%%%
X = params.data;
if update_way
    maj = params.nfacts:-1:1;
else
    maj = 1:params.nfacts;
end

for i = 1:params.niter
    for j = maj
        if strcmp(params.cons{j}{1},'const')
            facts{j} = handles_cell{j}(facts{j});
        else
            L = facts(1:j-1);
            R = facts(j+1:end);
            if isreal(X)
                [grad, LC] = grad_comp(L,facts{j},R,params.data,lambda);
            else
                [grad, LC] = grad_comp_cpx(L,facts{j},R,params.data,lambda);
            end
            c = LC*1.001;%1e4;%LC*1.001;%*
            cons = params.cons{j};
            if strcmp(cons{1},'l0pen')
                facts{j} = prox_l0pen(facts{j} - (1/c)*grad,sqrt(2*cons{2}/c));
            elseif strcmp(cons{1},'l1pen')
                facts{j} = prox_l1pen(facts{j} - (1/c)*grad,cons{2}/c);
            else
                facts{j} = handles_cell{j}(facts{j} - (1/c)*grad);
            end
        end
    end
    lambda = trace(mult_right(X',facts))/trace(mult_right(dvp(facts)',facts));

    %lambda = trace(X'*dvp(facts))/trace(dvp(facts)'*dvp(facts));
    if verbose
        disp(['Iter ' num2str(i) ', RMSE=' num2str(norm(X-lambda*dvp(facts),'fro')/sqrt(numel(X))) ])
    end
end
end
