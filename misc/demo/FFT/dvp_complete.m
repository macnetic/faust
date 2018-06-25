

function [prod,residus] = dvp_complete(facts,method)
%  Development of the input factorized matrix.
%  [Prod,Residus] = dvp_complete(Facts,Method) develops the cell-array of matrices Facts into the matrix Prod
%  which is the product of the matrices contained in Facts:
%  Prod = Facts{1}*Facts{2}*...*Facts{n}.
% 
% Residus is a cell-array of matrices correspond to the cumulated product of the Facts  
% if Method = 'L2R' i.e left to right
% Residus{i} = Facts{1} *  Facts{2} * ... * Facts{i}
% if Method = 'R2L' i.e right to left
% Residus{i} = Facts{n} *  Facts{n-1} * ... * Facts{n+1-i}
% if Method = 'Central' i.e right to left, middle_id is the index of the central factor 
% Residus{i} = Facts{middle_id - (i-1)} * Facts{middle_id} *  Facts{middle_id + (i-1)}
%
%
% For more information on the FAuST Project, please visit the website of
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, Luc Le Magoarou, Remi Gribonval
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




nb_factor=length(facts);

nRow=size(facts{1},1);
nCol=size(facts{end},2);

if (strcmp(method,'Central'))
    
    middle = ceil(nb_factor/2);
    prod = facts{middle};
    residus=cell(1,middle);
    residus{1} = prod;
    
    for i=1:middle-1
       prod = facts{middle-i}*prod*facts{middle+i};
       residus{i+1} = prod;
    end
else
    
    switch method
        case 'L2R'
            dim=nRow;
            id_fact=1:nb_factor;
            mult = @(x,y) x*y;
        case 'R2L'
            dim=nCol;
            id_fact=nb_factor:-1:1;
            mult = @(x,y) y*x;
        otherwise
            error('invalid method');
    end
    
    prod=speye(dim);
    residus=cell(nb_factor,1);
    
    for i=1:nb_factor
        prod= mult(prod,facts{id_fact(i)});
        residus{i}=prod;
    end
end


end

