%% Description quick_start
%
%  This demo shows that a Faust is handled as matlab builtin matrix,
% presenting  all the function that are overloaded for Faust class 
% (size,mtimes,transpose...)
% and ends with a little time comparison to illustrate 
% the speed-up of using a Faust for multiplication. 
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
%   Adrien Leman	: adrien.leman@inria.fr
%   Thomas Gautrais	: thomas.gautrais@inria.fr
%	Luc Le Magoarou	: luc.le-magoarou@inria.fr
%	Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%
%%
import matfaust.Faust


% loading a Faust A from saved-one 
A=Faust('faust_quick_start.mat')

% get the size of the faust
[dim1,dim2] = size(A)
 
%  transpose a faust  
A_trans = A'		 
 
% multiplication by A
x1 = rand(dim2,1);
y1 = A*x1;
 
% multiplication by A'
x2 = rand(dim1,5);
y2 = A'*x2;


% get the 2-norm (spectral norm) of the faust A
norm_A = norm(A) % equivalent to norm(A,2);
 
% convert Faust to full matrix
A_full=full(A);

% get the number of non-zeros coefficient
nz = nnz_sum(A)
 
% READING coefficient
coeff=A(3,4)
col_2=A(:,2);
submatrix_A=A(3:5,2:3)
submatrix_A=A(end-5:end,4980:end-1)

% WARNING : WRITING coefficient NOT ALLOWED 
% A(i,j)=3; will throw an error with this message :
%            'function not implemented for Faust class'



 
 
%% speed-up multiplication
nb_mult=100;
time_full=0;
time_faust=0;

for i=1:nb_mult
   tic
      y=A_full*x1;
   time_full=time_full+toc;
   
   tic
      y=A*x1;
   time_faust=time_faust+toc;   
end

disp('multiplication SPEED-UP using Faust');
disp(['Faust is ' num2str(time_full/time_faust) ' faster than a full matrix']);

