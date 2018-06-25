%% Description test_matlab_faust
% This functions test different method
% overloaded for a faust (constructor, size, mtimes, mtimes_trans ...)
% with a given Faust
%
% input parameter : dim1 number of row of the Faust F
%		    dim2 number of column of the Faust F
%		    dim3 number of column of the matrix that will be multiplied by the Faust
%                   nb_fact number of factor of the Faust
%
%
%  
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Nicolas Bellot, Adrien Leman, Thomas Gautrais, 
%			Luc Le Magoarou, Remi Gribonval
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
%   Luc Le Magoarou	: luc.le-magoarou@inria.fr
%   Remi Gribonval	: remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%%


function test_matlab_faust_mult(F,F_dense,x,x_F_trans)
%   test the result of the multiplication between F Faust and matrix x (F*x and F.'*x_F_trans)
%   Detailed explanation goes here

dim1F = size(F,1);
dim2F = size(F,2);

dim1F_dense = size(F_dense,1);
dim2F_dense = size(F_dense,2);

dim1x=size(x,1);
dim1x_F_trans=size(x_F_trans,1);

%% CHECK THE INPUT
if (dim1F ~= dim1F_dense) || (dim2F ~= dim2F_dense)
   error('F and F_dense : size mismatch'); 
end

if (dim1x ~= dim2F)
    error('x (input matrix ) dimension mismatch');
end

if (dim1x_F_trans ~= dim1F)
    error('x_F_trans (input matrix ) dimension mismatch');
end


%% TEST
% expected value for the different multiplication
Y_expected_trans = F_dense.'*x_F_trans;
Y_expected = F_dense*x;
Y_expected_ctrans = F_dense'*x_F_trans;


istransposed=1;
nontransposed=0;



Y_star = F*x;

if (~isequal(Y_expected,Y_star))
    error(['multiplication faust : invalid result  ' ]);
end



Y_star_trans = F.'*x_F_trans;
%(Y_expected_trans~= Y_star_trans)
if (~isequal(Y_expected_trans,Y_star_trans))
    error(['multiplication faust with transposition : invalid result  ' ]);
end


if(isreal(F))
	Y_star_ctrans = F'*x_F_trans;	
	if (~isequal(Y_expected_ctrans,Y_star_ctrans))
   		 error(['multiplication faust with conjugate-transposition : invalid result  ' ]);
	end
end

Y_mtimes_trans = mtimes_trans(F,x_F_trans,istransposed);
%(Y_expected_trans ~= Y_mtimes_trans)
if (~isequal(Y_expected_trans,Y_mtimes_trans))
    error(['multiplication faust with transposition : invalid result  '  ]);
end


Y_mtimes = mtimes_trans(F,x,nontransposed);
%(Y_expected ~= Y_mtimes)
if (~isequal(Y_expected,Y_mtimes))
    error(['multiplication faust : invalid result  '  ]);
end

Y_mtimes_trans_N = mtimes_trans(F.',x_F_trans,nontransposed);
%(Y_expected_trans ~= Y_mtimes_trans_N)
if (~isequal(Y_expected_trans,Y_mtimes_trans_N))
    error(['multiplication faust with transposition : invalid result  '  ]);
end


Y_mtimes_trans_T = mtimes_trans(F.',x,istransposed);
%(y_expected ~= y_mtimes_trans_T)
if (~isequal(Y_expected,Y_mtimes_trans_T))
    error(['multiplication faust : invalid result  '  ]);
end









end
