%% Description: setup_FAuST.m 
% Run script to set the useful paths in order to use the matlab wrapper of the C++ FAuST toolbox.
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.inria.fr>
%
%% License:
% Copyright (2016):	Bellot Nicolas, Adrien Leman, Luc Le Magoarou, Remi Gribonval
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


ROOT_DIR=[fileparts(mfilename('fullpath')),filesep];


             
fprintf('Welcome to the Matlab wrapper of the FAuST C++ toolbox.');
fprintf('FAuST root directory is %s\n',ROOT_DIR);
fprintf('Adding path %s and subdirectories to matlab path\n',ROOT_DIR);
%addpath(genpath(ROOT_DIR));
addpath(ROOT_DIR);%% avoiding to use genpath because if +matfaust is added to path eye.m, rand.m, version.m etc. will conflict with matlab eye, etc.
addpath([ROOT_DIR '/mex']);
if(startsWith(computer, 'PCWIN'))
	% executing the script on Windows
	addpath([ROOT_DIR '/mex/Release'])
end
addpath([ROOT_DIR '/data']);
addpath([ROOT_DIR '/tools']);
matfaust.enable_gpu_mod()







fprintf(['\n\n To get started with the FAuST Toolbox : \n'])
fprintf(['\n launch quick_start or run_all_demo.m \n'])


