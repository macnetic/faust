%% Description demo_fact_hadamard
%
%  This demo hierarchically factorizes the Hadamard dictionary and then
%  plots the results. This essentially reproduces figure 2 from [1].
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
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




n = 32;
M = log2(n);

%% Generating the data
[matrix,~] = hadamard_mat(M);





%% Setting of the parameters
params.nrow=n;
params.ncol=n;
params.nfacts = M;
params.cons = cell(2,M-1);

for j=1:M-1
    params.cons{1,j} = {'splincol',2,n,n};
    params.cons{2,j} = {'splincol',n/2^j,n,n};
end

params.niter1 = 30;
params.niter2 = 30;
params.update_way = 1;
params.verbose = 0;


hadamard_faust = faust_decompose(matrix,params);
Xhat = full(hadamard_faust);
relative_error = norm(matrix - Xhat)/norm(matrix);
fprintf(['\n\n relative error between hadamard matrix and its transform : ' num2str(relative_error) '\n']);

nb_mult= 500;
dense_times=zeros(nb_mult,1);
faust_times=zeros(nb_mult,1);

for i=1:nb_mult
   x=rand(n,1);
   y_faust=rand(n,1);
   y_X=rand(n,1);
   
   tic 
   y_X=matrix*x;
   t1=toc;
   
   tic 
   y_faust=hadamard_faust*x;
   t2=toc;
   
    dense_times(i)=t1;
    faust_times(i)=t2;
    
    
end


dense_t=mean(dense_times);
faust_t=mean(faust_times);
speed_up=dense_t/faust_t;

disp(['multiplication by hadamard matrix : ' num2str(dense_t) 'sec']);
disp(['multiplication by faust : ' num2str(faust_t) 'sec']);
disp(['multiplication speed-up using faust : ' num2str(speed_up) ]);


%% Plotting the results

% figure configuration
runPath=which(mfilename);
pathname = fileparts(runPath);
figure_dir = [pathname filesep '..' filesep 'Figures'];
format_fig='-dpng';


fighandle=figure;
hold on;
subplot(1,params.nfacts+1,1);
imagesc(Xhat); axis square
set(gca,'xtick',[],'ytick',[])

for kk = 1:params.nfacts
    subplot(1,params.nfacts+1,1+kk)
    imagesc(facts{kk}); axis square
    set(gca,'xtick',[],'ytick',[])
end

fighandle.Name =['Hadamard-factorisation_image'];
figure_name = [figure_dir filesep 'Hadamard-factorisation_image'];
print(figure_name, format_fig);


fighandle=figure;
hold on;
subplot(1,params.nfacts+1,1);
imagesc(Xhat); axis square
set(gca,'xtick',[],'ytick',[])


for kk = 1:params.nfacts
    subplot(1,params.nfacts+1,1+kk)
    spy(facts{kk}); axis square
    set(gca,'xtick',[],'ytick',[])
end
fighandle.Name =['Hadamard-factorisation_nnz_coeff'];
figure_name = [figure_dir filesep 'Hadamard-factorisation_nnz_coeff'];
print(figure_name, format_fig);




