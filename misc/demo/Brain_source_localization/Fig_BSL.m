%% Description Fig_BSL 
%  Brain source localization figures
%  This script builds the BSL figure (Fig. 9.) used in [1].
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
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
%   Nicolas Bellot : nicolas.bellot@inria.fr
%   Leman Adrien   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%
% [2]   A. Gramfort, M. Luessi, E. Larson, D. Engemann, D. Strohmeier, 
%	C. Brodbeck, L. Parkkonen, M. Hamalainen, MNE software for processing
%	MEG and EEG data <http://www.ncbi.nlm.nih.gov/pubmed/24161808>, 
%	NeuroImage, Volume 86, 1 February 2014, Pages 446-460, ISSN 1053-8119, 
%	[DOI] <http://dx.doi.org/10.1016/j.neuroimage.2013.10.027>
%%

runPath=which(mfilename);
pathname = fileparts(runPath);
matfile = fullfile(pathname, 'output/results_BSL_user.mat');
if (not(exist(matfile)))
    error('run BSL.m before Fig_BSL.m');
end
load(matfile);
%% convergence analysis

d1 = cat(4,resDist(:,1,1,:),resDist(:,1,2,:));
d2 = cat(4,resDist(:,2,1,:),resDist(:,2,2,:));
d3 = cat(4,resDist(:,3,1,:),resDist(:,3,2,:));
test2 = 100*[squeeze(d1);zeros(1,1000);squeeze(d2);zeros(1,1000);squeeze(d3)];


figure('color',[1 1 1]);
T = bplot(test2','linewidth',1.5);
legend(T)

ylabel('Distance between true and estimated sources (cm)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1;nb_approx_MEG+3:2*nb_approx_MEG+3;2*nb_approx_MEG+5:3*nb_approx_MEG+5];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',reshape(x_tick',1,numel(x_tick)));
set(ax,'xticklabel', [])
axi = axis;

axis([0 3*(nb_approx_MEG+2) -0.3 4.7])
HorizontalOffset = 10;
verticalOffset = 0.43;
verticalOffset2 = 0.7;
yTicks = get(ax,'ytick');
xTicks = get(ax, 'xtick');
minY = min(yTicks);

text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(2+nb_approx_MEG), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(1+2*(nb_approx_MEG+1)), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(1:(nb_approx_MEG+1)))), minY - verticalOffset2, '$1<d<5$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(nb_approx_MEG+2:2*(nb_approx_MEG+1)))), minY - verticalOffset2, '$5<d<8$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(2*(nb_approx_MEG+1)+1:3*(nb_approx_MEG+1)))), minY - verticalOffset2, '$d>8$','HorizontalAlignment','center','interpreter', 'latex');





for i=1:nb_approx_MEG
    text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');
    text(xTicks(i+2+nb_approx_MEG), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');
    text(xTicks(i+1+2*(nb_approx_MEG+1)), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');
    
end

title('Fig 9 : C++ wrapper faust');


%% MATLAB convergence analysis

d1 = cat(4,resDist_matlab(:,1,1,:),resDist_matlab(:,1,2,:));
d2 = cat(4,resDist_matlab(:,2,1,:),resDist_matlab(:,2,2,:));
d3 = cat(4,resDist_matlab(:,3,1,:),resDist_matlab(:,3,2,:));
test2 = 100*[squeeze(d1);zeros(1,1000);squeeze(d2);zeros(1,1000);squeeze(d3)];

figure('color',[1 1 1]);
title('MATLAB');
T = bplot(test2','linewidth',1.5);
legend(T)

ylabel('Distance between true and estimated sources (cm)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1;nb_approx_MEG+3:2*nb_approx_MEG+3;2*nb_approx_MEG+5:3*nb_approx_MEG+5];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',reshape(x_tick',1,numel(x_tick)));
set(ax,'xticklabel', [])
axi = axis;

axis([0 3*(nb_approx_MEG+2) -0.3 4.7])
HorizontalOffset = 10;
verticalOffset = 0.43;
verticalOffset2 = 0.7;
yTicks = get(ax,'ytick');
xTicks = get(ax, 'xtick');
minY = min(yTicks);

text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(2+nb_approx_MEG), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(1+2*(nb_approx_MEG+1)), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(1:(nb_approx_MEG+1)))), minY - verticalOffset2, '$1<d<5$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(nb_approx_MEG+2:2*(nb_approx_MEG+1)))), minY - verticalOffset2, '$5<d<8$','HorizontalAlignment','center','interpreter', 'latex');
text(xTicks(round(mean(2*(nb_approx_MEG+1)+1:3*(nb_approx_MEG+1)))), minY - verticalOffset2, '$d>8$','HorizontalAlignment','center','interpreter', 'latex');

for i=1:nb_approx_MEG
    text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');
    text(xTicks(i+2+nb_approx_MEG), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');
    text(xTicks(i+1+2*(nb_approx_MEG+1)), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) ' }$' ],'HorizontalAlignment','center','interpreter', 'latex');
    
end

title('Fig 9 : matlab faust');



%% time comparison

timeS = cat(3,compute_Times(:,1,:),compute_Times(:,2,:),compute_Times(:,3,:));


timeS = squeeze(1000*timeS)';
figure('color',[1 1 1]);
T = bplot(timeS,'linewidth',1.5);
legend(T)
ylabel('Computed Time (ms)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',x_tick');
set(ax,'xticklabel', [])
axi = axis;


HorizontalOffset = 10;
verticalOffset = 0.5;
verticalOffset2 = 1;
yTicks = get(ax,'ytick');
xTicks = get(ax, 'xtick');
minY = min(yTicks);

text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');


for i=1:nb_approx_MEG
    text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');
    
end
title(' C++ wrapper faust');



%% MATLAB time comparison

timeS_matlab = cat(3,compute_Times_matlab(:,1,:),compute_Times_matlab(:,2,:),compute_Times_matlab(:,3,:));


timeS_matlab = squeeze(1000*timeS_matlab)';

figure('color',[1 1 1]);
T = bplot(timeS_matlab,'linewidth',1.5);
legend(T)

ylabel('Computed Time (ms)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',x_tick');
set(ax,'xticklabel', [])
axi = axis;


HorizontalOffset = 10;
verticalOffset = 0.5;
verticalOffset2 = 1;
yTicks = get(ax,'ytick');
xTicks = get(ax, 'xtick');
minY = min(yTicks);

text(xTicks(1), minY - verticalOffset, '$\mathbf{M}$','HorizontalAlignment','center','interpreter', 'latex');


for i=1:nb_approx_MEG
    text(xTicks(i+1), minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');
    
end
title('MATLAB faust');




%% speed-up
mean_Times=mean(timeS);
mean_Times_matlab=mean(timeS_matlab);
dense_matrix_time=mean_Times(1);
real_RCG=dense_matrix_time./mean_Times;
real_RCG_matlab=dense_matrix_time./mean_Times_matlab;
figure,
ax = gca;

set(ax,'xticklabel', [])
set(ax,'Visible','off');
plot(1:nb_approx_MEG,RCG_approxS_MEG,'linewidth',1.5);
set(ax,'XTick',[]);
hold on
plot(1:nb_approx_MEG,real_RCG(2:end),'linewidth',1.5);
plot(1:nb_approx_MEG,real_RCG_matlab(2:end),'linewidth',1.5);
plot(1:nb_approx_MEG,ones(1,nb_approx_MEG),'linewidth',1.5);
verticalOffset=1.01;
for i=1:nb_approx_MEG
    text(i, minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');    
end
legend('theoretical RCG','speed up C++ wrapper faust','speed up MATLAB faust','neutral speed up');
title('speed up');

