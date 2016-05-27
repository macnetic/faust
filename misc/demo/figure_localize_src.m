close all
clear all
set_path
load 'results_localise_src'

%% convergence analysis

d1 = cat(4,resDist(:,1,1,:),resDist(:,1,2,:));
d2 = cat(4,resDist(:,2,1,:),resDist(:,2,2,:));
d3 = cat(4,resDist(:,3,1,:),resDist(:,3,2,:));
test2 = 100*[squeeze(d1);zeros(1,1000);squeeze(d2);zeros(1,1000);squeeze(d3)];
%boxPlot(test2')
figure('color',[1 1 1]);
T = bplot(test2','linewidth',1.5);
legend(T)
% title('OMP')
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


%% MATLAB convergence analysis

d1 = cat(4,resDist_matlab(:,1,1,:),resDist_matlab(:,1,2,:));
d2 = cat(4,resDist_matlab(:,2,1,:),resDist_matlab(:,2,2,:));
d3 = cat(4,resDist_matlab(:,3,1,:),resDist_matlab(:,3,2,:));
test2 = 100*[squeeze(d1);zeros(1,1000);squeeze(d2);zeros(1,1000);squeeze(d3)];
%boxPlot(test2')
figure('color',[1 1 1]);
title('MATLAB');
T = bplot(test2','linewidth',1.5);
legend(T)
% title('OMP')
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





%% time comparison

timeS = cat(3,compute_Times(:,1,:),compute_Times(:,2,:),compute_Times(:,3,:));


timeS = squeeze(1000*timeS)';
%boxPlot(test2')
figure('color',[1 1 1]);
T = bplot(timeS,'linewidth',1.5);
legend(T)
% title('OMP')
ylabel('Computed Time (ms)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',x_tick');
set(ax,'xticklabel', [])
axi = axis;

% axis([0 3*(nb_approx_MEG+2) -0.3 4.7])
% xlim([0 3*(nb_approx_MEG]);
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



%% MATLAB time comparison

timeS_matlab = cat(3,compute_Times_matlab(:,1,:),compute_Times_matlab(:,2,:),compute_Times_matlab(:,3,:));


timeS_matlab = squeeze(1000*timeS_matlab)';
%boxPlot(test2')
figure('color',[1 1 1]);
T = bplot(timeS_matlab,'linewidth',1.5);
legend(T)
% title('OMP')
ylabel('Computed Time (ms)')
box on
ax = gca;
x_tick=[1:nb_approx_MEG+1];
mean_x_tick = round(mean(x_tick,2));
set(ax,'XTick',x_tick');
set(ax,'xticklabel', [])
axi = axis;

% axis([0 3*(nb_approx_MEG+2) -0.3 4.7])
% xlim([0 3*(nb_approx_MEG]);
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
title('MATLAB');




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
verticalOffset=1.01;
for i=1:nb_approx_MEG
    text(i, minY - verticalOffset, ['$\widehat{\mathbf{M}}_{' int2str(RCG_approxS_MEG(i)) '}$' ],'HorizontalAlignment','center','interpreter', 'latex');    
end
legend('theoretical RCG','speed up faust','speed up MATLAB');

