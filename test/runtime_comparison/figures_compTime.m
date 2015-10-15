function f1=figures_compTime(M, n, densite_or_nnzrowcol, N_run_inside, t_dense_matlab, t_fact_matlab, t_dense_cpp, t_fact_cpp)

densite = zeros(size(M))-1.0;
nz = zeros(size(M))-1.0;
for k=1:numel(M)
    if densite_or_nnzrowcol(k) <= 1
        nz(k) = round(densite_or_nnzrowcol(k)*n(k));
        densite(k) = nz(k)/n(k);
    else
        nz(k) = densite_or_nnzrowcol(k);
        densite(k) = nz(k)/n(k);
    end
end

f1=figure('color',[1,1,1],'Visible','off');

%Marker Size
MS1 = 10;
MS2 = 20;

%courbe matlab
nb_pl = 8;
if size(unique(n)) == size(n)
    abscisse = n;
    abscisse_str = 'dimension n';
    data.ind_xlabel=1;
elseif size(unique(nz)) == size(nz)
    abscisse = nz;
    abscisse_str = 'nnz par ligne/colonne';
    data.ind_xlabel=2;
elseif size(unique(M)) == size(M)
    abscisse = M;
    abscisse_str = 'nombre de facteurs';
    data.ind_xlabel=3;
end


ax1=subplot(2,1,1, 'Position', [0.2 0.55 0.75 0.4]);
pl = cell(nb_pl,1);
pl{1} = loglog(abscisse,t_dense_matlab,'c.','markersize',MS1); hold on;
pl{2} = loglog(abscisse,t_dense_cpp,'.','markersize',MS1); 
set(pl{2}, 'Color', [237, 127, 16]/255); %orange
 
pl{3} = loglog(abscisse,t_fact_matlab,'m.','markersize',MS1); 
pl{4} = loglog(abscisse,t_fact_cpp,'y.','markersize',MS1); 

pl{5} = loglog(abscisse,mean(t_dense_matlab,1),'b.-','markersize',MS2);
pl{6} = loglog(abscisse,mean(t_dense_cpp,1),'k.-','markersize',MS2);

pl{7} = loglog(abscisse,mean(t_fact_matlab,1),'r.-','markersize',MS2);
pl{8} = loglog(abscisse,mean(t_fact_cpp,1),'g.-','markersize',MS2);

abscisse_sort = sort(abscisse);
set(ax1, 'XTick', abscisse_sort);
xTickLabels = cell(4,numel(n));
xLabels = {'Dimension n','nnz par col/ligne','nbre de facteurs','densite'};

for k=1:numel(n)
    idx = find(abscisse_sort(k)==abscisse);
    xTickLabels{1,k} = sprintf('%d',n(idx));
    xTickLabels{2,k} = sprintf('%d',nz(idx));
    xTickLabels{3,k} = sprintf('%d',M(idx));
    xTickLabels{4,k} = sprintf('%.2e',densite(idx));
end
set(ax1, 'XTickLabel', xTickLabels(data.ind_xlabel,:));
xlabel(ax1, xLabels{data.ind_xlabel},'fontsize',14);

set(ax1,'UserData',data);

ylabel('Time (s)','fontsize',14);
%title('Vector multiplication runtime comparison')
warning('off')
h1=legend([pl{1}(1),pl{2}(1), pl{3}(1),pl{4}(1), pl{5}(1),pl{6}(1), pl{7}(1),pl{8}(1)],...
    '(1):Dense matlab','(2):Dense cpp',...
    '(3):Faust matlab','(4):Faust cpp',...
    '(5):Dense matlab (mean)','(6):Dense cpp (mean)',...
    '(7):Faust matlab (mean)','(8):Faust cpp (mean)',...
    'location','northwest','fontsize',12);
warning('on')
set(h1,'fontsize',10,'location','southeast');
grid

%si toutes les densites sont identiques
titre=['Nb run = ' num2str(size(t_dense_matlab,1)) ' ; Nb run inside = ' num2str(N_run_inside)];
if min(densite_or_nnzrowcol==densite_or_nnzrowcol(1)) ...
        && densite_or_nnzrowcol(1)<=1
   titre=[titre ' ; densite demandee = ' num2str(densite_or_nnzrowcol(1))];
%si tous les nz sont identiques
elseif min(densite_or_nnzrowcol==densite_or_nnzrowcol(1)) ...
        && densite_or_nnzrowcol(1)>1
   titre=[titre ' ; nnz el par col/ligne = ' num2str(densite_or_nnzrowcol(1))];
end
if min(M==M(1))
   titre=[titre ' ; nbre de facteurs = ' num2str(M(1))];
end
if min(n==n(1))
   titre=[titre ' ; dimension = ' num2str(n(1))];
end
title(titre);

group = uibuttongroup(f1,'Position',[0.01 0.55 0.09 0.4]);
check = zeros(nb_pl,1)-1;
% espace entre deux bouttons check
a = 0.02;
%hauteur du check button
b = (1-a*(nb_pl+1))/nb_pl;
for k=1:nb_pl
    check(k) = uicontrol(group, 'Style','check', ...
        'String',['(' num2str(k) ')'], ...
        'Units','normalized' , 'Position',[0.1 (nb_pl-k)*(a+b)+a 0.8 b], 'Min',0,'Max',k,'Value',k);
end

setappdata(f1,'xTickLabels',xTickLabels);
setappdata(f1,'xLabels',xLabels);
setappdata(f1,'pl',pl);
setappdata(f1,'check',check);
setappdata(f1,'axe',ax1);
setappdata(f1,'legende',h1);
set(check,'Callback',@push_callback);
set(ax1,'ButtonDownFcn',@ax_callback);

RCtheo = nz.*M./n;

%courbe theorique
ax2=subplot(2,1,2, 'Position', [0.2 0.05 0.75 0.4]);
l1 = loglog(abscisse,1./RCtheo,'k.-');hold on
%courbe matlab
RCprac_matlab = mean(t_fact_matlab,1)./mean(t_dense_matlab,1);
l21 = loglog(abscisse,1./RCprac_matlab,'r.-');
%courbe cpp
RCprac_cpp = mean(t_fact_cpp,1)./mean(t_dense_cpp,1);
l22 = loglog(abscisse,1./RCprac_cpp,'g.-');

leg = legend([l1(1),l21(1),l22(1)],'RCG','$\widehat{\textrm{RCG\ matlab}}$','$\widehat{\textrm{RCG\ cpp}}$');
set(leg,'fontsize',10,'location','northwest','interpreter','latex');

set(ax2, 'XTick', sort(abscisse));
set(ax2, 'XTickLabel', xTickLabels(data.ind_xlabel,:));
xlabel(ax2,xLabels{data.ind_xlabel},'fontsize',14);
set(ax2,'UserData',data);
ylabel(ax2,'Computational gain','fontsize',14);%\frac{\text{Time gain}}{\text{RC}} 
grid

set(ax2,'ButtonDownFcn',@ax_callback);

function push_callback(hObj, eventData)
pl = getappdata(get(get(hObj,'Parent'),'Parent'),'pl');
% si le checkbox vient d etre decoche
if get(hObj,'Value') == get(hObj,'Min')
    set(pl{get(hObj,'Max')},'Visible','off')
% si le checkbox vient d etre coche
elseif get(hObj,'Value') == get(hObj,'Max')
    set(pl{get(hObj,'Max')},'Visible','on')
else
    error('incoherence entre proprietes Min, Max et Value du checkbox');
end

function ax_callback(hObj, eventData)
data = get(hObj,'UserData');
set(hObj,'UserData',data);
xTickLabels = getappdata(get(hObj,'Parent'),'xTickLabels');
xLabels = getappdata(get(hObj,'Parent'),'xLabels');
set(hObj, 'XTickLabel', xTickLabels(data.ind_xlabel+1,:));
xlabel(hObj,xLabels{data.ind_xlabel+1},'fontsize',14);
data.ind_xlabel=mod(data.ind_xlabel+1,4);
set(hObj,'UserData',data);
