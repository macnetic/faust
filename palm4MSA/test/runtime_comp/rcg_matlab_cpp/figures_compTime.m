function figures_compTime(M, n, densite_or_nnzperroworcol, t_dense_matlab, t_fact_matlab, t_dense_cpp, t_fact_cpp)

densite = zeros(size(M))-1.0;
nz = zeros(size(M))-1.0;
for k=1:numel(M)
    if densite_or_nnzperroworcol(k) < 1
        densite(k) = densite_or_nnzperroworcol(k);
        nz(k) = round(densite(k)*n(k));
    else
        nz(k) = densite_or_nnzperroworcol(k);
        densite(k) = nz(k)/n(k);
    end
end

f1=figure('color',[1,1,1]);%,'Visible','off');

%Marker Size
MS1 = 10;
MS2 = 20;

%courbe matlab
nb_pl = 8;
ax1=subplot(2,1,1, 'Position', [0.2 0.55 0.75 0.4]);
pl = cell(nb_pl,1);
pl{1} = loglog(n,t_dense_matlab,'c.','markersize',MS1); hold on;
pl{2} = loglog(n,t_dense_cpp,'.','markersize',MS1); 
set(pl{2}, 'Color', [237, 127, 16]/255); %orange
 
pl{3} = loglog(n,t_fact_matlab,'m.','markersize',MS1); 
pl{4} = loglog(n,t_fact_cpp,'y.','markersize',MS1); 

pl{5} = loglog(n,mean(t_dense_matlab),'b.-','markersize',MS2);
pl{6} = loglog(n,mean(t_dense_cpp),'k.-','markersize',MS2);

pl{7} = loglog(n,mean(t_fact_matlab),'r.-','markersize',MS2);
pl{8} = loglog(n,mean(t_fact_cpp),'g.-','markersize',MS2);

set(ax1, 'XTick', n);
xTickLabels = cell(4,numel(n));
xLabels = {'Dimension n','nbre de facteurs','nnz par col/ligne','densite'};

for k=1:numel(n)
    xTickLabels{1,k} = sprintf('%d',n(k));
    xTickLabels{2,k} = sprintf('%d',M(k));
    xTickLabels{3,k} = sprintf('%d',nz(k));
    xTickLabels{4,k} = sprintf('%.2e',densite(k));
end
set(ax1, 'XTickLabel', xTickLabels(1,:));
xlabel(ax1, xLabels{1},'fontsize',14);
data.ind_xlabel=1;
set(ax1,'UserData',data);

ylabel('Time (s)','fontsize',14);
%title('Vector multiplication runtime comparison')
warning('off')
h1=legend([pl{1}(1),pl{2}(1), pl{3}(1),pl{4}(1), pl{5}(1),pl{6}(1), pl{7}(1),pl{8}(1)],...
    '(1):Dense multiplication matlab','(2):Dense multiplication cpp',...
    '(3):Multi-layer sparse multiplication matlab','(4):Multi-layer sparse multiplication cpp',...
    '(5):Dense multiplication (mean) matlab','(6):Dense multiplication (mean) cpp',...
    '(7):Multi-layer sparse multiplication (mean) matlab','(8):Multi-layer sparse multiplication (mean) cpp',...
    'location','northwest','fontsize',12);
warning('on')
set(h1,'fontsize',12,'location','southeast');
grid

%si toutes les densites sont identiques
if min(densite==densite(1))
   title(['Nb run = ' num2str(size(t_dense_matlab,1)) ' ; densite de chaque facteur = ' num2str(densite(1))])
%si tous les nz sont identiques
elseif min(nz==nz(1))
   title(['Nb run = ' num2str(size(t_dense_matlab,1)) ' ; nnz el par col/ligne = ' num2str(nz(1))])
else
    title(['Nb run = ' num2str(size(t_dense_matlab,1))])
end

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
l1 = loglog(n,1./RCtheo,'k.-');hold on
%courbe matlab
RCprac_matlab = mean(t_fact_matlab)./mean(t_dense_matlab);
l21 = loglog(n,1./RCprac_matlab,'r.-');
%courbe cpp
RCprac_cpp = mean(t_fact_cpp)./mean(t_dense_cpp);
l22 = loglog(n,1./RCprac_cpp,'g.-');

leg = legend([l1,l21,l22],'RCG','$\widehat{\textrm{RCG\ matlab}}$','$\widehat{\textrm{RCG\ cpp}}$');
set(leg,'fontsize',12,'location','northwest','interpreter','latex');

set(ax2, 'XTick', n);
set(ax2, 'XTickLabel', xTickLabels(1,:));
xlabel(ax2,xLabels{1},'fontsize',14);
data.ind_xlabel=1;
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
data.ind_xlabel=mod(data.ind_xlabel+1,3);
set(hObj,'UserData',data);
