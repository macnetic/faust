clear 
close all
clc
load('results_comptime_c++_upto_8192');
t_dense_c = t_dense;
t_fact_c = t_fact;
clear t_dense;
clear t_fact;
load('results_comptime_matlab_upto8192');
M = 1:size(t_dense,2);
n = 2.^M;

figure('color',[1,1,1])
subplot_tight(2,1,1,0.08)
p1 = loglog(n,t_dense_c,'r.','markersize',8); p1 = p1(1);
hold on
p5 = loglog(n,mean(t_dense_c),'r','markersize',20);
p3 = loglog(n,t_fact_c,'m.','markersize',8); p3 = p3(1);
p6 = loglog(n,mean(t_fact_c),'m','markersize',20,'linewidth',2);

p11 = loglog(n,t_dense,'c.','markersize',8); p11 = p11(1);
p15 = loglog(n,mean(t_dense),'c.-','markersize',20);
p13 = loglog(n,t_fact_c,'b.','markersize',8); p13 = p13(1);
p16 = loglog(n,mean(t_fact),'b.-','markersize',20,'linewidth',2);




axis([1.9 17800 0 1])
ax=gca;
set(ax, 'XTick', 2.^M);
%title('Vector multiplication runtime comparison')
xlabel('Dimension n','fontsize',14,'fontweight','bold');
ylabel('Time (s)','fontsize',14,'fontweight','bold');
h=legend([p1,p5,p3,p6,p11,p15,p13,p16],'Dense multiplication C++','Dense multiplication (mean) C++','FAUST multiplication C++','FAUST multiplication (mean) C++','Dense multiplication matlab','Dense multiplication (mean) matlab','FAUST multiplication matlab','FAUST multiplication (mean) matlab','location','northwest','fontsize',14);
set(h,'fontsize',12,'location','northwest','fontweight','bold');
set(gca,'fontweight','bold');

subplot_tight(2,1,2,0.08)
RCtheo = 2*M./n;
RCprac_C = mean(t_fact_c)./mean(t_dense_c);
RCprac = mean(t_fact)./mean(t_dense);
l1 = loglog(n,1./RCprac_C,'r.-');
hold on
l2 = loglog(n,1./RCprac,'b.-');
l3 = loglog(n,1./RCtheo,'k.-');
leg = legend([l1,l2,l3],'$\widehat{\textrm{RCG C++}}$','$\widehat{\textrm{RCG Matlab}}$','RCG');
set(leg,'fontsize',12,'location','northwest','interpreter','latex','fontweight','bold');
axis([1.9 17800 0.05 1000])%axis([1.9 17800 0 35])
xlabel('Dimension n','fontsize',14,'fontweight','bold');
ylabel('Computational gain','fontsize',14,'fontweight','bold');%\frac{\text{Time gain}}{\text{RC}} 
ax=gca;
set(ax, 'XTick', 2.^M,'fontweight','bold');

