close all;
clear all;

colorS={'r','g','b','m'};
suite_nbr_factor=[1,2,4,8];
DIM=2048;
density_max=0.5;
NB_SAMPLE=8;
NB_RUN=100;
n_diff_nbr_factor = length(suite_nbr_factor);
inputfile=['time_comp_NB_FACT_'];
for i=1:n_diff_nbr_factor
    inputfile=[inputfile,int2str(suite_nbr_factor(i)),'_'];   
end
inputfile_root=[inputfile 'DIM_' int2str(DIM) '_density_max_' num2str(density_max)];
inputfile=[inputfile_root '.mat'];

load(inputfile);
t=density_max/NB_SAMPLE:density_max/NB_SAMPLE:density_max;
figure,
hold on
varS_dense = [];
A=[];
legend_var=cell(1,n_diff_nbr_factor+1);
for i=1:n_diff_nbr_factor
    
   var_fact=eval(['t_fact_',int2str(suite_nbr_factor(i))]);
   var_dense=eval(['t_dense_',int2str(suite_nbr_factor(i))]);
   
   p=plot(t,mean(var_fact),colorS{i},'linewidth',2);
   %%plot(repmat(t,NB_RUN,1),var_fact,[colorS{i} '*']);
   %plot(t,mean(var_dense),[colorS{i} '*']);
    varS_dense=cat(3,varS_dense,var_dense);
    A=[A,p];
    legend_var{i}=['nbr facteurs ' int2str(suite_nbr_factor(i))];
end
 p=plot(t,repmat(mean(varS_dense(:)),1,NB_SAMPLE),'k','linewidth',2);
 A=[A,p];
legend_var{end}='produit dense'; 
title(['Prod Faust Vector (DIM ' int2str(DIM) ')']);
xlabel('Density %','fontsize',14,'fontweight','bold');
ylabel('Time (s)','fontsize',14,'fontweight','bold');
h=legend(A,legend_var,'location','northwest','fontsize',14,'fontweight','bold');
set(h,'fontsize',12,'location','northwest');
set(gca,'fontweight','bold');
output_dir='figures/';
im_file=[output_dir inputfile_root '.jpeg'];
print(im_file,'-djpeg');
