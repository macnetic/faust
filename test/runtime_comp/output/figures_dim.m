close all;
clear all;

colorS={'r','g','b','m'};
nbr_factor=2;
DIM_dep=32;
DIM_fin=8192;
nbr_dim=(log(DIM_fin)-log(DIM_dep))/log(2)+1;
density_max=0.5;
NB_SAMPLE=8;
NB_RUN=100;
inputfile_root=['time_comp_DIM_' int2str(DIM_dep) '_' int2str(DIM_fin)];
inputfile_root=[inputfile_root '_NB_FACT_' int2str(nbr_factor) '_density_max_' num2str(density_max)];
inputfile=[inputfile_root '.mat'];

load(inputfile);
densities=density_max/NB_SAMPLE:density_max/NB_SAMPLE:density_max;
varS_dense = [];
var_fact=[];
var_dense=[];
log_DIMS=[log(DIM_dep)/log(2):1:log(DIM_fin)/log(2)];
DIMS=2.^log_DIMS;

for i=1:nbr_dim
    
   var_fact_tmp=eval(['t_fact_',int2str(DIMS(i))]);
   var_dense_tmp=eval(['t_dense_',int2str(DIMS(i))]);
   
   var_fact = cat(3,var_fact,var_fact_tmp);
   var_dense = cat(3,var_dense,var_dense_tmp);
   


end



var_fact_mean=squeeze(mean(var_fact));
var_dense_mean=repmat(squeeze(mean(mean(var_dense)))',NB_SAMPLE,1);


figure,
imagesc(log_DIMS,densities,var_fact_mean);
colorbar;
title(['Prod Faust Vector (nbr_factor ' int2str(nbr_factor) ')']);
set(gca,'YTick',densities);
set(gca,'Xscale','log');
xlabel('log(DIM)/log(2) ','fontsize',14);
ylabel('density % ','fontsize',14);



figure,
imagesc(log_DIMS,densities,var_dense_mean);
colorbar;
title(['Prod Dense Mat Vector']);
set(gca,'XTick',DIMS);
set(gca,'YTick',densities);
xlabel('log(DIM)/log(2) ','fontsize',14);
ylabel('density % ','fontsize',14);


tps_fact_min=(var_fact_mean<var_dense_mean);


figure,
imagesc(log_DIMS,densities,tps_fact_min);
colorbar;
title(['Comp Faust-vector vs Dense mat-vector']);
%set(gca,'XTick',DIMS);
set(gca,'YTick',densities);
xlabel('log(DIM)/log(2) ','fontsize',14);
ylabel('density % ','fontsize',14);










% xlabel('Density %','fontsize',14);
% ylabel('Time (s)','fontsize',14);
% h=legend(A,legend_var,'location','northwest','fontsize',14);
% set(h,'fontsize',12,'location','northwest');
% 
% output_dir='figures/';
% im_file=[output_dir inputfile_root '.jpeg'];
% print(im_file,'-djpeg');
