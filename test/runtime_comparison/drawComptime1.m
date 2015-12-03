%% this script must be run after devcpp/test/gen_artificial_faust.m
% which some faust example of Faust with different dimensions RCG and 
% number of factors, after run
% Z:\Devel\Faust\trunk\devcpp\build\testing\bin\comptime1 which some time
% comparison between dense mat-vector product and faust-vector product
%% this script displays the resulting RCG (speed-up) between the 2 different
% multiplication (faust vs dense matrix) for different RCG
close all;
clear all;

opt_subplot = 1;
liste_prefix = {'sp_','sp_col_','sp_row_'};
prefix_name = {'sp','spcol','sprow'};
nbr_test=length(liste_prefix);
nbr_essai=2;
inputdir=['..' filesep '..' filesep 'build' filesep 'testing' filesep 'output' filesep];
% inputdir=[ 'output_rowmajor' filesep];
for i=1:nbr_test
    tdense=[];
    tfaust=[];
    prefix = liste_prefix{i};
    inputdir1=[inputdir prefix];
    for h=1:nbr_essai
        
        tdense_h=load([ inputdir1 'temps_dense' int2str(h) '.dat']);
        tfaust_h=load([ inputdir1 'temps_faust' int2str(h) '.dat']);
        DIMS=load([ inputdir1 'DIMS' int2str(h) '.dat']);
        RCGS=load([ inputdir1 'RCGS' int2str(h) '.dat']);
        NB_FACTS=load([ inputdir1 'NB_FACTS' int2str(h) '.dat']);
        nDIMS=length(DIMS);
        nRCGS=length(RCGS);
        nNB_FACTS=length(NB_FACTS);
        tdense=cat(3,tdense,tdense_h);
        tfaust=cat(3,tfaust,tfaust_h);
    end
[dim1,dim2,dim3]=size(tdense);
mean_tfaust_int=mean(tfaust,3);
mean_tdense_int=mean(tdense,3);
mean_tfaust_int=mean(mean_tfaust_int);
mean_tdense_int=mean(mean_tdense_int);
if (length(mean_tfaust_int) ~= (nDIMS*nRCGS*nNB_FACTS))
   error('incompatible file'); 
end
mean_tfaust{i} = reshape(mean_tfaust_int,nRCGS,nDIMS,nNB_FACTS);
mean_tdense_int = reshape(mean_tdense_int,nRCGS,nDIMS,nNB_FACTS);
% mean_tdense{i} = repmat(mean(mean_tdense_int),nRCGS,1,1);
mean_tdense{i} = repmat(mean(mean(mean_tdense_int),3),nRCGS,1,nNB_FACTS);






% if (opt_subplot == 0)
%     for ll=1:length(NB_FACTS)
%         figure,
%         imagesc(RCG(:,:,ll));
%         colorbar;
%         title(['RCG ' prefix 'NB_facts' NB_FACTS(ll) ]);
%     end
%     
%     
% else
%     figure,
%     for ll=1:length(NB_FACTS)
%         subplot(1,ll,NB_FACTS);
%         imagesc(RCG(:,:,ll));
%         colorbar;
%         title(['RCG ' prefix 'NB_facts' NB_FACTS(ll) ]);
%     end
% end




real_RCG{i} = mean_tdense{i}./mean_tfaust{i};
color_axe=[min(real_RCG{i}(:)) max(real_RCG{i}(:))];


figure,

if (round(sqrt(nNB_FACTS)) == sqrt(nNB_FACTS))
    m_subplot = sqrt(nNB_FACTS);
    n_subplot = sqrt(nNB_FACTS);
else
    m_subplot = 1;
    n_subplot = nNB_FACTS; 
end

for ll=1:nNB_FACTS
    
subplot(m_subplot,n_subplot,ll);
imagesc(log(DIMS)./log(2),log(RCGS)/log(2),real_RCG{i}(:,:,ll));
title(['n_fact=' int2str(NB_FACTS(ll)) ' ' prefix]);
caxis(color_axe);

% set(gca,'YTick',RCGS);
% set(gca,'XTick',DIMS);
xlabel('log(DIM)','FontWeight','bold');
ylabel('log(theo RCG)','FontWeight','bold');
if (ll==nNB_FACTS)
    colorbar;
end
% bool_RCG = (real_RCG> 2);
% 
% figure,
% imagesc(bool_RCG);
% title('real RCG');
% colorbar;

end

max_speed_up = max(real_RCG{i}(:));
mean_speed_up = mean(real_RCG{i}(:));
min_speed_up = min(real_RCG{i}(:));

disp([prefix 'SPEED_UP : max ' num2str(max_speed_up) ' moy ' num2str(mean_speed_up) ' min ' num2str(min_speed_up) ]); 







end

% figure,
% for i=1:nbr_test
% loglog(DIMS,squeeze(mean_tdense{i}(1,:,1)));
% hold on
% end
% title('t_dense comp');
% legend(liste_prefix);


for i=1:nbr_test
    comp=mean_tfaust{i}./mean_tfaust{mod((i+1)-1,nbr_test)+1};
    color_axe=[min(comp(:)) max(comp(:))];
            figure,
    for ll=1:nNB_FACTS

        subplot(m_subplot,n_subplot,ll);
        imagesc(log(DIMS)./log(2),log(RCGS)/log(2),comp(:,:,ll));
        title(['nb fact=' int2str(NB_FACTS(ll))],'fontsize',10);
        caxis(color_axe);
        
        % set(gca,'YTick',RCGS);
        % set(gca,'XTick',DIMS);
%         xlabel('log(DIM)','FontWeight','bold');
%         ylabel('log(theo RCG)','FontWeight','bold');
        if (ll==nNB_FACTS)
            colorbar;
        end
    end
    %% pour avoir acces a suplabel
    addpath('Z:\Devel\Faust\trunk\matlab_tools');
    [ax,h1]=suplabel(['tps ' prefix_name{i} '/' prefix_name{mod((i+1)-1,nbr_test)+1}] ,'t');
    [ax,h2]=suplabel('log(DIM)' ,'x');
    [ax,h3]=suplabel('log(RCG)' ,'y');
end











