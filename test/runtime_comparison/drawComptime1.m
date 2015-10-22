close all;
clear all;
inputdir=['..' filesep 'output' filesep];
tdense=load([ inputdir 'temps_dense.dat']);
tfaust=load([ inputdir 'temps_faust.dat']);
DIMS=load([ inputdir 'DIMS.dat']);
RCGS=load([ inputdir 'RCGS.dat']);
NB_FACTS=load([ inputdir 'NB_FACTS.dat']);
nDIMS=length(DIMS);
nRCGS=length(RCGS);
nNB_FACTS=length(NB_FACTS);
[dim1,dim2]=size(tdense);
mean_tfaust=mean(tfaust);
mean_tdense=mean(tdense);
if (length(mean_tfaust) ~= (nDIMS*nRCGS*nNB_FACTS))
   error('incompatible file'); 
end
mean_tfaust = reshape(mean_tfaust,nRCGS,nDIMS,nNB_FACTS);
mean_tdense = reshape(mean_tdense,nRCGS,nDIMS,nNB_FACTS);
mean_tdense = repmat(mean(mean_tdense),nRCGS,1,1);















real_RCG = mean_tdense./mean_tfaust;
color_axe=[min(real_RCG(:)) max(real_RCG(:))];



for i=1:nNB_FACTS
    
figure,
imagesc(log(DIMS)./log(2),RCGS,real_RCG(:,:,i));
title(['real RCG (nb fact = ' int2str(NB_FACTS(i)) ')']);
caxis(color_axe);
colorbar;
set(gca,'YTick',RCGS);
xlabel('log(DIM)','FontWeight','bold');
ylabel('theo RCG','FontWeight','bold');

% bool_RCG = (real_RCG> 2);
% 
% figure,
% imagesc(bool_RCG);
% title('real RCG');
% colorbar;

end
