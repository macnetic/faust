close all;
clear all;

tdense=load('temps_dense_Faust_example.dat');
tfaust=load('temps_faust_Faust_example.dat');
DIMS=load('Faust_example.mat','Dims');
RCGS=load('Faust_example.mat','RCGs');
DIMS=DIMS.Dims;
RCGS=RCGS.RCGs;
nDIMS=length(DIMS);
nRCGS=length(RCGS);
[dim1,dim2]=size(tdense);
mean_tfaust=mean(tfaust);
mean_tdense=mean(tdense);
if (length(mean_tfaust) ~= (nDIMS*nRCGS))
   error('incompatible file'); 
end
mean_tfaust = reshape(mean_tfaust,nRCGS,nDIMS);
mean_tdense = reshape(mean_tdense,nRCGS,nDIMS);
mean_tdense = repmat(mean(mean_tdense),nRCGS,1);

% imagesc(mean_tfaust);
% colorbar;
figure,
imagesc(mean_tdense);
colorbar;










real_RCG = mean_tdense./mean_tfaust

figure,
imagesc(real_RCG);
title('real RCG');
colorbar;

bool_RCG = (real_RCG> 2);

figure,
imagesc(bool_RCG);
title('real RCG');
colorbar;


