clear all;
close all;

DIM1 = 5;
DIM2 = 2;
data = zeros(DIM1,DIM2);
data(:)=1:DIM1*DIM2;
params.data=data;
params.nfacts = 3;
params.cons{1,1} = {'sp',3,DIM1,DIM1};%;%{'splin',P_REP_tab(i),N_SAMP,DIM};%params.cons{1,1} = %
params.cons{2,1} = {'sp',4,DIM1,DIM2};
params.cons{1,2} = {'sp',5,DIM1,DIM2};%;%{'splin',P_REP_tab(i),N_SAMP,DIM};%params.cons{1,1} = %
params.cons{2,2} = {'sp',6,DIM2,DIM2};
params.niter1 = 100;
params.niter2 = 751;
[mexlambda,mexfact]=mexHierarchical_fact(params);
disp(['MEX LAMBDA ' num2str(mexlambda)]);
[lambda,fact]=hierarchical_fact(params);
disp([' LAMBDA ' num2str(lambda)]);

sum_coeff = 0;
threshold = 0.0000001;
for i=1:params.nfacts
   comp=abs(mexfact{i} - fact{i})> threshold;
   nb_coeff_diff = sum(comp(:));
   sum_coeff = nb_coeff_diff + sum_coeff; 
    
end
disp([' nb coeff different ' int2str(sum_coeff)]);