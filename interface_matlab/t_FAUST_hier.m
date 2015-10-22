clear all;
close all;
addpath('../build/mex')
cd '../../Code_Luc';
set_path;
cd '../devcpp/interface_matlab/';

DIM1 = 5;
DIM2 = 2;
data = zeros(DIM1,DIM2);
data(:)=1:DIM1*DIM2;
params.data=data;
params.nfacts = 5;
params.cons{1,1} = {'sp',3,DIM1,DIM1};%;%{'splin',P_REP_tab(i),N_SAMP,DIM};%params.cons{1,1} = %
params.cons{2,1} = {'splin',1,DIM1,DIM2};
params.cons{1,2} = {'spcol',3,DIM1,DIM2};%;%{'splin',P_REP_tab(i),N_SAMP,DIM};%params.cons{1,1} = %
params.cons{2,2} = {'normcol',6,DIM2,DIM2};
params.niter1 = 100;
params.niter2 = 751;
[mexlambda,mexfact]=mexHierarchical_fact(params);
disp(['MEX LAMBDA ' num2str(mexlambda)]);
[lambda,fact]=hierarchical_fact(params);
disp([' LAMBDA ' num2str(lambda) ' MEX_LAMBDA ' num2str(mexlambda)]);

sum_coeff = 0;
threshold = 0.00001;
disp('\n\n COMPARAISON ');
for i=1:params.nfacts
   comp=abs(mexfact{i} - fact{i})> threshold;
   nb_coeff_diff = sum(comp(:));
   disp([int2str(i) ' nb_coeff diff : ' int2str(nb_coeff_diff)]);
   sum_coeff = nb_coeff_diff + sum_coeff; 
    
end
disp([' nb coeff different ' int2str(sum_coeff)]);


mex_error = norm(params.data- lambda*dvp(fact));
error = norm(params.data - mexlambda*dvp(mexfact));