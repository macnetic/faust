clear all;
close all;
addpath('../build/interface_matlab')
addpath('')
cd '../../Code_Luc'
set_path;
cd '../devcpp/interface_matlab/';

DIM1 = 5;
DIM2 = 2;
data = zeros(DIM1,DIM2);
data(:)=1:DIM1*DIM2;
params.data=data;
params.nfacts = 3;
params.cons{1,1} = {'sp',3,DIM1,DIM1};
params.cons{1,2} = {'sp',4,DIM1,DIM2};
params.cons{1,3} = {'sp',5,DIM2,DIM2};

params.niter = 400;
init_facts=cell(1,params.nfacts);
for i=1:params.nfacts-1
init_facts{i}=eye(params.cons{i}{3},params.cons{i}{4});
end
params.init_facts=init_facts
[mexlambda,mexfact]=mexPalm4MSA(params);
disp(['MEX LAMBDA ' num2str(mexlambda)]);
[lambda,fact]=palm4MSA(params);
disp([' LAMBDA ' num2str(lambda)]);

sum_coeff = 0;
threshold = 0.0000001;
for i=1:params.nfacts
   comp=abs(mexfact{i} - fact{i})> threshold;
   nb_coeff_diff = sum(comp(:));
   sum_coeff = nb_coeff_diff + sum_coeff; 
    
end
disp([' nb coeff different ' int2str(sum_coeff)]);