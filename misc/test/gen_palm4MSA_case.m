clear all;
close all;




filename='cas_test5.mat';
N_FACTS = 4; % Desired number of factors
DIM1=5;
DIM2=8;
DIMS_INTER = round(randi( [min([DIM1+1,DIM2+1]) max([DIM1-1,DIM2-1])] ,[1 ,N_FACTS-1]));
params.nfacts = N_FACTS;
params.niter = 200;
params.verbose = 1;
params.updateway = 0;% parameter are updated from right to left
params.init_lambda = 1; % defaut value equal 1
M=rand(DIM1,DIM2);
params.data=M;



% Constraints
for i=2:N_FACTS-1
if(mod(i,2)==0)
    name_const = 'sp';
else
    name_const = 'sppos'
end
params.cons{1,i} = {name_const,min(DIMS_INTER(i-1:i)),DIMS_INTER(i-1),DIMS_INTER(i)};
params.init_facts{1,i} = eye(DIMS_INTER(i-1),DIMS_INTER(i));
end

params.cons{1,1} = {'sp',min([DIMS_INTER(1),DIM1]),DIM1,DIMS_INTER(1)};
params.init_facts{1,1} = eye(DIM1,DIMS_INTER(1));



params.cons{1,N_FACTS} = {'sp',min([DIMS_INTER(end),DIM2]),DIMS_INTER(N_FACTS-1),DIM2};
params.init_facts{1,N_FACTS} = eye(DIMS_INTER(end),DIM2);




    

    [lambda, facts] = palm4MSA(params);

    
   
    
    



%%profile viewer
%%p = profile('info');
%%profsave(p,['comp_normR_profiler FACT' int2str(N_FACTS) 'DIM : ' int2str(DIMS(1)) ' ' int2str(DIMS(end))]) ;


save(filename);

