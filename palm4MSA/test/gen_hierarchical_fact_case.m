clear all;
close all;

filename='cas_hierarchical_factside1.mat';
N_FACTS = 4; % Desired number of factors
DIM1=2;
DIM2=20;
DIMS_INTER = round(randi( [min([DIM1+1,DIM2+1]) max([DIM1-1,DIM2-1])] ,[1 ,N_FACTS-1]));
params.nfacts = N_FACTS;
params.niter1 = 14;
params.niter2 = 15;
params.verbose = 1;
params.fact_side = 1;
params.updateway = 0;% parameter are updated from right to left
params.init_lambda = 100; % defaut value equal 1
M=rand(DIM1,DIM2);
params.data=M;




%% Constraints
if (params.fact_side == 0)
    
    for i=2:N_FACTS-1

        params.cons{1,i} = {'sp',i,DIMS_INTER(i-1),DIMS_INTER(i)};
        params.cons{2,i} = {'sppos',i*i,DIMS_INTER(i),DIM2};

    end

    params.cons{1,1} = {'normcol',0.5,DIM1,DIMS_INTER(1)};
    params.cons{2,1} = {'const',ones(DIMS_INTER(1),DIM2),DIMS_INTER(1),DIM2};

else
    
    for i=N_FACTS-1:-1:2
        
        params.cons{1,N_FACTS-i+1} = {'sppos',i*i,DIM1,DIMS_INTER(i-1)};
        params.cons{2,N_FACTS-i+1} = {'sp',i,DIMS_INTER(i-1),DIMS_INTER(i)};
        

    end

    
    params.cons{1,1} = {'const',ones(DIM1,DIMS_INTER(N_FACTS-1)),DIM1,DIMS_INTER(N_FACTS-1)};
    params.cons{2,1} = {'normcol',0.5,DIMS_INTER(N_FACTS-1),DIM2};
end






    

 [lambda, facts, errors] = hierarchical_fact(params);

    
   
    
    



%%profile viewer
%%p = profile('info');
%%profsave(p,['comp_normR_profiler FACT' int2str(N_FACTS) 'DIM : ' int2str(DIMS(1)) ' ' int2str(DIMS(end))]) ;


save(filename);