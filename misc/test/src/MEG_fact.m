%% test the mexhierarchical_fact function in the MEG configuration

%addpath(['..' filesep 'tools']);
%set_path;

%load config_MEG;
%[mexlambda,mexfact,fc]=launch_hierarchical_fact(params);


function [testPass]=MEG_fact(expectedLambda, expectedLambdaPrecision)

addpath(['..' filesep 'tools']);
set_path;

load config_MEG;
[mexlambda,mexfact,fc]=launch_hierarchical_fact(params);

if (mexlambda >= (expectedLambda - expectedLambdaPrecision)  &&  (mexlambda <= (expectedLambda + expectedLambdaPrecision)) )
    testPass=0; % le résultats est OK
     disp('');
    disp('Test is successful');
else 
    testPass=1; % le resutat n'est pas bon.
    disp('');
    disp('Test is FAILED !!!');
    exit (FAILURE);
end