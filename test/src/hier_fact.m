%% test the mexhierarchical_fact function in the hier fact configuration

addpath(['..' filesep 'tools']);
set_path;

load config_compared_hierarchical_fact;
[mexlambda,mexfact,fc]=launch_hierarchical_fact(params);