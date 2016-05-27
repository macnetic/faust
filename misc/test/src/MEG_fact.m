%% test the mexhierarchical_fact function in the MEG configuration

addpath(['..' filesep 'tools']);
set_path;

load config_MEG;
[mexlambda,mexfact,fc]=launch_hierarchical_fact(params);