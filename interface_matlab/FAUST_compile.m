ROOT_DIR=pwd;
PALM_DIR = [ROOT_DIR filesep '..' filesep 'palm4MSA'];
MAT_DIR = [ROOT_DIR filesep '..' filesep 'matrix'];
FAUST_CORE_DIR = [ROOT_DIR filesep '..' filesep 'faust_core'];
EIGEN_DIR = getenv('EIGEN_ROOT_DIR');
TOOLS_DIR = 'tools';

PALM_FILE = dir([PALM_DIR filesep '*.cpp']);
MAT_FILE = dir([MAT_DIR filesep '*.cpp']);
FAUST_CORE_FILE = dir([FAUST_CORE_DIR filesep '*.cpp']);
TOOLS_FILE = dir([TOOLS_DIR filesep '*.cpp']); 
PALM_SRC= cell(0,0);
MAT_SRC = cell(1,length(MAT_FILE)-1);
FAUST_CORE_SRC = cell(0,0);

for i=1:length(PALM_FILE)
    PALM_SRC = [PALM_SRC , [PALM_DIR filesep PALM_FILE(i).name]];
end
cpt = 1;
for i=1:length(MAT_FILE)
    if (strcmp(MAT_FILE(i).name,'faust_timer.cpp'))
    else    
        MAT_SRC{cpt} = [MAT_DIR filesep MAT_FILE(i).name];
        cpt = cpt + 1;
    end
end

for i=1:length(FAUST_CORE_FILE)
    FAUST_CORE_SRC = [FAUST_CORE_SRC , [FAUST_CORE_DIR filesep   FAUST_CORE_FILE(i).name]];
end

for i=1:length(TOOLS_FILE)
   TOOLS_SRC{i} =  [TOOLS_DIR filesep TOOLS_FILE(i).name];
 
end

 S1 = ['mex_functions' filesep 'mexHierarchical_fact.cpp'];
SRC_FILE = [S1,PALM_SRC,MAT_SRC,FAUST_CORE_SRC,TOOLS_SRC];
INCLUDE = {['-I' EIGEN_DIR],['-I' MAT_DIR],['-I' FAUST_CORE_DIR],['-I' PALM_DIR],['-I' TOOLS_DIR]};
OTHER_OPT = {'-v'};
INPUT_MEX = [INCLUDE,OTHER_OPT,SRC_FILE];
mex(INPUT_MEX{:});



