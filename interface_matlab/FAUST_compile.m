ROOT_DIR=pwd;
PALM_DIR = [ROOT_DIR filesep '..' filesep 'palm4MSA'];
MAT_DIR = [ROOT_DIR filesep '..' filesep 'matrix'];
FAUST_CORE_DIR = [ROOT_DIR filesep '..' filesep 'faust_core'];
EIGEN_DIR = getenv('EIGEN_ROOT_DIR');

PALM_FILE = dir([PALM_DIR filesep '*.cpp']);
MAT_FILE = dir([MAT_DIR filesep '*.cpp']);
FAUST_CORE_FILE = dir([FAUST_CORE_DIR filesep '*.cpp']);
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


 S1 = 'mexHierarchical_fact.cpp';
SRC_FILE = [S1,PALM_SRC,MAT_SRC,FAUST_CORE_SRC];
INCLUDE = {['-I' EIGEN_DIR],['-I' MAT_DIR],['-I' FAUST_CORE_DIR],['-I' PALM_DIR]};
OTHER_OPT = {'-v'};
INPUT_MEX = [INCLUDE,OTHER_OPT,SRC_FILE];
mex(INPUT_MEX{:});
% V='-v';
% I_EIG=['-I' EIGEN_DIR];
% I_MAT=['-I' MAT_DIR];
% I_FAUST=['-I' FAUST_CORE_DIR];
% I_PALM=['-I' PALM_DIR];
% cell_tab=cell(1,20);
% cell_tab{1}=V;
% cell_tab{2}=I_EIG;
% cell_tab{3}=I_MAT;
% cell_tab{4}=I_FAUST;
% cell_tab{5}=I_PALM;
% cell_tab{6}=S1;
% cell_tab{7}=P1;
% cell_tab{8}=P2;
% cell_tab{9}=P3;
% cell_tab{10}=P4;
% cell_tab{11}=P5;
% cell_tab{12}=P6;
% cell_tab{13}=P7;
% cell_tab{14}=P8;
% cell_tab{15}=M1;
% cell_tab{16}=M2;
% cell_tab{17}=M3;
% cell_tab{18}=M4;
% cell_tab{19}=M5;
% cell_tab{20}=F1;
% %%mex(V,I_EIG,I_MAT,I_FAUST,I_PALM,S1,P1,P2,P3,P4,P5,P6,P7,P8,M1,M2,M3,M4,M5,F1);
% mex(cell_tab{:});


