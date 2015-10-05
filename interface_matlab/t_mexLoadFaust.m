% nbRow = 5;
% nbCol = 10;
% density = 0.3;
clear all;
close all;
% 
% S=sprand(nbRow,nbCol,density);
%
addpath('../build/interface_matlab');
addpath('tools/');
addpath([getenv('MKLDIR') '/lib/intel64']);
addpath([getenv('MKL_COMPILER_DIR') '/lib/intel64']);

getenv('LD_LIBRARY_PATH')
%setenv('LD_LIBRARY_PATH',[getenv('')])
S{1}=2*eye(3,5);
% S{2}=3*eye(5,7);
% S{3}=randint(7,7);

% Sparrow = zeros(7,3);
% Sparrow(2,1) = 1;
% Sparrow(5,1) = 1;
% Sparrow(3,2) = 1;
% Sparrow(2,3) = 2;
% Sparrow(5,3) = 1;
% Sparrow(6,3) = 1;
% Sparrow = sparse(Sparrow);
% S=Sparrow;
% S=sparse(S);
fc=matlab_faust(S);
PROD=S{1};

for i=2:length(S)
	PROD=PROD*S{i};
end

PROD
x=ones(size(S{end},2),1);
y=fc*x;
y
