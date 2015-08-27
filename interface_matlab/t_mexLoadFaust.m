% nbRow = 5;
% nbCol = 10;
% density = 0.3;
clear all;
close all;
% 
% S=sprand(nbRow,nbCol,density);
% 
S{1}=2*eye(3,5);
S{2}=3*eye(5,7);
S{3}=randint(7,7);

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
mexLoadFaust(S);
