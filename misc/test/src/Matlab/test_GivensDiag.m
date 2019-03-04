% This matlab test is a reference for C++ port of GivensFGFT
% the C++ test in misc/test/src/C++/GivensFGFT.cpp.in is using the same file
% so the results can be compared


p = mfilename('fullpath');
[filepath, ~, ~ ] = fileparts(p);

addpath([filepath '../../../../wrapper/matlab/tools'])
load([filepath '/../../../data/mat/test_GivensDiag_Lap_U_J.mat']) % Lap, U, J
%J=5472
%J=57
[facts_givens,D,err,L,choices] = diagonalization_givens(Lap,J);

[U,Dref] = eig(Lap);
norm_U = norm(U,'fro')

for i=1:length(facts_givens)
	norm(full(facts_givens{i}), 'fro');
end

fourier_diag = eye(size(Lap,1));
for j=1:numel(facts_givens)
	fourier_diag = fourier_diag*facts_givens{j};
end

fourier_diag_full = full(fourier_diag);
%fourier_diag_fronorm = norm(fourier_diag_full, 'fro')
spectrum = diag(D);
D_fronorm = norm(D, 'fro')
[Dhat,I] = sort(spectrum);
Uhat = fourier_diag_full(:,I);

err0 = norm(Uhat*full(diag(Dhat))*Uhat' - Lap, 'fro')/norm(Lap, 'fro')
%err2 = norm(Uhat'*Lap*Uhat - diag(diag(Uhat'*Lap*Uhat)),'fro')/norm(Lap,'fro') % equals err0 ?
err1 = norm(U-Uhat, 'fro')/norm(U, 'fro')
Dhat = full(diag(Dhat));

save('~/faust/misc/data/mat/test_GivensDiag_Lap_U_J.mat', '-v7', 'Uhat', 'Dhat', 'choices', 'U', 'Lap', 'J', 'err')
