% This matlab test is a reference for C++ port of EigTJParallel
% the C++ test in misc/test/src/C++/EigTJParallel.cpp.in is using the same file
% so the results can be compared

p = mfilename('fullpath');
[filepath, ~, ~ ] = fileparts(p);
load([filepath, '/../../../data/mat/test_GivensDiagParallel_Lap_U_J_choices.mat']) % Lap, U, J
ref_choices = choices;


[facts_givens,D,err,L,choices] = diagonalization_givens_parall(Lap,J,size(Lap,1)/2);

norm_U = norm(U,'fro')

for i=1:length(facts_givens)
	norm(full(facts_givens{i}), 'fro')
end

fourier_diag = eye(size(Lap,1));
for j=1:numel(facts_givens)
	fourier_diag = fourier_diag*facts_givens{j};
end

fourier_diag_full = full(fourier_diag);
%fourier_diag_fronorm = norm(fourier_diag_full, 'fro')
spectrum = diag(D);
D_fronorm = norm(D, 'fro')
[sorted_spectrum,I] = sort(spectrum);
%sorted_spectrum
Uhat_givens = fourier_diag_full(:,I);

disp("Error: norm(Uhat*Dhat*Uhat'- Lap, 'fro')/norm(Lap, 'fro')")
err0 = norm(Uhat_givens*full(diag(sorted_spectrum))*Uhat_givens' - Lap, 'fro')/norm(Lap, 'fro')
%err2 = norm(Uhat_givens'*Lap*Uhat_givens - diag(diag(Uhat_givens'*Lap*Uhat_givens)),'fro')/norm(Lap,'fro') % equals err0
%
err1 = norm(U-Uhat_givens, 'fro')/norm(U, 'fro')

disp("Verifying that reference choices for pivots are respected by this exec.")
disp("ref_choices == choices:")
all(all(ref_choices == choices))
Dhat = full(diag(sorted_spectrum));
Uhat = Uhat_givens
save([filepath '/../../../data/mat/test_GivensDiagParallel_Lap_U_J_choices.mat'], '-v7', 'Uhat', 'Dhat', 'choices', 'U', 'Lap', 'J', 'err', 'err0', 'err1')
