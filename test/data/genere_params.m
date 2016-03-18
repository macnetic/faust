function genere_params(dim1,dim2,nfacts)

nb_row = max(dim1,dim2);
nb_col = min(dim1,dim2);

rho = 0.7;

params.data = rand(nb_row,nb_col);
params.nfacts = nfacts;

params.cons{1,1}={'splin', round(2+nb_col/15), nb_row, nb_col};
for k=2:nfacts-1
    %params.cons{1,k}={'sp', round(nb_col*nb_col/25.5), nb_col, nb_col};
    params.cons{1,k}={'sp', round(rho^(nfacts)*nb_col*nb_col), nb_col, nb_col};
end
for k=1:nfacts-1
    %params.cons{2,k}={'sp', round(nb_col*nb_col*(1-(0.08*(k-1)))) , nb_col, nb_col};
    params.cons{2,k}={'sp', round(rho^(k-1)*nb_col*nb_col) , nb_col, nb_col};
end

params.niter1 = 200;
params.niter2 = 200;
params.verbose = 0;

save(sprintf('%s/config_%drows_%dcols_%dfacts.mat',pwd(),nb_row,nb_col,nfacts));



