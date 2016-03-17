function genere_params(dim1,dim2,nfacts)

nb_row = max(dim1,dim2);
nb_col = min(dim1,dim2);


params.data = rand(nb_row,nb_col);
params.nfacts = nfacts;

params.cons{1,1}={'splin', 15, nb_row, nb_col};
for k=2:nfacts-1
    params.cons{1,k}={'sp', nb_row*nb_row/25.5, nb_row, nb_row};
end
for k=1:nfacts-1
    params.cons{2,k}={'sp', nb_row*nb_row*(1-(0.08*(k-1))) , nb_row, nb_row};
end

params.niter1 = 200;
params.niter2 = 200;
params.verbose = 0;

save(sprintf('%s/config_%drows_%dcols_%dfacts.mat',pwd(),nb_row,nb_col,nfacts));



