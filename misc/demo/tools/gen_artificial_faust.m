%% script pour generer des Faust avec un RCG, une Dimension et un nombre de facteur fixe
% lancer l'executable multiply_comptime ensuite pour effectuer les tests de performances
% puis  display_multiply_comptime pour afficher les resultats
function facts=gen_artficial_faust(Dim,RCG,nb_fact,constraint)
%% cas 1
% RCGs = [2 4 6 8 10];
% Dims = [128 256 512];
% nb_facts = [2,4,6];

%% cas 2

% RCGs = [2 4 8 16];
% Dims = [128 256 512 1024 2048];

% nb_facts = [2,4,8,16];

% RCGs = [2 4];
% Dims = [32 64 128];
% nb_facts = [2,4];


% constraint='sp_row';% choix possible 'sp_row' et 'sp_col_'











if ((strcmp(constraint,'sp') + strcmp(constraint,'sp_row') + strcmp(constraint,'sp_col')) == 0)
    error('the constraint parameter must be equal to sp or sp_row or sp_col');
end


dim1=Dim;
dim2=Dim;
densite_per_fact = 1/(nb_fact*RCG);


facts=cell(1,nb_fact);


if (strcmp(constraint,'sp_row'))
    nb_elt_per_row = round(dim2*densite_per_fact);
    nb_elt = nb_elt_per_row * dim1;
    id_i=reshape(repmat((1:dim1),nb_elt_per_row,1),dim1*nb_elt_per_row,1);
    id_j=zeros(nb_elt,1);
    value=zeros(nb_elt,1);
    
end
if (strcmp(constraint,'sp_col'))
    nb_elt_per_col = round(dim1*densite_per_fact);
    nb_elt = nb_elt_per_col * dim2;
    id_j=reshape(repmat((1:dim2),nb_elt_per_col,1),dim2*nb_elt_per_col,1);
    id_i=zeros(nb_elt,1);
    value=zeros(nb_elt,1);
end
for k=1:nb_fact
    if (strcmp(constraint,'sp'))
        facts{k} = sprand(dim1,dim2,densite_per_fact);
    end
    if (strcmp(constraint,'sp_row'))
        value=rand(nb_elt,1);
        for ll=0:dim1-1
            id_j(ll*nb_elt_per_row+1:(ll+1)*nb_elt_per_row)=randperm(dim2,nb_elt_per_row)';
        end
        facts{k}=sparse(id_i,id_j,value,dim1,dim2);
    end
    if (strcmp(constraint,'sp_col'))
        value=rand(nb_elt,1);
        for ll=0:dim2-1
            id_i(ll*nb_elt_per_col+1:(ll+1)*nb_elt_per_col)=randperm(dim1,nb_elt_per_col)';
        end
        facts{k}=sparse(id_i,id_j,value,dim1,dim2);
    end
    
    
    
    
   
    
    
    
    
    
    
end


end
