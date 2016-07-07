%% Description runtime_comparison.m
%  Runtime comparison
%
%  This script performs runtime comparison between faust multiplication
%  and dense matrix multiplication for various configuration of faust
%  (dimension of the faust, number of factors, Relative Complexity Gain
%  (RCG), fix type of sparsity (sp, spcol,splin))
%
%	
%
% For more information on the FAuST Project, please visit the website of 
% the project :  <http://faust.gforge.inria.fr>
%
%% License:
% Copyright (2016):	Luc Le Magoarou, Remi Gribonval
%			INRIA Rennes, FRANCE
%			http://www.inria.fr/
%
% The FAuST Toolbox is distributed under the terms of the GNU Affero 
% General Public License.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published 
% by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
% See the GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% Contacts:	
%   Nicolas Bellot : nicolas.bellot@inria.fr
%   Leman Adrien   : adrien.leman@inria.fr
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
%	approximations of matrices and applications", Journal of Selected 
%	Topics in Signal Processing, 2016.
%	<https://hal.archives-ouvertes.fr/hal-01167948v1>
%
% [2]   A. Gramfort, M. Luessi, E. Larson, D. Engemann, D. Strohmeier, 
%	C. Brodbeck, L. Parkkonen, M. Hamalainen, MNE software for processing
%	MEG and EEG data <http://www.ncbi.nlm.nih.gov/pubmed/24161808>, 
%	NeuroImage, Volume 86, 1 February 2014, Pages 446-460, ISSN 1053-8119, 
%	[DOI] <http://dx.doi.org/10.1016/j.neuroimage.2013.10.027>
%%


params.RCGs = [2 4 8];
params.Dims = [128 256 512];
params.nb_facts = [2,4,8];

params.constraint = 'sp';
%params.constraint = 'sp_col';
%params.constraint = 'sp_row';
params.matrix_or_vector='vector';
%matrix_or_vector='matrix';
params.Nb_mult=500;
[RCGs,Dims,nb_facts,constraint,Nb_mult,matrix_or_vector]=deal(params.RCGs,params.Dims,params.nb_facts,params.constraint,params.Nb_mult,params.matrix_or_vector);




NDims=length(Dims);
NRCGs=length(RCGs);
Nnb_facts=length(nb_facts);








list_faust=cell(NDims,NRCGs,Nnb_facts);
list_dense=cell(NDims,1);

%% loading all the different faust and dense matrix
for j=1:NDims
       dim=Dims(j);
       A=rand(dim); 
       list_dense{j}=A; 
       

       
       for k=1:NRCGs
           RCG=RCGs(k);
           
           for l=1:Nnb_facts
               nfact=nb_facts(l);
               fact=gen_artificial_faust(dim,RCG,nfact,constraint);
               faust_transform=matlab_faust(fact);
               list_faust{j,k,l}=faust_transform;
               taillefaust=size(faust_transform);
               if((taillefaust(1) ~= dim) + (taillefaust(2) ~= dim))
                   error('invalid faust');
               end
               

           end
       end
end



%% time comparison

t_faust=zeros(Nb_mult+1,NDims,NRCGs,Nnb_facts);
t_dense=zeros(Nb_mult+1,NDims);

h = waitbar(0,'runtime comparison (faust vs dense matrix) for various configuration ...');
for i=1:Nb_mult+1
    waitbar(i/(Nb_mult+1));
    
   for j=1:NDims
       dim=Dims(j);
       
       
       
       for k=1:NRCGs
           RCG=RCGs(k);
           
           for l=1:Nnb_facts
               nfact=nb_facts(l);
               fact=gen_artificial_faust(dim,RCG,nfact,constraint);
               faust_transform=matlab_faust(fact);
               if strcmp(matrix_or_vector,'matrix')
                  x=rand(dim);
                  y=zeros(dim);
                  yfaust=zeros(dim); 
                   
               elseif strcmp(matrix_or_vector,'vector')
                  x=rand(dim,1);
                  y=zeros(dim,1);
                  yfaust=zeros(dim,1); 
                   
                   
               else
                   error('matrix_or_vector string must be equal to matrix or vector');
               end
               
               taillefaust=size(faust_transform);
               if((taillefaust(1) ~= dim)+(taillefaust(2) ~= dim))
                   error('invalid faust');
               end
               
               if(k==1)&&(l==1)
                   A=list_dense{j};
                   tic;
                   y=A*x;
                   t=toc;
                   t_dense(i,j)=t;
               end
               
               faust_transform=list_faust{j,k,l};
               tic;
               yfaust=faust_transform*x;
               t=toc;
               t_faust(i,j,k,l)=t; 
           end
       end
   end
end
close(h);

t_faust(1,:,:,:)=[];
t_dense(1,:)=[];

runPath=which(mfilename);
pathname = fileparts(runPath);
matfile = fullfile(pathname, 'runtime_comparison.mat');
save(matfile,'t_faust','t_dense','params');

