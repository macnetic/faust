function [facts,D,err,L,coord_choices] = diagonalization_givens_parall(Lap,J,t)
%% Description diagonalization_givens_parall 
%  Greedy matrix diagonalization with prallel Givens rotations.
%  [facts,D,err,L,coord_choices] = diagonalization_givens_parall(Lap,J,t)
%  runs the matrix diagonalization algorithm (Algorithm in figure 4 of [1]) 
%  on the specified input matrix "Lap" with "J" Givens rotations with "t" 
%  Givens rotations chosen in parallel at each step. It returns the 
%  obtained Givens rotations in the cell array "facts", the approximate 
%  eigenvalues matrix in "D", the error in "err", the corresponding
%  approximately diagonalized matrix in "L" and the coordinates of the
%  chosen Givens rotations in "coord_choices".
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
%	Luc Le Magoarou: luc.le-magoarou@inria.fr
%	Remi Gribonval : remi.gribonval@inria.fr
%
%% References:
% [1]	Le Magoarou L., Gribonval R. and Tremblay N., "Approximate fast 
%   graph Fourier transforms via multi-layer sparse approximations", 
%   submitted to IEEE Transactions on Signal and Information Processing
%   over Networks.
%	<https://hal.inria.fr/hal-01416110>
%%


facts = cell(1,round(J/t));
n=size(Lap,1);
L=sparse(Lap);
C = zeros(n);
err=zeros(1,round(J/t));
%coord_choices = zeros(2,J);
coord_choices = []
%N_edges = (nnz(Lap)-n)/2;




for j=1:round(J/t)
    
    L_low = abs(tril(L,-1));
    
    %%%%%%%%% Un tri
    ind_nnz = find(L_low);
    %disp(['nnz in L :' num2str(numel(ind_nnz))])
	[~,ind_sorted] = sort(L_low(ind_nnz),'descend');
    [Rvec, Svec] = ind2sub([n,n],ind_nnz(ind_sorted));
    RSmat = [Rvec, Svec];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% t max
    %     RSmat=[];
    %     Nrot=0;
    %     %ind_nnz = find(L_low);
    %     while(Nrot<t && any(L_low(:)))
    %         %ind_nnz_old = ind_nnz;
    %         %ind_nnz = ind_nnz_start(L_low(ind_nnz_start)~=0);
    %         ind_nnz = find(L_low);
    %         %ind_nnz = ind_nnz_old(ind_nnz);
    %         [~,ind_max] = max(L_low(ind_nnz));
    %         [r, s] = ind2sub([n,n],ind_nnz(ind_max));
    %         RSmat = [RSmat ; r, s];
    %         L_low(r,:)=0; L_low(:,r)=0;
    %         L_low(s,:)=0; L_low(:,s)=0;
    %
    %         %          [rows,cols] = ind2sub([n,n],ind_nnz);
    %         %          nnz_sub = [rows,cols];
    %         %          other = xor(nnz_sub==r,nnz_sub==s);
    %         %          nnz_sub( any(other,2), : ) = [];
    %         %          ind_nnz = sub2ind([n,n],nnz_sub(:,1),nnz_sub(:,2));
    %
    %         %         ind_nnz(ind_nnz>= (r-1)*n+1 & ind_nnz<= r*n)=[]; % suppress rth column
    %         %         ind_nnz(ind_nnz>= (s-1)*n+1 & ind_nnz<= s*n)=[]; % suppress sth column
    %         %         ind_nnz(mod(ind_nnz,n)==r)=[]; % suppress rth row
    %         %         ind_nnz(mod(ind_nnz,n)==s)=[]; % suppress sth row
    %         Nrot=Nrot+1;
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    S = eye(n);
    Nrot=0;
    k=1;
    % chosen = [];
    while(Nrot<t && k<= size(RSmat,1))
        %for k=1:size(RSmat,1)
        %if RSmat(k,1) &&  RSmat(k,2)
            p = RSmat(1,1);
            q = RSmat(1,2);
            other = RSmat==p | RSmat==q;
            tokill = any(other,2);
            %RSmat( tokill, : ) = [];
            RSmatnew = RSmat(~tokill,:);
            RSmat = RSmatnew;
            % if  ~any(abs(p-chosen)<0.5) && ~any(abs(q-chosen)<0.5)%numel(find(S(p,:)))==1 && numel(find(S(:,q)))==1
            %coord_choices(1,j) = p;
            %coord_choices(2,j) = q;
			coord_choices = [ coord_choices, [p;q]];
            theta = atan2(L(q,q) - L(p,p),2*L(p,q))/2 + pi/4;
            S(p,p) = cos(theta); S(p,q) = -sin(theta);
            S(q,p) = sin(theta); S(q,q) = cos(theta);
            Nrot = Nrot+1;
            k=k+1;
            %             chosen = [chosen, p, q];
%         else
%             k=k+1;
%         end
    end
    L = S'*L*S;
    facts{j} = sparse(S);
	D = spdiag(diag(L)); 
    % D = sparse(diag(diag(L)));
    %    errchg = err - norm(D-L,'fro')^2
    if mod(j,1)==0
        %err(j) = norm(D-L,'fro')^2/norm(L,'fro')^2;
        err(j) = norm(D-L,'fro')^2/norm(Lap,'fro')^2;
        %err(j) = norm(D-L)/norm(Lap);
        disp(['Iter ' num2str(j) ', error = ' num2str(err(j))])
        %    disp(['Number of edges: ' num2str(N_edges)])
    end
end
%----------------- Sort the spectrum
spectrum = diag(D);
[~,Ispec] = sort(spectrum);
D = D(Ispec,Ispec);
facts{round(J/t)} = facts{round(J/t)}(:,Ispec);
end
