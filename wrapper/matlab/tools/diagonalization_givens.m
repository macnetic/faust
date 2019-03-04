function [facts,D,err,L,coord_choices] = diagonalization_givens(Lap,J)
	%% Description diagonalization_givens 
	%  Greedy matrix diagonalization with Givens rotations.
	%  [facts,D,err,L,coord_choices] = diagonalization_givens(Lap,J) runs the 
	%  matrix diagonalization algorithm (Algorithm in figure 2 of [1]) on the 
	%  specified input matrix "Lap" with "J" Givens rotations. It returns the 
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


	%facts = zeros(J,3);
	facts = cell(1,J);
	n=size(Lap,1);
	L=Lap;
	C = 15*ones(n);
	err=zeros(1,J);
	coord_choices = zeros(2,J);
	%N_edges = (nnz(Lap)-n)/2;


	for r=1:n
		for s=r+1:n
			C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
		end
	end


	[C_min_row, q_candidates] = min(C,[],2);


	for j=1:J

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%     A_old = A;B_old=B;C_old=C;
		%     for r=1:n
		%     for s=r+1:n
		%         A(r,s) = 0.5*(L(s,s) - L(r,r))^2 - 2*L(r,s)^2;
		%         B(r,s) = L(r,s)*(L(s,s) - L(r,r));
		%         C(r,s) = A(r,s)/2 - sqrt(A(r,s)^2/4 + B(r,s)^2);
		%     end
		% end
		% [C_max_row, q_candidates] = min(C,[],2);

		%     figure;imagesc(A - A_old);
		%     figure;imagesc(B - B_old);
		%     figure;imagesc(C - C_old);
		%     pause

		%     [~,I] = nanmin(C(:));
		%     [p,q] = ind2sub(size(C),I);
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		[~,p] = min(C_min_row);
		q = q_candidates(p);



		coord_choices(1,j) = p;
		coord_choices(2,j) = q;
		theta1 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 ;
		err_theta1 = (L(p,q)*cos(2*theta1) + 0.5*(L(q,q) - L(p,p))*sin(2*theta1))^2;
		theta2 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 + pi/4;
		err_theta2 = (L(p,q)*cos(2*theta2) + 0.5*(L(q,q) - L(p,p))*sin(2*theta2))^2;
		if err_theta1<err_theta2
			theta=theta1;
		else
			theta=theta2;
		end




		S = eye(n);
		S(p,p) = cos(theta); S(p,q) = -sin(theta);
		S(q,p) = sin(theta); S(q,q) = cos(theta);
		S = sparse(S);
		facts{j} = S;
%		if(j == 105)
%			input B
%		end
		%   facts(j,:) = [p,q,theta];



		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%verif = A*sin(2*theta)^2 + 2*B*sin(2*theta)*cos(2*theta);
		%verif(p,q)
		%
		%     p
		%     q
		%     theta1
		%     err_theta1
		%     theta2
		%     err_theta2
		%     theta
		%         angle = -2*pi:0.001:2*pi;
		%         figure
		%         plot(angle,(L(p,q)*cos(2*angle) + 0.5*(L(q,q) - L(p,p))*sin(2*angle)).^2,'b');
		% %         hold on
		% %         plot(angle,-A(p,q)/2*cos(4*angle) + B(p,q)*sin(4*angle),'r--');
		% %         hold on
		% %         plot(angle,sqrt(A(p,q)^2/4 + B(p,q)^2)*cos(4*angle - atan2(-2*B(p,q),A(p,q)) + pi),'g:');
		%         pause
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		%     c=cos(theta);
		%     s=sin(theta);
		%     L_old = L;
		%     L(p,:) = c * L_old(p,:)  + s * L_old(q,:);
		%     L(q,:) = -s * L_old(p,:)  + c * L_old(q,:);
		%     L_old = L;
		%     L(:,p) = c * L_old(:,p)  + s * L_old(:,q);
		%     L(:,q) = -s * L_old(:,p)  + c * L_old(:,q);
		% N_edges = (nnz(L)-n)/2;
		L = S'*L*S;



		D = sparse(diag(diag(L)));
		%    errchg = err - norm(D-L,'fro')^2
		if mod(j,100)==0
			%err(j) = norm(D-L,'fro')^2/norm(L,'fro')^2;
			err(j) = norm(D-L,'fro')^2/norm(Lap,'fro')^2;
			%err(j) = norm(D-L)/norm(Lap);
			disp(['Iter ' num2str(j) ', error = ' num2str(err(j))])
			%    disp(['Number of edges: ' num2str(N_edges)])
		end

		for r=[p,q]
			for s=r+1:n
				C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
			end
			[C_min_row(r), q_candidates(r)] = min(C(r,:));
		end


		for s=[p,q]
			for r=1:s-1
				C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
				if C(r,s)<C_min_row(r)
					C_min_row(r) = C(r,s);
					q_candidates(r) = s;
				elseif q_candidates(r) == s
					[C_min_row(r), q_candidates(r)] = min(C(r,:));
				end
			end
		end



	end
end
