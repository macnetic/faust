using namespace Faust;

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::next_step()
{

}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::choose_pivot()
{
//Matlab ref. code:
//        [~,p] = min(C_min_row);
//        q = q_candidates(p);
//        coord_choices(1,j) = p;
//        coord_choices(2,j) = q;
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::sort_L_in_C()
{
// Matlab ref. code:
//**************  at initialization
//    for r=1:n
//        for s=r+1:n
//            C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
//        end
//    end
//    [C_min_row, q_candidates] = min(C,[],2);
//
/************ at end of ite.
*       for r=[p,q]
*           for s=r+1:n
*               C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
*           end
*           [C_min_row(r), q_candidates(r)] = min(C(r,:));
*       end
*
* 		for s=[p,q]
*            for r=1:s-1
*                C(r,s) = -abs(L(r,s));%-2*L(r,s)^2;
*                if C(r,s)<C_min_row(r)
*                    C_min_row(r) = C(r,s);
*                    q_candidates(r) = s;
*                elseif q_candidates(r) == s
*                    [C_min_row(r), q_candidates(r)] = min(C(r,:));
*                end
*            end
*        end
*/
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::calc_theta()
{
// Matlab ref. code:
//		theta1 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 ;
//        err_theta1 = (L(p,q)*cos(2*theta1) + 0.5*(L(q,q) - L(p,p))*sin(2*theta1))^2;
//        theta2 = atan2(L(q,q) - L(p,p),2*L(p,q))/2 + pi/4;
//        err_theta2 = (L(p,q)*cos(2*theta2) + 0.5*(L(q,q) - L(p,p))*sin(2*theta2))^2;
//        if err_theta1<err_theta2
//            theta=theta1;
//        else
//            theta=theta2;
//        end
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_fact()
{
// Matlab ref. code:
//        S = eye(n);
//        S(p,p) = cos(theta); S(p,q) = -sin(theta);
//        S(q,p) = sin(theta); S(q,q) = cos(theta);
//        S = sparse(S);
//        facts{j} = S;
//
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_L()
{
// L = S'*L*S
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::update_D()
{
// D = spdiag(diag(L))
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::calc_err()
{
// Matlab ref. code:
//        if mod(j,100)==0
//            %err(j) = norm(D-L,'fro')^2/norm(L,'fro')^2;
//            err(j) = norm(D-L,'fro')^2/norm(Lap,'fro')^2;
//            %err(j) = norm(D-L)/norm(Lap);
//            disp(['Iter ' num2str(j) ', error = ' num2str(err(j))])
//            %    disp(['Number of edges: ' num2str(N_edges)])
//        end
//
}

template<typename FPP, Device DEVICE, typename FPP2>
void GivensFGFT<FPP,DEVICE,FPP2>::compute_facts()
{

}

template<typename FPP, Device DEVICE, typename FPP2>
GivensFGFT<FPP,DEVICE,FPP2>::GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, faust_unsigned_int J) : Lap(Lap), facts(J), D(Lap.getNbRow(), Lap.getNbCol()), C(Lap.getNbRow(), Lap.getNbRow()), errs(J), coord_choices(J)
{
	/** Matlab ref. code:
	 *     facts = cell(1,J);
	 *     n=size(Lap,1);
	 *     L=Lap;
	 *     C = 15*ones(n);
	 *     err=zeros(1,J);
	 *     coord_choices = zeros(2,J);
	 *
	 */
}



