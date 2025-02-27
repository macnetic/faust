#include <cstring>

	template <typename FPP, FDevice DEVICE, typename FPP2>
Palm4MSAFGFT<FPP,DEVICE,FPP2>::Palm4MSAFGFT(const ParamsPalmFGFT<FPP, DEVICE, FPP2>& params, const bool isGlobal) : Palm4MSA<FPP,DEVICE,FPP2>(params, isGlobal), D(MatSparse<FPP,Cpu>(params.init_D))
{
}

template<typename FPP,FDevice DEVICE,typename FPP2>
Palm4MSAFGFT<FPP,DEVICE,FPP2>::Palm4MSAFGFT(const MatDense<FPP,DEVICE>& Lap, const ParamsFGFT<FPP,DEVICE,FPP2> & params, const bool isGlobal) : Palm4MSA<FPP,DEVICE,FPP2>(Lap, params, isGlobal), D(MatSparse<FPP,Cpu>(params.init_D))
{
}



template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::compute_grad_over_c()
{

	if(!this->isCComputed)
	{
		throw std::logic_error("this->c must be set before computing grad/this->c");
	}

	/*! \brief There are 4 ways to compute gradient : <br>
	 * (0) : lambda*(L'*(lambda*(L*this->S)*R - X))*R' : complexity = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2; <br>
	 * (1) : lambda*L'*((lambda*(L*this->S)*R - X)*R') : complexity = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2; <br>
	 * (2) : lambda*(L'*(lambda*L*(this->S*R) - X))*R' : complexity = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2; <br>
	 * (3) : lambda*L'*((lambda*L*(this->S*R) - X)*R') : complexity = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2; <br>
	 *  with L of size L1xL2 <br>
	 *       this->S of size L2xS2 <br>
	 *       R of size S2xR2 <br>
	 */
	unsigned long long int L1, L2, R2, S2;
	if (!this->isUpdateWayR2L)
	{
		L1 = (unsigned long long int) this->LorR.getNbRow();
		L2 = (unsigned long long int) this->LorR.getNbCol();
		R2 = (unsigned long long int) this->RorL[this->m_indFact].getNbCol();
	}
	else
	{
		L1 = (unsigned long long int) this->RorL[this->m_indFact].getNbRow();
		L2 = (unsigned long long int) this->RorL[this->m_indFact].getNbCol();
		R2 = (unsigned long long int) this->LorR.getNbCol();
	}
	S2 = (unsigned long long int) this->S[this->m_indFact].getNbCol();
	vector<unsigned long long int > complexity(4,0);
	complexity[0] = L1*L2*S2 + L1*S2*R2 + L2*L1*R2 + L2*R2*R2;
	complexity[1] = L1*L2*S2 + L1*S2*R2 + L1*R2*S2 + L2*L1*S2;
	complexity[2] = L2*S2*R2 + L1*L2*R2 + L2*L1*R2 + L2*R2*S2;
	complexity[3] = L2*S2*R2 + L1*L2*R2 + L1*R2*S2 + L2*L1*S2;

	int idx = distance(complexity.begin(), min_element(complexity.begin(), complexity.end()));
	this->error = this->data;
	Faust::MatDense<FPP,DEVICE> tmp1,tmp2,tmp3;

	if (idx==0 || idx==1) // computing L*this->S first, then (L*this->S)*R and finally the this->error
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp1 = L*this->S
			multiply(this->LorR, this->S[this->m_indFact], tmp1);
			// tmp2 = L*this->S*R
			multiply(tmp1, this->RorL[this->m_indFact], tmp2);
		}
		else
		{
			// tmp1 = L*this->S
			multiply(this->RorL[this->m_indFact], this->S[this->m_indFact], tmp1);
			// tmp2 = L*this->S*R
			multiply(tmp1, this->LorR, tmp2);
		}
	}
	else // computing this->S*R first, then L*(this->S*R)
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp1 = this->S*R
			multiply(this->S[this->m_indFact], this->RorL[this->m_indFact], tmp1);
			// tmp2 = L*this->S*R
			multiply(this->LorR, tmp1, tmp2);
		}
		else
		{
			// tmp1 = this->S*R
			multiply(this->S[this->m_indFact], this->LorR, tmp1);
			// tmp2 = L*this->S*R
			multiply(this->RorL[this->m_indFact], tmp1, tmp2);
		}
	}
	// tmp1 = L*this->S*R*D
	tmp1 = tmp2;
	tmp1 *= D;
	// this->error = lambda*tmp1*lambda*tmp2'-data // this->error is data before this call
	gemm(tmp1, tmp2, this->error, FPP(this->m_lambda*this->m_lambda), (FPP)-1.0, 'N', this->TorH);

	if (idx==0 || idx==2) // computing L'*this->error first, then (L'*this->error)*R'
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp3 = this->m_lambda*L'*this->error (= this->m_lambda*L' * (this->m_lambda*L*this->S*R - data) )
			gemm(this->LorR, this->error, tmp3, FPP(this->m_lambda),(FPP) 0.0, this->TorH, 'N');
			// tmp2 = lambda*L*this->S*R*D*R'
			gemm(tmp1, this->RorL[this->m_indFact], tmp2, FPP(this->m_lambda), (FPP) 0, 'N', this->TorH);
		}
		else
		{
			// tmp3 = this->m_lambda*L'*this->error (= this->m_lambda*L' * (this->m_lambda*L*this->S*R - data) )
			gemm(this->RorL[this->m_indFact], this->error, tmp3, FPP(this->m_lambda), (FPP) 0.0, this->TorH, 'N');
			// tmp2 = lambda*L*this->S*R*D*R'
			gemm(tmp1, this->LorR, tmp2, FPP(this->m_lambda), (FPP) 0, 'N', this->TorH);
		}
		// grad_over_c = 1/this->c*tmp3*tmp2
		gemm(tmp3, tmp2, this->grad_over_c, (FPP) (1.0/this->c), (FPP) 0.0,'N','N');

	}
	else // computing this->error*R' first, then L'*(this->error*lambda*LSRD*R')
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp2 = lambda*tmp1*R' = lambda*LSRD*R'
			gemm(tmp1, this->RorL[this->m_indFact], tmp2, (FPP) this->m_lambda, (FPP) 0, 'N', this->TorH);
			// tmp3 = this->m_lambda*this->error*tmp2
			gemm(this->error, tmp2, tmp3, FPP(this->m_lambda), (FPP) 0.0, 'N', 'N');
			// grad_over_c = 1/this->c*L'*tmp3
			gemm(this->LorR, tmp3, this->grad_over_c,(FPP) (1.0/this->c), (FPP) 0.0,this->TorH,'N');
		}
		else
		{
			// tmp2 = lambda*tmp1*R' = lambda*LSRD*R'
			gemm(tmp1, this->LorR, tmp2, (FPP) this->m_lambda, (FPP) 0, 'N', this->TorH);
			// tmp3 = this->m_lambda*this->error*tmp2
			gemm(this->error, tmp2, tmp3, FPP(this->m_lambda), (FPP) 0.0, 'N', 'N');
			// grad_over_c = 1/this->c*L'*tmp3
			gemm(this->RorL[this->m_indFact], tmp3, this->grad_over_c, (FPP) (1.0/this->c), (FPP) 0.0,this->TorH,'N');
		}

	}


	this->isGradComputed = true;

#ifdef __COMPILE_TIMERS__
	t_global_compute_grad_over_c.stop();
	t_local_compute_grad_over_c.stop();
#endif
}


template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::compute_lambda()
{
	// override parent's method
	// Xhat = (S[0]*...*S[nfact-1])*D*(S[0]*...*S[nfact-1])'
	// Xhat = LorR*D*LorR' //  LorR equals the prod of all factors after their update iterations (in loop of next_step())
	MatDense<FPP,Cpu> tmp;
	// tmp = D*LorR'
	Faust::spgemm(this->D, this->LorR, tmp, (FPP) 1.0, (FPP) 0.0, 'N', this->TorH);
	// LorR = LorR*tmp
	gemm(this->LorR, tmp, D_grad_over_c, (FPP) 1.0, (FPP) 0.0, 'N', 'N');
	//NOTE: D_grad_over_c has nothing to do here but is equal to LorR*D*LorR'
	//		this product is thus not re-computed in compute_D_grad_over_c()
	// at this stage we can rely on parent function to compute lambda
	Palm4MSA<FPP,DEVICE,FPP2>::compute_lambda(D_grad_over_c);
	// reset LorR at the factor product to continue next iterations
	//then we finish the lambda computation with a sqrt() (Fro. norm)
	this->m_lambda = std::sqrt(/*Faust::abs(*/this->m_lambda);//);
	// (that's an additional operation in Palm4MSAFGFT)
}

template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::next_step()
{
	Palm4MSA<FPP, Cpu, FPP2>::next_step();
	// besides to what the parent has done
	// we need to update D
	this->compute_D();
}


template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::compute_D()
{
	// besides to what the parent has done
	// we need to update D
	compute_D_grad_over_c();
//#pragma omp parallel for schedule(static)
	for(faust_unsigned_int i = 0; i < D.getNbCol();i++)
		D.getValuePtr()[i] = D.getValuePtr()[i]-D_grad_over_c.getData()[i*D.getNbCol()+i];
}

template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::compute_D_grad_over_c()
{
	// Uhat = lambda*LorR
	// grad = 0.5*Uhat'*(Uhat*D*Uhat' - X)*Uhat
	MatDense<FPP, Cpu> tmp;
	//compute_lambda has already computed D_grad_over_c = LorR*D*LorR'
	D_grad_over_c.scalarMultiply(FPP(this->m_lambda*this->m_lambda));
	D_grad_over_c -= this->data;
	// tmp = Uhat'*(Uhat*D*Uhat' - X)
	gemm(this->LorR, D_grad_over_c, tmp, (FPP) this->m_lambda, (FPP) 0., this->TorH, 'N');
	// D_grad_over_c = 0.5*Uhat'*(Uhat*D*Uhat' - X)*Uhat
	gemm(tmp, this->LorR, D_grad_over_c, (FPP) FPP(.5*this->m_lambda/this->c), (FPP) 0., 'N', 'N');
}

template <typename FPP, FDevice DEVICE, typename FPP2>
const MatSparse<FPP, DEVICE>& Palm4MSAFGFT<FPP,DEVICE,FPP2>::get_D()
{
	return this->D;
}

template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::get_D(FPP* diag_data)
{
	memcpy(diag_data, this->D.getValuePtr(), sizeof(FPP)*this->D.getNbCol());
}

template <typename FPP, FDevice DEVICE, typename FPP2>
void Palm4MSAFGFT<FPP,DEVICE,FPP2>::compute_c()
{
	//do nothing because the Palm4MSAFGFT has always a constant step size
	this->isCComputed = true;
}

