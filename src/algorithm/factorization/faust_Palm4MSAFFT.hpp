#include <cstring>

	template <typename FPP, Device DEVICE, typename FPP2>
Palm4MSAFFT<FPP,DEVICE,FPP2>::Palm4MSAFFT(const ParamsPalmFFT<FPP, DEVICE, FPP2>& params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal) : Palm4MSA<FPP,DEVICE,FPP2>(params, blasHandle, isGlobal), D(params.init_D)
{
	//TODO: is there something to check additionally to what parent's ctor checks ?
}

template<typename FPP,Device DEVICE,typename FPP2>
Palm4MSAFFT<FPP,DEVICE,FPP2>::Palm4MSAFFT(const MatDense<FPP,DEVICE>& Lap, const ParamsFFT<FPP,DEVICE,FPP2> & params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal) : Palm4MSA<FPP,DEVICE,FPP2>(Lap, params, blasHandle, isGlobal), D(params.init_D)
{

}

template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_grad_over_c()
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
			multiply(this->LorR, this->S[this->m_indFact], tmp1, this->blas_handle);
			// tmp2 = L*this->S*R
			multiply(tmp1, this->RorL[this->m_indFact], tmp2, this->blas_handle);
		}
		else
		{
			// tmp1 = L*this->S
			multiply(this->RorL[this->m_indFact], this->S[this->m_indFact], tmp1, this->blas_handle);
			// tmp2 = L*this->S*R
			multiply(tmp1, this->LorR, tmp2, this->blas_handle);
		}
	}
	else // computing this->S*R first, then L*(this->S*R)
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp1 = this->S*R
			multiply(this->S[this->m_indFact], this->RorL[this->m_indFact], tmp1, this->blas_handle);
			// tmp2 = L*this->S*R
			multiply(this->LorR, tmp1, tmp2, this->blas_handle);
		}
		else
		{
			// tmp1 = this->S*R
			multiply(this->S[this->m_indFact], this->LorR, tmp1, this->blas_handle);
			// tmp2 = L*this->S*R
			multiply(this->RorL[this->m_indFact], tmp1, tmp2, this->blas_handle);
		}
	}
	// tmp1 = L*this->S*R*D //TODO: review the mul with D being MatSparse
	//TODO: optimize by determining the best product order regarding computation time
	multiply(tmp2, D, tmp1, this->blas_handle);
	// this->error = lambda*tmp1*lambda*tmp2'-data // this->error is data before this call
	gemm(tmp1, tmp2, this->error, this->m_lambda*this->m_lambda, (FPP)-1.0, 'N', 'T', this->blas_handle);

	//false is for disabling evaluation (because the transpose does it later)
	this->LorR.conjugate(false);
	this->RorL[this->m_indFact].conjugate(false);

	if (idx==0 || idx==2) // computing L'*this->error first, then (L'*this->error)*R'
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp3 = this->m_lambda*L'*this->error (= this->m_lambda*L' * (this->m_lambda*L*this->S*R - data) )
			gemm(this->LorR, this->error, tmp3, this->m_lambda,(FPP) 0.0, 'T', 'N', this->blas_handle);
			// tmp2 = lambda*L*this->S*R*D*R'
			gemm(tmp1, this->RorL[this->m_indFact], tmp2, this->m_lambda, (FPP) 0, 'N', 'T', this->blas_handle);
		}
		else
		{
			// tmp3 = this->m_lambda*L'*this->error (= this->m_lambda*L' * (this->m_lambda*L*this->S*R - data) )
			gemm(this->RorL[this->m_indFact], this->error, tmp3, this->m_lambda, (FPP) 0.0, 'T', 'N', this->blas_handle);
			// tmp2 = lambda*L*this->S*R*D*R'
			gemm(tmp1, this->LorR, tmp2, this->m_lambda, (FPP) 0, 'N', 'T', this->blas_handle);
		}
		// grad_over_c = 1/this->c*tmp3*tmp2
		gemm(tmp3, tmp2, this->grad_over_c, (FPP) 1.0/this->c, (FPP) (FPP) 0.0,'N','N', this->blas_handle);

	}
	else // computing this->error*R' first, then L'*(this->error*lambda*LSRD*R')
	{
		if (!this->isUpdateWayR2L)
		{
			// tmp2 = lambda*tmp1*R' = lambda*LSRD*R'
			gemm(tmp1, this->RorL[this->m_indFact], tmp2, (FPP) this->m_lambda, (FPP) 0, 'N', 'T', this->blas_handle);
			// tmp3 = this->m_lambda*this->error*tmp2
			gemm(this->error, tmp2, tmp3, this->m_lambda, (FPP) 0.0, 'N', 'N', this->blas_handle);
			// grad_over_c = 1/this->c*L'*tmp3
			gemm(this->LorR, tmp3, this->grad_over_c,(FPP) 1.0/this->c, (FPP) 0.0,'T','N', this->blas_handle);
		}
		else
		{
			// tmp2 = lambda*tmp1*R' = lambda*LSRD*R'
			gemm(tmp1, this->LorR, tmp2, (FPP) this->m_lambda, (FPP) 0, 'N', 'T', this->blas_handle);
			// tmp3 = this->m_lambda*this->error*tmp2
			gemm(this->error, tmp2, tmp3, this->m_lambda, (FPP) 0.0, 'N', 'N', this->blas_handle);
			// grad_over_c = 1/this->c*L'*tmp3
			gemm(this->RorL[this->m_indFact], tmp3, this->grad_over_c, (FPP) 1.0/this->c, (FPP) 0.0,'T','N', this->blas_handle);
		}

	}

	//TODO: avoid type checking by adding another template function (for complex and real types) or function pointer
	if(typeid(this->data.getData()[0]) == typeid(complex<float>) || typeid(this->data.getData()[0]) == typeid(complex<double>))
	{
		this->LorR.conjugate();
		this->RorL[this->m_indFact].conjugate();
	}
	this->isGradComputed = true;

#ifdef __COMPILE_TIMERS__
	t_global_compute_grad_over_c.stop();
	t_local_compute_grad_over_c.stop();
#endif
}


template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_lambda()
{
	//TODO: override parent's method
	// Xhat = (S[0]*...*S[nfact-1])*D*(S[0]*...*S[nfact-1])'
	// Xhat = LorR*D*LorR' //  LorR equals the prod of all factors after their update iterations (in loop of next_step())
	MatDense<FPP,Cpu> tmp;
	// tmp = D*LorR'
	gemm(this->D, this->LorR, tmp, (FPP) 1.0, (FPP) 0.0, 'N', 'T', this->blas_handle);
	// LorR = LorR*tmp
	gemm(this->LorR, tmp, D_grad_over_c, (FPP) 1.0, (FPP) 0.0, 'N', 'N', this->blas_handle);
	tmp = this->LorR;
	this->LorR = D_grad_over_c;
	//NOTE: D_grad_over_c has nothing to do here but is equal to LorR*D*LorR*D'
	//		this product is thus not re-computed in compute_D_grad_over_c()
	//TODO: avoid all these copies
	// at this stage we can rely on parent function to compute lambda
	Palm4MSA<FPP,DEVICE,FPP2>::compute_lambda();
	// reset LorR at the factor product to continue next iterations
	this->LorR = tmp;
	//then we finish the lambda computation with a sqrt() (Fro. norm)
	this->m_lambda = std::sqrt(this->m_lambda);
	// (that's an additional operation in Palm4MSAFFT)
}

template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::next_step()
{
	Palm4MSA<FPP, Cpu, FPP2>::next_step();
	// besides to what the parent has done
	// we need to update D
	this->compute_D();
}

template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_D()
{
	// besides to what the parent has done
	// we need to update D
	compute_D_grad_over_c();
	D_grad_over_c.scalarMultiply(this->m_lambda/this->c);
	D -= D_grad_over_c;
	//TODO: optimize MatSparse + no-copy (Eigen::DiagonalMatrix ?)
	FPP * data = new FPP[D.getNbRow()*D.getNbCol()];
	memset(data, 0, sizeof(FPP)*D.getNbRow()*D.getNbCol());
	for(faust_unsigned_int i = 0; i < D.getNbCol();i++)
		data[i*D.getNbCol()+i] = D[i*D.getNbCol()+i];
	D = MatDense<FPP,Cpu>(data, D.getNbRow(), D.getNbCol());
}

template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_D_grad_over_c()
{
	// grad = 0.5*LorR'*(LorR*D*LorR' - X)*LorR
	MatDense<FPP, Cpu> tmp;
	//compute_lambda has already compute D_grad_over_c = LorR*D*LorR'
	D_grad_over_c -= this->data;
	//TODO: opt. by determining best order of product
	// tmp = LorR'*(LorR*D*LorR' - X)
	gemm(this->LorR, D_grad_over_c, tmp, (FPP) 1., (FPP) 0., 'T', 'N', this->blas_handle);
	// D_grad_over_c = LorR'*(LorR*D*LorR' - X)*LorR
	gemm(tmp, this->LorR, D_grad_over_c, (FPP) 1., (FPP) 0., 'N', 'N', this->blas_handle);
}

template <typename FPP, Device DEVICE, typename FPP2>
const MatDense<FPP, DEVICE>& Palm4MSAFFT<FPP,DEVICE,FPP2>::get_D()
{
	return this->D;
}

