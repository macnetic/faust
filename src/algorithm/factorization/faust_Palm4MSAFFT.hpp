
	template <typename FPP, Device DEVICE, typename FPP2>
Palm4MSAFFT<FPP,DEVICE,FPP2>::Palm4MSAFFT(const ParamsPalmFFT<FPP, DEVICE, FPP2>& params, const BlasHandle<DEVICE> blasHandle, const bool isGlobal) : Palm4MSA<FPP,DEVICE,FPP2>(params, blasHandle, isGlobal)
{
	//TODO: manage init_D ?

}


template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_grad_over_c()
{
	//TODO: override parent's method
	Palm4MSA<FPP,DEVICE,FPP2>::compute_grad_over_c();
}


template <typename FPP, Device DEVICE, typename FPP2>
void Palm4MSAFFT<FPP,DEVICE,FPP2>::compute_lambda()
{
	//TODO: override parent's method
	Palm4MSA<FPP,DEVICE,FPP2>::compute_lambda();
}
