template<typename FPP,Device DEVICE,typename FPP2>
const char * Faust::HierarchicalFactFFT<FPP,DEVICE,FPP2>::m_className="Faust::HierarchicalFactFFT";


template<typename FPP,Device DEVICE,typename FPP2>
void Faust::HierarchicalFactFFT<FPP,DEVICE,FPP2>::get_D(FPP* out_diag) const
{
	return dynamic_cast<Palm4MSAFFT<FPP,Cpu,FPP2>*>(this->palm_global)->get_D(out_diag);
}

template<typename FPP,Device DEVICE,typename FPP2>
 const Faust::MatDense<FPP, DEVICE>& Faust::HierarchicalFactFFT<FPP,DEVICE,FPP2>::get_D()
{
	return dynamic_cast<Palm4MSAFFT<FPP,Cpu,FPP2>*>(this->palm_global)->get_D();
}
