template<typename FPP,Device DEVICE,typename FPP2>
const char * Faust::HierarchicalFactFGFT<FPP,DEVICE,FPP2>::m_className="Faust::HierarchicalFactFGFT";


template<typename FPP,Device DEVICE,typename FPP2>
void Faust::HierarchicalFactFGFT<FPP,DEVICE,FPP2>::get_D(FPP* out_diag) const
{
	return dynamic_cast<Palm4MSAFGFT<FPP,Cpu,FPP2>*>(this->palm_global)->get_D(out_diag);
}

template<typename FPP,Device DEVICE,typename FPP2>
 const Faust::MatDense<FPP, DEVICE>& Faust::HierarchicalFactFGFT<FPP,DEVICE,FPP2>::get_D()
{
	return dynamic_cast<Palm4MSAFGFT<FPP,Cpu,FPP2>*>(this->palm_global)->get_D();
}
