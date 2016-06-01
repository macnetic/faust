#ifndef __FAUST_GPU_CORE_2_CPU_CORE_HPP__
#define __FAUST_GPU_CORE_2_CPU_CORE_HPP__


template<typename FPP, typename FPP1>
void faust_gpu2cpu(Faust::Transform<FPP,Cpu>& fcore, const Faust::Transform<FPP1,Gpu>& cu_fcore)
{

	std::vector< Faust::MatSparse<FPP,Cpu> > list_spmat(cu_fcore.size());
	std::vector< Faust::MatSparse<FPP1,Gpu> > list_cu_spmat;
	cu_fcore.get_facts(list_cu_spmat);

	for (int i=0;i<cu_fcore.size();i++)
	{
		faust_gpu2cpu(list_spmat[i],list_cu_spmat[i]);
	}
	Faust::Transform<FPP,Cpu> fcorebis(list_spmat);
	fcore=fcorebis;
}

#endif

