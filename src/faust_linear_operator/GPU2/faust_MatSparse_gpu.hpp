template<typename FPP>
template<typename FPP2>
Faust::MatSparse<Real<FPP2>, GPU2> Faust::MatSparse<FPP,GPU2>::to_real() const
{
	//TODO: should be done directly in GPU memory in gpu_mod (optimization)
	MatSparse<FPP, Cpu> cpu_smat;
	auto ffp2_cpu_smat = cpu_smat.template to_real<Real<FPP2>>();
	MatSparse<Real<FPP2>, GPU2> ffp2_gpu_smat(ffp2_cpu_smat);
	return ffp2_gpu_smat;
}

