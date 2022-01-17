template<typename FPP>
template<typename FPP2>
Faust::MatDense<Real<FPP2>, GPU2> Faust::MatDense<FPP,GPU2>::to_real() const
{
	//TODO: should be done directly in GPU memory in gpu_mod (optimization)
	MatDense<FPP, Cpu> cpu_dmat;
	auto ffp2_cpu_dmat = cpu_dmat.template to_real<Real<FPP2>>();
	MatDense<Real<FPP2>, GPU2> ffp2_gpu_dmat(ffp2_cpu_dmat);
	return ffp2_gpu_dmat;
}
