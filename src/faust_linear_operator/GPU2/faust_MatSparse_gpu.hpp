namespace Faust {
	template<typename FPP>
		MatGeneric<FPP,GPU2>* MatSparse<FPP,GPU2>::Clone(const bool isOptimize /*default value=false*/) const
		{
			if (isOptimize)
			{
				MatDense<FPP,GPU2> M((*this));
				return optimize(M,(*this));
			} else
			{
				return new MatSparse<FPP,GPU2>((*this));
			}
		}

	template<typename FPP>
		template<typename FPP2>
		MatSparse<Real<FPP2>, GPU2> MatSparse<FPP,GPU2>::to_real() const
		{
			//TODO: should be done directly in GPU memory in gpu_mod (optimization)
			MatSparse<FPP, Cpu> cpu_smat;
			auto ffp2_cpu_smat = cpu_smat.template to_real<Real<FPP2>>();
			MatSparse<Real<FPP2>, GPU2> ffp2_gpu_smat(ffp2_cpu_smat);
			return ffp2_gpu_smat;
		}


	template<typename FPP>
		MatDense<FPP, GPU2> MatSparse<FPP,GPU2>::to_dense() const
		{
			return MatDense<FPP, GPU2>(*this);
		}
}

