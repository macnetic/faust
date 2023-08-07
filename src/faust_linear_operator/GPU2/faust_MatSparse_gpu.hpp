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
			this->tocpu(cpu_smat);
			MatSparse<FPP2, Cpu> ffp2_cpu_smat;
			ffp2_cpu_smat.resize(cpu_smat.getNonZeros(), cpu_smat.getNbRow(), cpu_smat.getNbCol());
			//ffp2_cpu_smat = cpu_smat.template to_real<Real<FPP2>>();
			for(int i=0;i<cpu_smat.getNonZeros();i++)
			{
				ffp2_cpu_smat.getValuePtr()[i] = static_cast<FPP2>(std::real(std::complex<Real<FPP>>(cpu_smat.getValuePtr()[i])));
				ffp2_cpu_smat.getColInd()[i] = cpu_smat.getColInd()[i]; //TODO: memcpy
			}
			memcpy(ffp2_cpu_smat.getRowPtr(), cpu_smat.getRowPtr(), sizeof(int)*cpu_smat.getNbRow()+1);
			MatSparse<Real<FPP2>, GPU2> ffp2_gpu_smat(ffp2_cpu_smat);
			return ffp2_gpu_smat;
		}

	template<typename FPP>
		template<typename FPP2>
		MatSparse<FPP2, GPU2> MatSparse<FPP,GPU2>::cast() const
		{
			//TODO: should be done directly in GPU memory in gpu_mod (optimization)
			MatSparse<FPP, Cpu> cpu_smat;
			this->tocpu(cpu_smat);
			MatSparse<FPP2, Cpu> ffp2_cpu_smat;
			ffp2_cpu_smat.resize(cpu_smat.getNonZeros(), cpu_smat.getNbRow(), cpu_smat.getNbCol());
			for(int i=0;i<cpu_smat.getNonZeros();i++)
			{
				ffp2_cpu_smat.getValuePtr()[i] = static_cast<FPP2>(cpu_smat.getValuePtr()[i]);
				ffp2_cpu_smat.getColInd()[i] = cpu_smat.getColInd()[i]; //TODO: memcpy
			}
			memcpy(ffp2_cpu_smat.getRowPtr(), cpu_smat.getRowPtr(), sizeof(int)*cpu_smat.getNbRow()+1);
			MatSparse<FPP2, GPU2> ffp2_gpu_smat(ffp2_cpu_smat);
			return ffp2_gpu_smat;
		}

	template<typename FPP>
		MatDense<FPP, GPU2> MatSparse<FPP,GPU2>::to_dense() const
		{
			return MatDense<FPP, GPU2>(*this);
		}


	template<typename FPP>
		MatSparse<FPP, Cpu> MatSparse<FPP,GPU2>::tocpu() const
		{
			MatSparse<FPP, Cpu> cpu_smat(this->getNbRow(), this->getNbCol());
			tocpu(cpu_smat);
			return cpu_smat;
		}
}

