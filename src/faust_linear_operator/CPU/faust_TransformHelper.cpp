#include <faust_TransformHelper.h>

namespace Faust
{

	template<>
		template<>
		TransformHelper<double, Cpu>* TransformHelper<std::complex<double>, Cpu>::real<double>()
		{
			std::vector<MatGeneric<double,Cpu>*> real_data;
			MatSparse<std::complex<double>, Cpu> *curfac_sp;
			MatDense<std::complex<double>, Cpu> *curfac_ds;
			for(auto curfac: this->transform->data)
			{
				if(curfac_ds = dynamic_cast<MatDense<std::complex<double>, Cpu>*>(curfac))
				{
					auto real_fac = new MatDense<double,Cpu>(curfac->getNbRow(), curfac->getNbCol());
					curfac_ds->real(*real_fac);
					real_data.push_back(real_fac);
				}
				else if(curfac_sp = dynamic_cast<MatSparse<std::complex<double>, Cpu>*>(curfac))
				{
					auto real_fac = new MatSparse<double,Cpu>(curfac->getNbRow(), curfac->getNbCol());
					curfac_sp->real(*real_fac);
					real_data.push_back(real_fac);
				}
				else
				{
					throw std::runtime_error("real() failed because a factor is neither a MatDense or a MatSparse");
				}
			}
			return new TransformHelper<double, Cpu>(real_data, 1.0, false, false, true);
		}

	template<>
		template<>
		TransformHelper<double, Cpu>* TransformHelper<double, Cpu>::real<double>()
		{
			return this->clone();
		}
}
