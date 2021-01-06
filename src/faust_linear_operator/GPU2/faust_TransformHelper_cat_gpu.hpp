namespace Faust
{
	template<typename FPP>
		TransformHelper<FPP,GPU2>* cat(const std::vector<TransformHelper<FPP,GPU2>*> & THs, int dim /* 0 vert, 1 hor */)
		{
			std::vector<TransformHelper<FPP,Cpu>*> cpuTHs;
			TransformHelper<FPP,Cpu> *cpu_th_out = nullptr;
			for(auto t: THs)
				cpuTHs.push_back(t->tocpu());
			if(dim == 0)
				cpu_th_out = vertcat(cpuTHs);
			else // dim == 1
				cpu_th_out = horzcat(cpuTHs);
			auto gpu_th_out = new TransformHelper<FPP,GPU2>(*cpu_th_out, -1, nullptr);
			delete cpu_th_out;
			for(auto t: cpuTHs)
				delete t;
			return gpu_th_out;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* vertcat(const std::vector<TransformHelper<FPP,GPU2>*> & THs)
		{
			return cat(THs, 0);
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* horzcat(const std::vector<TransformHelper<FPP,GPU2>*> & THs)
		{
			return cat(THs, 1);
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::vertcat(const TransformHelper<FPP,GPU2>* G)
		{
			TransformHelper<FPP,Cpu> th;
			TransformHelper<FPP,Cpu> thG;
			this->tocpu(th);
			G->tocpu(thG);
			auto th_out = th.vertcat(&thG);
			auto gpu_th_out = new TransformHelper<FPP,GPU2>(*th_out, -1, nullptr);
			delete th_out;
			return gpu_th_out;
		}

	template<typename FPP>
		TransformHelper<FPP,GPU2>* TransformHelper<FPP,GPU2>::horzcat(const TransformHelper<FPP,GPU2>* G)
		{
			TransformHelper<FPP,Cpu> th;
			TransformHelper<FPP,Cpu> thG;
			this->tocpu(th);
			G->tocpu(thG);
			auto th_out = th.horzcat(&thG);
			auto gpu_th_out = new TransformHelper<FPP,GPU2>(*th_out, -1, nullptr);
			delete th_out;
			return gpu_th_out;
		}
}
