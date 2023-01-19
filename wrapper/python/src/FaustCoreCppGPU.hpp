#ifdef USE_GPU_MOD



template<typename FPP>
FaustCoreCpp<FPP, Cpu>* clone_gpu2cpu(FaustCoreCpp<FPP, GPU2>* gpu_t)
{
    Faust::TransformHelper<FPP,Cpu>* th = gpu_t->transform->tocpu();
    FaustCoreCpp<FPP, Cpu>* core = new FaustCoreCpp<FPP, Cpu>(th);
    return core;
}

template<typename FPP>
FaustCoreCpp<FPP, GPU2>* clone_cpu2gpu(FaustCoreCpp<FPP, Cpu>* cpu_t)
{
	auto th = new Faust::TransformHelper<FPP,GPU2>(*cpu_t->transform);
	FaustCoreCpp<FPP, GPU2>* core = new FaustCoreCpp<FPP, GPU2>(th);
	return core;
}
#endif
