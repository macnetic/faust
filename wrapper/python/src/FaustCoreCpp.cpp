#include "FaustCoreCpp.h"

void* _enable_gpu_mod(const char* libpath, const bool silent)
{
#ifdef USE_GPU_MOD
    return Faust::enable_gpu_mod(libpath, silent);
#else
    return nullptr;
#endif
}

bool _is_gpu_mod_enabled()
{
#ifdef USE_GPU_MOD
    return Faust::is_gpu_mod_enabled();
#else
    return false;
#endif
}
