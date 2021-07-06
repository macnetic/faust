# pxd for functions available/making sense only for GPU
cdef extern from "FaustCoreCpp.h":
    FaustCoreCppCPU[FPP]* clone_gpu2cpu[FPP](FaustCoreCppGPU[FPP]*)
    FaustCoreCppGPU[FPP]* clone_cpu2gpu[FPP](FaustCoreCppCPU[FPP]*)
