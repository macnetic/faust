cdef class FaustAlgoCplxDbl:

    @staticmethod
    def fourierFaust(n, norma):
        if(n>31):
            raise ValueError("Faust doesn't handle a DFT of order larger than "
                             "2**31")
        core = FaustCoreGenCplxDbl(core=True)
        core.core_faust_cplx = \
                FaustCoreCy.FaustCoreCpp[complex].fourierFaust(n, norma)
        if(core.core_faust_cplx == NULL):
            raise MemoryError()
        # fourier is always a complex Faust
        return core

