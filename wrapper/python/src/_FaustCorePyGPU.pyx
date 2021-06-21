#TODO: this class should disappear by defining a generic class for both CPU and
# GPU implementations. The easiest way is to use cmake and one or two variables
# (class name, etc.) to generate FaustCore and FaustCoreGPU
cdef class FaustCoreGPU:
    #TODO: refactor with FaustCore using cmake
    #### ATTRIBUTE ########
    # classe Cython
    cdef FaustCoreCy.FaustCoreCppGPU[double]* core_faust_dbl
    cdef FaustCoreCy.FaustCoreCppGPU[complex]* core_faust_cplx
    cdef bool _isReal
    #### CONSTRUCTOR ####
    #def __cinit__(self,np.ndarray[double, mode="fortran", ndim=2] mat):
    def  __cinit__(self,list_factors=None, alpha=1.0, core=False,
                   optimizedCopy=False):
        cdef double [:,:] data
        cdef double [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef complex [:,:] data_cplx
        cdef complex [:] data1d_cplx #only for csr mat factor
        cdef unsigned int nbrow
        cdef unsigned int nbcol
        optimizedCopy=False # TODO: for the moment the auto-conversion of
        #  factors (for opt. of mul by vec) is not enabled, re-enable it later
        # by getting from constructor set by a call from a Faust constructor
        if(alpha != 1.0):
            print("WARNING: the constructor argument for multiplying the Faust"
                  " by a scalar is DEPRECATED and might not be supported in next"
                  " versions of FAµST.")
        if(list_factors is not None):
            if(match('.*int', repr(list_factors[0].dtype))):
               list_factors[0] = list_factors[0].astype(np.float)
            list_factors[0] *= alpha
            self._isReal = True
            for i,factor in enumerate(list_factors):
                # Faust uses row-major order for sparse matrices
                # and col-major order for dense matrices
                # but libmatio uses col-major order for sparse matrices
                if(isinstance(factor, sparse.csc.csc_matrix)):
                    factor = list_factors[i] = factor.tocsr()
                    #print("FaustCorePy.pyx __cinit__(),toarray() factor:", factor)
                if(not isinstance(factor, np.ndarray) and
                   not isinstance(factor, sparse.csr.csr_matrix)):
                   #print("FaustCorePy.pyx __cinit__(), factor:",factor)
                   raise ValueError("Faust factors must be numpy.ndarray or "
                                    " scipy.sparse.csr.csr_matrix")
                if(isinstance(factor[0,0], np.complex)):
                    # str(factor.dtype) == complex128/64/complex_
                    self._isReal = False
                    break
            if(self._isReal):
                self.core_faust_dbl = new FaustCoreCy.FaustCoreCppGPU[double]()
            else:
                self.core_faust_cplx = new FaustCoreCy.FaustCoreCppGPU[complex]()
            for factor in list_factors:
                nbrow=factor.shape[0];
                nbcol=factor.shape[1];
                if(self._isReal):
                   if(isinstance(factor, np.ndarray)):
                      data=factor.astype(float,'F')
                      self.core_faust_dbl.push_back(&data[0,0], nbrow, nbcol,
                                                    optimizedCopy)
                   else: #csr.csr_matrix
                      data1d=factor.data.astype(float,'F')
                      indices=factor.indices.astype(np.int32, 'F')
                      indptr=factor.indptr.astype(np.int32, 'F')
                      self.core_faust_dbl.push_back(&data1d[0], &indptr[0],
                                                    &indices[0], factor.nnz,
                                                    nbrow, nbcol, optimizedCopy)
                else:
                    if(isinstance(factor, sparse.csc.csc_matrix)):
                        #TODO: understand how is it possible to have a sparse
                        # mat here and fix it (because it should have been
                        # converted above already)
                        factor = list_factors[i] = factor.tocsr()
                    #print('FaustCorePy.pyx type factor=',type(factor))
                    #print("FaustCorePy.pyx factor=",factor)
                    if(isinstance(factor, np.ndarray)):
                        data_cplx=factor.astype(np.complex128,'F')
                        self.core_faust_cplx.push_back(&data_cplx[0,0], nbrow,
                                                       nbcol, optimizedCopy)
                    else:
                        #print("FaustCore, factor dims:", nbrow, nbcol)
                        data1d_cplx = factor.data.astype(np.complex128, 'F')
                        indices=factor.indices.astype(np.int32, 'F')
                        indptr=factor.indptr.astype(np.int32, 'F')
                        self.core_faust_cplx.push_back(&data1d_cplx[0], &indptr[0],
                                                    &indices[0], factor.nnz,
                                                       nbrow, nbcol,
                                                       optimizedCopy)
        elif(core): # trick to initialize a new FaustCoreCpp from C++ (see
        # transpose, conj and adjoint)
            pass
        #else:
        #TODO: raise error for undefined object here

    def to_string(self):
        cdef const char* c_str
        if(self._isReal):
            c_str = self.core_faust_dbl.to_string()
        else:
            c_str = self.core_faust_cplx.to_string()
        cdef length = strlen(c_str)
        #printf("%s", c_str[:length])
        #py_str = str(c_str[:length], 'UTF-8')
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        free(<void*>c_str)
        return py_str

    def display(self):
        if(self._isReal):
            self.core_faust_dbl.Display()
        else:
            self.core_faust_cplx.Display()

    def multiply(self,M):
        if(isinstance(M, FaustCoreGPU)):
            return self.multiply_faust(M)
        if not isinstance(M, np.ndarray):
            raise ValueError('input M must be a numpy.ndarray')
        if(self._isReal):
           M=M.astype(float,'F')
           if not M.dtype=='float':
               raise ValueError('input M must be a double array')
        else:
           M=M.astype(complex,'F')
           if(M.dtype not in ['complex', 'complex128', 'complex64'] ): #could fail if complex128 etc.
               raise ValueError('input M must be complex array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError('input M must be Fortran contiguous (column major '
                            'order)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError('input M invalid number of dimensions')
        check_matrix(self._isReal, M)
        ndim_M=M.ndim
        cdef unsigned int nbrow_x=M.shape[0]
        cdef unsigned int nbcol_x #can't be assigned because we don't know yet if the input vector is 1D or 2D

        dimThis=self.shape()
        cdef unsigned int nbRowThis=dimThis[0];
        cdef unsigned int nbColThis=dimThis[1];


        cdef unsigned int nbrow_y=nbRowThis
        cdef unsigned int nbcol_y

        cdef double[:] xview_1D
        cdef double[:,:] xview_2D
        cdef complex[:] xview_1D_cplx
        cdef complex[:,:] xview_2D_cplx

        if ndim_M == 1:
            nbcol_x=1
            if(self._isReal):
                xview_1D=M
            else:
                xview_1D_cplx=M
        else:
            nbcol_x=M.shape[1]
            if(self._isReal):
                xview_2D=M
            else:
                xview_2D_cplx=M

            if (nbrow_x != nbColThis):
                raise ValueError('y=F*M multiplication with Faust: invalid dimension of the input matrix M');

        #void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose);
        nbcol_y = nbcol_x;

        cdef y
        cdef double[:,:] yview
        cdef complex[:,:] yview_cplx
        if(self._isReal):
            y = np.zeros([nbrow_y,nbcol_y], dtype='d',order='F')
            yview = y
        else:
            y = np.zeros([nbrow_y, nbcol_y], dtype='complex', order='F')
            yview_cplx = y

        if ndim_M == 1:
            if(self._isReal):
                self.core_faust_dbl.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_1D[0],nbrow_x,nbcol_x)
            else:
                self.core_faust_cplx.multiply(&yview_cplx[0,0], nbrow_y,
                                              nbcol_y, &xview_1D_cplx[0],
                                              nbrow_x,nbcol_x)
            y = np.squeeze(y) # we want a single dim. (but we created two
            # above)
        else:
            if(self._isReal):
                self.core_faust_dbl.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_2D[0,0],nbrow_x,nbcol_x)
            else:
                self.core_faust_cplx.multiply(&yview_cplx[0,0],nbrow_y,nbcol_y,&xview_2D_cplx[0,0],nbrow_x,nbcol_x)

        return y

    def shape(self):
        cdef unsigned int nbrow = 0
        cdef unsigned int nbcol = 0
        if(self._isReal):
            nbrow = self.core_faust_dbl.getNbRow();
            nbcol = self.core_faust_dbl.getNbCol();
        else:
            nbrow = self.core_faust_cplx.getNbRow();
            nbcol = self.core_faust_cplx.getNbCol();
        return (nbrow,nbcol)

    @staticmethod
    def randFaust(faust_nrows, faust_ncols, t, field, min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density=0.1, per_row=True, gpu=0):
        core = FaustCoreGPU(core=True)
        if(field == 3):
            core.core_faust_dbl = \
                    FaustCoreCy.FaustCoreCppGPU[double].randFaustGPU(faust_nrows,
                                                                     faust_ncols,
                                                                     t, min_num_factors,
                                                                     max_num_factors, min_dim_size,
                                                                     max_dim_size, density, per_row)
            core._isReal = True
            if(core.core_faust_dbl == NULL): raise MemoryError()
        elif(field == 4):
            core.core_faust_cplx = \
                    FaustCoreCy.FaustCoreCppGPU[complex].randFaustGPU(faust_nrows,
                                                                      faust_ncols,
                                                                      t, min_num_factors,
                                                                      max_num_factors, min_dim_size,
                                                                      max_dim_size, density, per_row)
            if(core.core_faust_cplx == NULL): raise MemoryError()
            core._isReal = False
        else:
            raise ValueError("FaustCorePy.randFaust(): field must be 3 for real or"
                             " 4 for complex")

        return core

    @staticmethod
    def hadamardFaust(n, norma):
        if(n>31):
            raise ValueError("Faust doesn't handle a Hadamard of order larger than "
                             "2**31")
        core = FaustCoreGPU(core=True)
        core.core_faust_dbl = \
        FaustCoreCy.FaustCoreCppGPU[double].hadamardFaustGPU(n, norma)
        if(core.core_faust_dbl == NULL):
            raise MemoryError()
        # hadamard is always a real Faust
        core._isReal = True
        return core

    @staticmethod
    def fourierFaust(n, norma):
        if(n>31):
            raise ValueError("Faust doesn't handle a FFT of order larger than "
                             "2**31")
        core = FaustCoreGPU(core=True)
        core.core_faust_cplx = \
                FaustCoreCy.FaustCoreCppGPU[complex].fourierFaustGPU(n, norma)

        if(core.core_faust_cplx == NULL):
            raise MemoryError()

        # fourier is always a complex Faust
        core._isReal = False
        return core

    @staticmethod
    def eyeFaust(n, m, t='real'):
        core = FaustCoreGPU(core=True)
        if(t == 'real'):
            core.core_faust_dbl = \
            FaustCoreCy.FaustCoreCppGPU[double].eyeFaustGPU(n, m)
            core._isReal = True
        elif(t == 'complex'):
            core.core_faust_cplx = \
                    FaustCoreCy.FaustCoreCppGPU[complex].eyeFaustGPU(n,
                                                                  m)
            core._isReal = False
        return core

    def nbytes(self):
        if(self._isReal):
            nbytes = self.core_faust_dbl.getNBytes();
        else:
            nbytes = self.core_faust_cplx.getNBytes();
        return nbytes

    def get_product(self):
        cdef double [:,:] y_data
        cdef complex [:,:] y_data_cplx

        if(self._isReal):
            y_arr = np.empty((self.core_faust_dbl.getNbRow(), self.core_faust_dbl.getNbCol()), dtype='d',
                             order='F')
            y_data = y_arr
            self.core_faust_dbl.get_product(&y_data[0,0], y_arr.shape[0],
                                            y_arr.shape[1])
        else:
            y_arr = np.empty((self.core_faust_cplx.getNbRow(), self.core_faust_cplx.getNbCol()),
                             dtype=np.complex,
                             order='F')
            y_data_cplx = y_arr
            self.core_faust_cplx.get_product(&y_data_cplx[0,0], y_arr.shape[0],
                                            y_arr.shape[1])
        return y_arr

    cdef multiply_faust(self, F):
        if(isinstance(F, FaustCoreGPU)):
            core = FaustCoreGPU(core=True)
            core._isReal = self._isReal
            if(self._isReal):
                if(F.isReal()):
                    core.core_faust_dbl = \
                            self.core_faust_dbl.mul_faust_gpu((<FaustCoreGPU?>F).core_faust_dbl)
                else:
                    core = \
                   <FaustCoreGPU?> (<FaustCoreGPU?> self._ascomplex()).multiply_faust(F)
            else:
                if(F.isReal()):
                    core.core_faust_cplx = \
                             self.core_faust_cplx.mul_faust_gpu((<FaustCoreGPU?>((<FaustCoreGPU?>((<FaustCoreGPU?>F)._ascomplex())))).core_faust_cplx)
                else:
                    core.core_faust_cplx = \
                            self.core_faust_cplx.mul_faust_gpu((<FaustCoreGPU?>F).core_faust_cplx)
            return core
        raise ValueError("F must be a Faust object")

    cdef _ascomplex(self, scalar=1.0):
        cplx_facs = [self.get_fact_opt(i).astype(np.complex) for i in \
                              range(0,self.get_nb_factors())]
        cplx_facs[0] *= scalar # instead of passing the scal to the
        # construc. It avoids disp of
        # deprecation warning
        core = FaustCoreGPU(cplx_facs)#, alpha=scalar)
        core._isReal = False
        return core

    def multiply_scal(self, scalar):
        core = FaustCoreGPU(core=True)
        core._isReal = self._isReal
        if(isinstance(scalar, int)):
            scalar = float(scalar)
        scalar_type_err = TypeError("The mul. scalar must be a real or a"
                                    " complex number")
        if(self._isReal):
            if(isinstance(scalar, float)):
                core.core_faust_dbl = \
                        self.core_faust_dbl.mul_scal_gpu(scalar)
            elif(isinstance(scalar, np.complex)):
                core = <FaustCoreGPU?>self._ascomplex(scalar)
                core._isReal = False
            else:
                raise scalar_type_err
        else:
            if(isinstance(scalar, np.complex) or isinstance(scalar,
                                                            float)):
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_scal_gpu(scalar)
            else:
                raise scalar_type_err
        return core

    def isReal(self):
        return self._isReal

    def nnz(self):
        cdef unsigned long long nnz = 0
        if(self._isReal):
            nnz = self.core_faust_dbl.nnz()
        else:
            nnz = self.core_faust_cplx.nnz()
        return nnz

    def norm(self, ord, **kwargs):
        cdef double norm
        cdef double threshold
        if(str(ord).lower() not in ["1","2","fro", "inf"]):
            raise ValueError("FaustCorePy.norm() invalid type of norm asked.")
        threshold = .001
        max_num_its = 100
        if('threshold' in kwargs.keys()):
            threshold = kwargs['threshold']
        if('max_num_its' in kwargs.keys()):
            max_num_its = kwargs['max_num_its']
        if(self._isReal):
            if(isinstance(ord,int)):
                norm = self.core_faust_dbl.norm(ord, threshold, max_num_its)
            elif(ord == np.inf):
                norm = self.core_faust_dbl.normInf()
            else:
                norm = self.core_faust_dbl.normFro()
        else:
            if(isinstance(ord,int)):
                norm = self.core_faust_cplx.norm(ord, threshold, max_num_its)
            elif(ord == np.inf):
                norm = self.core_faust_cplx.normInf()
            else:
                norm = self.core_faust_cplx.normFro()
        return norm

    def power_iteration(self, threshold, max_num_its):
        cdef double[:] _lambda_dbl_view
        cdef complex[:]  _lambda_cplx_view
        if(self._isReal):
            _lambda = np.empty((1,), dtype='float')
            _lambda_dbl_view = _lambda
            self.core_faust_dbl.power_iteration(&_lambda_dbl_view[0], threshold, max_num_its)
        else:
            _lambda = np.empty((1,), dtype=np.complex)
            _lambda_cplx_view = _lambda
            self.core_faust_cplx.power_iteration(&_lambda_cplx_view[0], threshold, max_num_its)
        return _lambda[0]

    def normalize(self, ord):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.normalize_gpu(ord)
        else:
            core.core_faust_cplx = self.core_faust_cplx.normalize_gpu(ord)
        core._isReal = self._isReal
        return core

    def get_nb_factors(self):
        cdef int nb_factors
        if(self._isReal):
            nb_factors = int(self.core_faust_dbl.get_nb_factors())
        else:
            nb_factors = int(self.core_faust_cplx.get_nb_factors())
        return nb_factors

    def get_fact_opt(self, i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef double[:,:] fact_dbl_view
        cdef complex[:,:] fact_cplx_view
        cdef rowptr, col_ids, elts, nnz
        cdef int[:] rowptr_view, col_ids_view
        if(self._isReal):
            dtype = 'd'
            is_fact_sparse = self.core_faust_dbl.is_fact_sparse(i)
            is_transposed = self.core_faust_dbl.isTransposed()
        else:
            dtype = 'complex'
            is_fact_sparse = self.core_faust_cplx.is_fact_sparse(i)
            is_transposed = self.core_faust_cplx.isTransposed()
        if(is_fact_sparse):
            if(self._isReal):
                # is_transposed = False # uncomment to disable the trick which
                # uses csc representation instead of csr transpose
                # to optimize copy
                nnz = self.core_faust_dbl.get_fact_nnz(i)
                col_ids = np.ndarray([nnz], dtype=np.int32)
                elts = np.ndarray([1,nnz], dtype=dtype)
                col_ids_view = col_ids
                shape = [self.core_faust_dbl.get_fact_nb_rows(i),
                 self.core_faust_dbl.get_fact_nb_cols(i)]
                if(is_transposed):
                    rowptr_sz = shape[1]+1
                else:
                    rowptr_sz = shape[0]+1

                rowptr = np.ndarray([rowptr_sz], dtype=np.int32)
                rowptr_view = rowptr
                fact_dbl_view = elts
    #                print("len(rowptr)=", len(rowptr))
    #                print("len(col_ids)=", len(col_ids))
    #                print("len(elts[0]=", len(elts[0,:]))
                self.core_faust_dbl.get_fact_sparse(i, &rowptr_view[0],
                                                        &col_ids_view[0],
                                                        &fact_dbl_view[0,0],
                                                       is_transposed)
            else:
                # is_transposed = False # uncomment to disable the trick which
                # uses csc representation instead of csr transpose
                # to optimize copy
                nnz = self.core_faust_cplx.get_fact_nnz(i)
                col_ids = np.ndarray([nnz], dtype=np.int32)
                elts = np.ndarray([1,nnz], dtype=dtype)
                col_ids_view = col_ids
                shape = [self.core_faust_cplx.get_fact_nb_rows(i),
                 self.core_faust_cplx.get_fact_nb_cols(i)]
                if(is_transposed):
                    rowptr_sz = shape[1]+1
                else:
                    rowptr_sz = shape[0]+1

                rowptr = np.ndarray([rowptr_sz], dtype=np.int32)
                rowptr_view = rowptr
                fact_cplx_view = elts
                self.core_faust_cplx.get_fact_sparse(i, &rowptr_view[0],
                                                     &col_ids_view[0],
                                                     &fact_cplx_view[0,0],
                                                     is_transposed)
                #print("(rowptr)=", (rowptr))
#                print("(col_ids)=", (col_ids))
#                print("(elts[0,:]=", (elts))
            if(is_transposed):
                fact = csc_matrix((elts[0,:], col_ids, rowptr), shape=shape)
            else:
                fact = csr_matrix((elts[0,:], col_ids, rowptr), shape=shape)
        else: # dense matrix
            if(is_transposed):
                order = 'C'
                # C contiguous repr. (row-major order ) is used to optimized the
                # request of transpose factor (no need to reorder data as we
                # should in col-major/Fortran repr.)
            else:
                order = 'F'
            if(self._isReal):
                fact = np.ndarray([self.core_faust_dbl.get_fact_nb_rows(i),
                                   self.core_faust_dbl.get_fact_nb_cols(i)], dtype=dtype,
                                  order=order)
                fact_dbl_view = fact
                self.core_faust_dbl.get_fact_dense(i, &fact_dbl_view[0, 0],
                                               <unsigned int*>NULL,
                                               <unsigned int*>NULL,
                                               is_transposed)
            else:
                fact = np.ndarray([self.core_faust_cplx.get_fact_nb_rows(i),
                             self.core_faust_cplx.get_fact_nb_cols(i)], dtype=dtype,
                            order=order)
                fact_cplx_view = fact
                self.core_faust_cplx.get_fact_dense(i, &fact_cplx_view[0, 0],
                                                   <unsigned int*>NULL,
                                                   <unsigned int*>NULL,
                                                   is_transposed)

        return fact

    def left(self, id):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.left_gpu(id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.left_gpu(id)
        core._isReal = self._isReal
        return core

    def right(self, id):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.right_gpu(id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.right_gpu(id)
        core._isReal = self._isReal
        return core

    def transpose(self):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.transpose_gpu()
        else:
            core.core_faust_cplx = self.core_faust_cplx.transpose_gpu()
        core._isReal = self._isReal
        return core

    def conj(self):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.conjugate_gpu()
        else:
            core.core_faust_cplx = self.core_faust_cplx.conjugate_gpu()
        core._isReal = self._isReal
        return core

    def zpruneout(self, nnz_tres, npasses, only_forward):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.zpruneout_gpu(nnz_tres,
                                                                npasses,
                                                                only_forward)
        else:
            core.core_faust_cplx = self.core_faust_cplx.zpruneout_gpu(nnz_tres,
                                                                  npasses,
                                                                  only_forward)
        core._isReal = self._isReal
        return core

    def getH(self):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.adjoint_gpu()
        else:
            core.core_faust_cplx = self.core_faust_cplx.adjoint_gpu()
        core._isReal = self._isReal
        return core

    def clone(self, dev='cpu'):
        core_gpu = FaustCoreGPU(core=True)
        core_cpu = FaustCore(core=True)
        if(dev.startswith('gpu')):
            if(self._isReal):
                core_gpu.core_faust_dbl = self.core_faust_dbl.clone_gpu()
            else:
                core_gpu.core_faust_cplx = self.core_faust_cplx.clone_gpu()
            core_gpu._isReal = self._isReal
            return core_gpu
        elif(dev == 'cpu'):
            if(self._isReal):
                core_cpu.core_faust_dbl = self.core_faust_dbl.clone_cpu()
                core_cpu._isReal = self._isReal
            else:
                core_cpu.core_faust_cplx = self.core_faust_cplx.clone_cpu()
            return core_cpu
        else:
            raise ValueError('dev='+str(dev)+' is not a valid device')

    def vertcat(self,F):
        if(F.isReal() and not self.isReal()):
            return self._vertcat((<FaustCoreGPU?>F)._ascomplex())
        return self._vertcat(F)

    cdef _horzcat(self, F):
         #TODO: refactor with _vertcat(), maybe by passing func ref to a cat func
         core = FaustCoreGPU(core=True)
         #TODO/ F must be a FaustCore
         if(self._isReal):
             if(not F.isReal()):
                 self = self._ascomplex()
                 core.core_faust_cplx = \
                 self.core_faust_cplx.horzcat_gpu((<FaustCoreGPU?>F).core_faust_cplx)
             else:
                 core.core_faust_dbl = \
                         self.core_faust_dbl.horzcat_gpu((<FaustCoreGPU?>F).core_faust_dbl)
         else:
             core.core_faust_cplx = \
                     self.core_faust_cplx.horzcat_gpu((<FaustCoreGPU?>F).core_faust_cplx)
         core._isReal = core.core_faust_dbl != NULL
         return core

    def horzcat(self,F):
        if(F.isReal() and not self.isReal()):
            return self._horzcat((<FaustCoreGPU?>F)._ascomplex())
        return self._horzcat(F)

    cdef _vertcat(self, F):
         core = FaustCoreGPU(core=True)
         #TODO/ F must be a FaustCore
         if(self._isReal):
             if(not F.isReal()):
                 self = self._ascomplex()
                 core.core_faust_cplx = \
                 self.core_faust_cplx.vertcat_gpu((<FaustCoreGPU?>F).core_faust_cplx)
             else:
                 core.core_faust_dbl = \
                         self.core_faust_dbl.vertcat_gpu((<FaustCoreGPU?>F).core_faust_dbl)
         else:
             core.core_faust_cplx = \
                     self.core_faust_cplx.vertcat_gpu((<FaustCoreGPU?>F).core_faust_cplx)
         core._isReal = core.core_faust_dbl != NULL
         return core

    cdef _vertcatn(self, Fs):
        cdef FaustCoreCy.FaustCoreCppGPU[double]** _Fs
        cdef FaustCoreCy.FaustCoreCppGPU[complex]** _Fs_cplx

        if self.isReal():
            _Fs = <FaustCoreCy.FaustCoreCppGPU[double]**> PyMem_Malloc(sizeof(void*) *
                                                               len(Fs))
        else:
            _Fs_cplx = <FaustCoreCy.FaustCoreCppGPU[complex]**> PyMem_Malloc(sizeof(void*) *
                                                               len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                _Fs[i] = (<FaustCoreGPU?>F).core_faust_dbl
            else:
                _Fs_cplx[i] = (<FaustCoreGPU?>F).core_faust_cplx
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.vertcatn_gpu(_Fs, len(Fs))
        else:
            core.core_faust_cplx = self.core_faust_cplx.vertcatn_gpu(_Fs_cplx,
                                                                 len(Fs))
        core._isReal = core.core_faust_dbl != NULL
        return core

    def vertcatn(self, *args):
        Fs = []
        i = 0
        any_complex = not self.isReal()
        while i < len(args) and not any_complex:
            any_complex = not args[i].isReal()
            i+=1
        for F in args:
            if F.isReal() and any_complex:
               F = (<FaustCoreGPU?>F)._ascomplex()
            Fs += [F]
        if any_complex and self.isReal():
            self = (<FaustCoreGPU?>self)._ascomplex()
        return self._vertcatn(Fs)

    cdef _horzcatn(self, Fs):
        cdef FaustCoreCy.FaustCoreCppGPU[double]** _Fs
        cdef FaustCoreCy.FaustCoreCppGPU[complex]** _Fs_cplx

        if self.isReal():
            _Fs = <FaustCoreCy.FaustCoreCppGPU[double]**> PyMem_Malloc(sizeof(void*) *
                                                               len(Fs))
        else:
            _Fs_cplx = <FaustCoreCy.FaustCoreCppGPU[complex]**> PyMem_Malloc(sizeof(void*) *
                                                               len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                _Fs[i] = (<FaustCoreGPU?>F).core_faust_dbl
            else:
                _Fs_cplx[i] = (<FaustCoreGPU?>F).core_faust_cplx
        core = FaustCoreGPU(core=True)
        # print("self.isReal():", self.isReal())
        if self.isReal():
            core.core_faust_dbl = self.core_faust_dbl.horzcatn_gpu(_Fs, len(Fs))
        else:
            core.core_faust_cplx = self.core_faust_cplx.horzcatn_gpu(_Fs_cplx,
                                                                 len(Fs))
        core._isReal = core.core_faust_dbl != NULL
        return core

    def horzcatn(self, *args):
        Fs = []
        i = 0
        any_complex = not self.isReal()
        # print("any_complex:", any_complex)
        while i < len(args) and not any_complex:
            any_complex = not args[i].isReal()
            i+=1
        for F in args:
            if F.isReal() and any_complex:
               F = (<FaustCoreGPU?>F)._ascomplex()
            Fs += [F]
        # print("any_complex:", any_complex)
        if any_complex and self.isReal():
            self = (<FaustCoreGPU?>self)._ascomplex()
        return self._horzcatn(Fs)

    def swap_cols(self, id1, id2, permutation, inplace):
        if(inplace):
            if(self._isReal):
                self.core_faust_dbl.swap_cols_gpu(id1, id2,
                                              permutation, inplace)
            else:
                self.core_faust_cplx.swap_cols_gpu(id1, id2,
                                               permutation,
                                               inplace)

            return self
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.swap_cols_gpu(id1, id2,
                                                                permutation, inplace)
        else:
            core.core_faust_cplx = self.core_faust_cplx.swap_cols_gpu(id1, id2,
                                                                  permutation,
                                                                  inplace)
        core._isReal = self._isReal
        return core

    def swap_rows(self, id1, id2, permutation, inplace):
        if(inplace):
            if(self._isReal):
                self.core_faust_dbl.swap_rows_gpu(id1, id2,
                                              permutation, inplace)
            else:
                self.core_faust_cplx.swap_rows_gpu(id1, id2,
                                               permutation,
                                               inplace)

            return self
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.swap_rows_gpu(id1, id2,
                                                                permutation, inplace)
        else:
            core.core_faust_cplx = self.core_faust_cplx.swap_rows_gpu(id1, id2,
                                                                  permutation,
                                                                  inplace)
        core._isReal = self._isReal
        return core

    def device(self):
        cdef char c_str[256]
        if(self._isReal):
            self.core_faust_dbl.device_gpu(c_str)
        else:
            self.core_faust_cplx.device_gpu(c_str)
        cdef length = strlen(c_str)
        #printf("%s", c_str[:length])
        #py_str = str(c_str[:length], 'UTF-8')
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        return py_str

    def slice(self, indices):
        # TODO: rename this function or cut in two: slice and fancy indexing
        core = FaustCoreGPU(core=True)
        start_row_id, end_row_id, start_col_id, end_col_id = (indices[0].start,
                                                              indices[0].stop,
                                                              indices[1].start,
                                                              indices[1].stop)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.slice_gpu(start_row_id,
                                                            end_row_id,
                                                            start_col_id,
                                                            end_col_id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.slice_gpu(start_row_id,
                                                              end_row_id,
                                                              start_col_id,
                                                              end_col_id)

        core._isReal = self._isReal
        return core

    def fancy_idx(self, indices):
        cdef unsigned long int[:] row_indices_view
        cdef unsigned long int[:] col_indices_view
        core = FaustCoreGPU(core=True)
        # fancy indexing
        #convert possible slice (on index 0 or 1 of out_indices) to
        # an array of indices
        for i in range(0,2):
            if(isinstance(indices[i], slice)):
               indices[i] = list(range(indices[i].start, indices[i].stop,
                                       indices[i].step))
        # it's possible that on certain architectures unsigned long int is a
        # 4-bytes integer
        # TODO: move to a cross-platform size type (like uint32/64_t from stdlib.h)
        if(sizeof(unsigned long int) == 8):
            dtype = np.uint64
        elif(sizeof(unsigned long int) == 4):
            dtype = np.uint32
        row_indices = np.array(indices[0], dtype=dtype)
        col_indices = np.array(indices[1], dtype=dtype)
        row_indices_view = row_indices
        col_indices_view = col_indices
#        print("FaustCorePy.fancy_idx(), row_indices=", row_indices, " size=",
#              row_indices.size)
#        print("FaustCorePy.fancy_idx(), col_indices=", col_indices," size=",
#              col_indices.size)
        if(self._isReal):
            core.core_faust_dbl = \
            self.core_faust_dbl.fancy_idx_gpu(&row_indices_view[0], row_indices.size,
                                          &col_indices_view[0], col_indices.size)
        else:
            core.core_faust_cplx = \
            self.core_faust_cplx.fancy_idx_gpu(&row_indices_view[0],
                                           row_indices.size,
                                           &col_indices_view[0],
                                           col_indices.size)

        core._isReal = self._isReal
        return core

    def save_mat_file(self,filepath):
        cdef char * cfilepath = <char*> PyMem_Malloc(sizeof(char) *
                                                     (len(filepath)+1))
        fparr = bytearray(filepath, "UTF-8");
        for i in range(0,len(filepath)):
            cfilepath[i] = fparr[i]
        cfilepath[i+1] = 0
        if(self._isReal):
            ret = self.core_faust_dbl.save_mat_file(cfilepath)
        else:
            ret = self.core_faust_cplx.save_mat_file(cfilepath)
        if(not ret):
            raise Exception("Failed to save the file: "+filepath)
        PyMem_Free(cfilepath)

    def multiply_csr_mat(self, X):
        cdef double [:] x_data1d
        cdef int [:] x_indices
        cdef int [:] x_indptr
        cdef complex [:] x_data_cplx
        cdef double [:,:] y_data
        cdef complex [:,:] y_data_cplx
        # X is supposed to be a csr_matrix
        x_indices = X.indices
        x_indptr = X.indptr
        x_nnz = X.nnz
        nbcol = X.shape[1]
        e = Exception("Dimensions must agree")
        if(X.dtype in [ 'double', 'float64']):
            x_data1d = X.data
            nbrow = self.core_faust_dbl.getNbRow()
            if(self.core_faust_dbl.getNbCol() != X.shape[0]): raise e
            y_data_arr = np.empty((nbrow,nbcol), dtype=np.double, order='F') # we don't know beforehand Y nnz
            y_data = y_data_arr
            if(self._isReal):
                self.core_faust_dbl.multiply_gpu(&y_data[0,0], nbrow, nbcol,
                                                 &x_data1d[0], &x_indptr[0],
                                                 &x_indices[0],
                                                 x_nnz, X.shape[0], X.shape[1])
            else:
                raise("x must be real if Faust is.")
                # shouldn't happen normally (avoided by calling function)
        else:
            x_data_cplx = X.data
            nbrow = self.core_faust_cplx.getNbRow()
            if(self.core_faust_cplx.getNbCol() != X.shape[0]): raise e
            y_data_arr = np.empty((nbrow,nbcol), dtype=np.complex, order='F') # we don't know beforehand Y nnz
            y_data_cplx = y_data_arr
            if(not self._isReal):
                self.core_faust_cplx.multiply_gpu(&y_data_cplx[0,0], nbrow, nbcol,
                                          &x_data_cplx[0], &x_indptr[0],
                                          &x_indices[0],
                                          x_nnz, X.shape[0], X.shape[1])
            else:
                raise("x must be complex if Faust is")
                # shouldn't happen normally (avoided by calling function)
        return y_data_arr

    def optimize_storage(self, time=False):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.optimize_storage_gpu(time)
        else:
            core.core_faust_cplx = self.core_faust_cplx.optimize_storage_gpu(time)
        core._isReal = self._isReal
        return core

    def optimize(self, transp=False):
        core = FaustCoreGPU(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.optimize_gpu(transp)
        else:
            core.core_faust_cplx = self.core_faust_cplx.optimize_gpu(transp)
        core._isReal = self._isReal
        return core

    def optimize_time(self, transp=False, inplace=False, nsamples=1):
        if(inplace):
            if(self._isReal):
                self.core_faust_dbl.optimize_time_gpu(transp, inplace, nsamples)
            else:
                self.core_faust_cplx.optimize_time(transp, inplace, nsamples)
        else:
            core = FaustCoreGPU(core=True)
            if(self._isReal):
                core.core_faust_dbl = self.core_faust_dbl.optimize_time_gpu(transp,
                                                                      inplace,
                                                                      nsamples)
            else:
                core.core_faust_cplx = self.core_faust_cplx.optimize_time_gpu(transp,
                                                                         inplace,
                                                                         nsamples)
            core._isReal = self._isReal
            return core


