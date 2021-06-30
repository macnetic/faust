cdef class FaustCoreCplx(FaustCore):

    #cdef FaustCoreCy.FaustCoreCpp[complex]* core_faust_cplx # TODO: uncomment
    # when all uses in parent are moved here

    def  __cinit__(self,list_factors=None, alpha=1.0, core=False,
                   optimizedCopy=False):
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef complex [:,:] data_cplx
        cdef complex [:] data1d_cplx #only for csr mat factor
        cdef unsigned int nbrow
        cdef unsigned int nbcol
        optimizedCopy=False # TODO: so far the auto-conversion of
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
            self.core_faust_cplx = new FaustCoreCy.FaustCoreCpp[complex]()
            for factor in list_factors:
                nbrow=factor.shape[0];
                nbcol=factor.shape[1];
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

    def clone(self, *args, **kwargs):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.clone()
        return core

    def nbytes(self):
        nbytes = self.core_faust_cplx.getNBytes();
        return nbytes

    def shape(self):
        cdef unsigned int nbrow = 0
        cdef unsigned int nbcol = 0
        nbrow = self.core_faust_cplx.getNbRow();
        nbcol = self.core_faust_cplx.getNbCol();
        return (nbrow,nbcol)

    def multiply_csr_mat(self, X):
        cdef int [:] x_indices
        cdef int [:] x_indptr
        cdef complex [:] x_data_cplx
        cdef complex [:,:] y_data_cplx
        # X is supposed to be a csr_matrix
        x_indices = X.indices
        x_indptr = X.indptr
        x_nnz = X.nnz
        nbcol = X.shape[1]
        e = Exception("Dimensions must agree")
        x_data_cplx = X.data
        nbrow = self.core_faust_cplx.getNbRow()
        if(self.core_faust_cplx.getNbCol() != X.shape[0]): raise e
        y_data_arr = np.empty((nbrow,nbcol), dtype=np.complex, order='F') # we don't know beforehand Y nnz
        y_data_cplx = y_data_arr
        self.core_faust_cplx.multiply(&y_data_cplx[0,0], nbrow, nbcol,
                                  &x_data_cplx[0], &x_indptr[0],
                                  &x_indices[0],
                                  x_nnz, X.shape[0], X.shape[1])
        return y_data_arr

    def get_product(self):
        cdef complex [:,:] y_data_cplx

        y_arr = np.empty((self.core_faust_cplx.getNbRow(), self.core_faust_cplx.getNbCol()),
                         dtype=np.complex,
                         order='F')
        y_data_cplx = y_arr
        self.core_faust_cplx.get_product(&y_data_cplx[0,0], y_arr.shape[0],
                                        y_arr.shape[1])
        return y_arr

    def multiply_faust(self, F):
        if(isinstance(F, FaustCore)):
            core = FaustCoreCplx(core=True)
            if(F.isReal()):
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_faust((<FaustCoreCplx?>((<FaustCoreCplx?>((<FaustCore?>F)._ascomplex())))).core_faust_cplx)
            else:
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_faust((<FaustCoreCplx?>F).core_faust_cplx)
            return core
        raise ValueError("F must be a Faust object")

    def isReal(self):
        return False

    def device(self):
        # always returns cpu but calls the cpp code just in case
        cdef char c_str[256]
        self.core_faust_cplx.device(c_str)
        cdef length = strlen(c_str)
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        return py_str

    def multiply_scal(self, scalar):
        core = FaustCoreCplx(core=True)
        if(isinstance(scalar, int)):
            scalar = float(scalar)
        scalar_type_err = TypeError("The mul. scalar must be a real or a"
                                    " complex number")
        if(isinstance(scalar, np.complex) or isinstance(scalar,
                                                        float)):
            core.core_faust_cplx = \
                    self.core_faust_cplx.mul_scal(scalar)
        else:
            raise scalar_type_err
        return core

    cdef _vertcat(self, F):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = \
                self.core_faust_cplx.vertcat((<FaustCoreCplx?>F).core_faust_cplx)
        return core

    def vertcat(self,F):
        # TODO: merge with _vertcat 
        if F.isReal():
            return self._vertcat((<FaustCoreCplx?>F)._ascomplex())
        return self._vertcat(F)

    cdef _horzcat(self, F):
         #TODO: refactor with _vertcat(), maybe by passing func ref to a cat func
         core = FaustCoreCplx(core=True)
         #TODO/ F must be a FaustCore
         core.core_faust_cplx = \
                 self.core_faust_cplx.horzcat((<FaustCoreCplx?>F).core_faust_cplx)
         return core

    def horzcat(self,F):
        if F.isReal():
            return self._horzcat((<FaustCoreCplx?>F)._ascomplex())
        return self._horzcat(F)

    cdef _vertcatn(self, Fs):
        cdef FaustCoreCy.FaustCoreCpp[complex]** _Fs_cplx

        _Fs_cplx = <FaustCoreCy.FaustCoreCpp[complex]**> PyMem_Malloc(sizeof(void*) *
                                                           len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                # can't happen if vertcatn is the caller
                raise TypeError("cat Faust with different dtype")
            else:
                _Fs_cplx[i] = (<FaustCoreCplx?>F).core_faust_cplx

        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.vertcatn(_Fs_cplx,
                                                                 len(Fs))
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
               F = (<FaustCore?>F)._ascomplex()
            Fs += [F]
        if any_complex and self.isReal():
            self = (<FaustCore?>self)._ascomplex()
        return self._vertcatn(Fs)

    cdef _horzcatn(self, Fs):
        cdef FaustCoreCy.FaustCoreCpp[complex]** _Fs_cplx

        _Fs_cplx = <FaustCoreCy.FaustCoreCpp[complex]**> PyMem_Malloc(sizeof(void*) *
                                                               len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                # can't happen if vertcatn is the caller
                raise TypeError("cat Faust with different dtype")
            else:
                _Fs_cplx[i] = (<FaustCoreCplx?>F).core_faust_cplx
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.horzcatn(_Fs_cplx,
                                                             len(Fs))
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
               F = (<FaustCore?>F)._ascomplex()
            Fs += [F]
        # print("any_complex:", any_complex)
        if any_complex and self.isReal():
            self = (<FaustCore?>self)._ascomplex()
        return self._horzcatn(Fs)

    def polyCoeffs(self, coeffs):
        core = FaustCoreCplx(core=True)
        cdef complex[:] cview_cplx
        cview_cplx = coeffs
        core.core_faust_cplx = self.core_faust_cplx.polyCoeffs(&cview_cplx[0])
        return core

    def mulPolyCoeffs(self, coeffs, X, out=None):
        cdef complex[:] cview_cplx
        cdef complex[:,:] Xview_cplx
        cdef complex[:,:] Yview_cplx

        #TODO: all arguments must be of the same scalar type as self

        X_1dim = False
        d = X.shape[0]
        if(X.ndim > 1):
            X = np.asfortranarray(X)
            n = X.shape[1]
        else:
            n = 1
            X = X.reshape(d, 1)
            X_1dim = True

        if not isinstance(out, type(None)):
            Y = out
            if Y.ndim == 1:
                raise ValueError('out must have 2 dimensions.')
            if Y.shape != (d,n):
                raise ValueError('out shape isn\'t valid.')
            dtype_err = ValueError('out dtype isn\'t valid.')
            if not Y.flags['F_CONTIGUOUS']:
                raise ValueError('the array must be in fortran/column continous order.')
            if Y.dtype != 'complex':
                raise dtype_err
            Yview_cplx = Y
            Xview_cplx = X
            cview_cplx = coeffs
        else:
            cview_cplx = coeffs
            Y = np.empty((X.shape[0], n), dtype=np.complex, order='F')
            Yview_cplx = Y
            Xview_cplx = X


        self.core_faust_cplx.mulPolyCoeffs(&Xview_cplx[0,0],
                                           n,
                                           &Yview_cplx[0,0],
                                           &cview_cplx[0])
        if X_1dim:
            # X is a vector, Y must be one too
            Y = np.squeeze(Y)

        return Y

    def polyNext(self):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.polyNext()
        return core

    # Left-Multiplication by a Faust F
    # y=multiply(F,M) is equivalent to y=F*M
    def multiply(self,M):
        if isinstance(M, FaustCore):
            return self.multiply_faust(M)
        if not isinstance(M, np.ndarray):
            raise ValueError('input M must be a numpy.ndarray')
        if(self.isReal()):
           M=M.astype(float,'F')
           if not M.dtype=='float':
               raise ValueError('input M must be a double array')
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
        check_matrix(self.isReal(), M)
        ndim_M=M.ndim
        cdef unsigned int nbrow_x=M.shape[0]
        cdef unsigned int nbcol_x #can't be assigned because we don't know yet if the input vector is 1D or 2D

        dimThis=self.shape()
        cdef unsigned int nbRowThis=dimThis[0];
        cdef unsigned int nbColThis=dimThis[1];


        cdef unsigned int nbrow_y=nbRowThis
        cdef unsigned int nbcol_y

        cdef complex[:] xview_1D_cplx
        cdef complex[:,:] xview_2D_cplx

        if ndim_M == 1:
            nbcol_x=1
            xview_1D_cplx=M
        else:
            nbcol_x=M.shape[1]
            xview_2D_cplx=M

            if (nbrow_x != nbColThis):
                raise ValueError('y=F*M multiplication with Faust: invalid dimension of the input matrix M');

        #void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose);
        nbcol_y = nbcol_x;

        cdef y
        cdef complex[:,:] yview_cplx
        y = np.empty([nbrow_y, nbcol_y], dtype='complex', order='F')
        yview_cplx = y

        if ndim_M == 1:
            self.core_faust_cplx.multiply(&yview_cplx[0,0], nbrow_y,
                                          nbcol_y, &xview_1D_cplx[0],
                                          nbrow_x,nbcol_x)
            y = np.squeeze(y) # we want a single dim. (but we created two
            # above)
        else:
            self.core_faust_cplx.multiply(&yview_cplx[0,0],nbrow_y,nbcol_y,&xview_2D_cplx[0,0],nbrow_x,nbcol_x)

        return y

    def set_FM_mul_mode(self, mode):
        self.core_faust_cplx.set_FM_mul_mode(mode)

    def set_Fv_mul_mode(self, mode):
        self.core_faust_cplx.set_Fv_mul_mode(mode)


    # print information about the faust (size, number of factor, type of factor (dense/sparse) ...)
    def display(self):
        self.core_faust_cplx.Display()

    def to_string(self):
        cdef const char* c_str
        c_str = self.core_faust_cplx.to_string()
        cdef length = strlen(c_str)
        #printf("%s", c_str[:length])
        #py_str = str(c_str[:length], 'UTF-8')
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        free(<void*>c_str)
        return py_str

    def nnz(self):
        cdef unsigned long long nnz = 0
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
        _lambda = np.empty((1,), dtype=np.complex)
        _lambda_cplx_view = _lambda
        self.core_faust_cplx.power_iteration(&_lambda_cplx_view[0], threshold, max_num_its)
        return _lambda[0]

    def normalize(self, ord):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.normalize(ord)
        return core

    def get_nb_factors(self):
        cdef int nb_factors
        nb_factors = int(self.core_faust_cplx.get_nb_factors())
        return nb_factors

    def get_fact(self,i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef complex[:,:] fact_cplx_view
        fact = np.empty([self.core_faust_cplx.get_fact_nb_rows(i),
                          self.core_faust_cplx.get_fact_nb_cols(i)],
                             dtype='complex', order='F')
        fact_cplx_view = fact
        self.core_faust_cplx.get_fact(i, &fact_cplx_view[0, 0])
        return fact

    def get_fact_opt(self, i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef complex[:,:] fact_cplx_view
        cdef rowptr, col_ids, elts, nnz
        cdef int[:] rowptr_view, col_ids_view
        dtype = 'complex'
        is_fact_sparse = self.core_faust_cplx.is_fact_sparse(i)
        is_transposed = self.core_faust_cplx.isTransposed()
        if(is_fact_sparse):
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
            fact = np.ndarray([self.core_faust_cplx.get_fact_nb_rows(i),
                         self.core_faust_cplx.get_fact_nb_cols(i)], dtype=dtype,
                        order=order)
            fact_cplx_view = fact
            self.core_faust_cplx.get_fact_dense(i, &fact_cplx_view[0, 0],
                                               <unsigned int*>NULL,
                                               <unsigned int*>NULL,
                                               is_transposed)
        return fact

    def slice(self, indices):
        # TODO: rename this function or cut in two: slice and fancy indexing
        core = FaustCoreCplx(core=True)
        start_row_id, end_row_id, start_col_id, end_col_id = (indices[0].start,
                                                              indices[0].stop,
                                                              indices[1].start,
                                                              indices[1].stop)
        core.core_faust_cplx = self.core_faust_cplx.slice(start_row_id,
                                                          end_row_id,
                                                          start_col_id,
                                                          end_col_id)

        return core

    def left(self, id):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.left(id)
        return core

    def right(self, id):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.right(id)
        return core

    def fancy_idx(self, indices):
        cdef unsigned long int[:] row_indices_view
        cdef unsigned long int[:] col_indices_view
        core = FaustCoreCplx(core=True)
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
        core.core_faust_cplx = \
        self.core_faust_cplx.fancy_idx(&row_indices_view[0],
                                       row_indices.size,
                                       &col_indices_view[0],
                                       col_indices.size)

        return core

    def save_mat_file(self,filepath):
        cdef char * cfilepath = <char*> PyMem_Malloc(sizeof(char) *
                                                     (len(filepath)+1))
        fparr = bytearray(filepath, "UTF-8");
        for i in range(0,len(filepath)):
            cfilepath[i] = fparr[i]
        cfilepath[i+1] = 0
        ret = self.core_faust_cplx.save_mat_file(cfilepath)
        if(not ret):
            raise Exception("Failed to save the file: "+filepath)
        PyMem_Free(cfilepath)

    def transpose(self):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.transpose()
        return core

    def swap_cols(self, id1, id2, permutation, inplace):
        if(inplace):
            self.core_faust_cplx.swap_cols(id1, id2,
                                           permutation,
                                           inplace)

            return self
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.swap_cols(id1, id2,
                                                              permutation,
                                                              inplace)
        return core

    def swap_rows(self, id1, id2, permutation, inplace):
        if(inplace):
            self.core_faust_cplx.swap_rows(id1, id2,
                                           permutation,
                                           inplace)

            return self
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.swap_rows(id1, id2,
                                                              permutation,
                                                              inplace)
        return core

    def optimize_storage(self, time=False):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.optimize_storage(time)
        return core

    def optimize(self, transp=False):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.optimize(transp)
        return core

    def optimize_time(self, transp=False, inplace=False, nsamples=1):
        if(inplace):
            self.core_faust_cplx.optimize_time(transp, inplace, nsamples)
        else:
            core = FaustCoreCplx(core=True)
            core.core_faust_cplx = self.core_faust_cplx.optimize_time(transp,
                                                                      inplace,
                                                                      nsamples)
            return core

    def conj(self):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.conjugate()
        return core

    def getH(self):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.adjoint()
        return core

    def zpruneout(self, nnz_tres, npasses, only_forward):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = self.core_faust_cplx.zpruneout(nnz_tres,
                                                              npasses,
                                                              only_forward)
        return core


    def __dealloc__(self):
        del self.core_faust_cplx



