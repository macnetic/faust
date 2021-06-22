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

    @staticmethod
    def eyeFaust(n, m):
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = FaustCoreCy.FaustCoreCpp[complex].eyeFaust(n,
                                                                          m)
        return core

    @staticmethod
    def fourierFaust(n, norma):
        if(n>31):
            raise ValueError("Faust doesn't handle a FFT of order larger than "
                             "2**31")
        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = \
                FaustCoreCy.FaustCoreCpp[complex].fourierFaust(n, norma)
        if(core.core_faust_cplx == NULL):
            raise MemoryError()
        # fourier is always a complex Faust
        return core

    @staticmethod
    def polyBasis(L, K, T0=None, on_gpu=False):
        cdef int[:] colind_view, rowptr_view
        cdef complex[:] cplx_vals_view
        cdef int[:] T0_colind_view, T0_rowptr_view
        cdef complex[:] T0_cplx_vals_view
        core = FaustCoreCplx(core=True)
        if not isinstance(T0, type(None)):
            if not isinstance(T0, csr_matrix):
                raise TypeError("T0 must be a csr_matrix")
            else:
                if T0.shape[0] != L.shape[0]:
                    raise ValueError("T0 must agree.")
                if T0.dtype != L.dtype:
                    raise TypeError("T0 and L must have the same dtype.")
                T0_cplx_vals_view = T0.data
                T0_colind_view = T0.indices
                T0_rowptr_view = T0.indptr
        colind_view = L.indices
        rowptr_view = L.indptr
        cplx_vals_view = L.data
        if isinstance(T0, type(None)):
            core.core_faust_cplx = \
            FaustCoreCy.FaustCoreCpp[complex].polyBasis(L.shape[0], L.shape[1],
                                                       &rowptr_view[0],
                                                       &colind_view[0],
                                                        &cplx_vals_view[0],
                                                        L.nnz,
                                                        K,
                                                        on_gpu)
        else:
            core.core_faust_cplx = \
                    FaustCoreCy.FaustCoreCpp[complex].polyBasis_ext(L.shape[0], L.shape[1],
                                                                &rowptr_view[0],
                                                                &colind_view[0],
                                                                &cplx_vals_view[0],
                                                                L.nnz,
                                                                K,
                                                                &T0_rowptr_view[0],
                                                                &T0_colind_view[0],
                                                                &T0_cplx_vals_view[0],
                                                                T0.nnz,
                                                                T0.shape[1],
                                                                on_gpu)
        return core

    def __dealloc__(self):
        del self.core_faust_cplx

cdef class FaustFactCplx(FaustFact):

    @staticmethod
    def fact_palm4msa_fft(Lap, p):
        return FaustFactCplx.fact_palm4msa_gen(Lap, p, p.init_D)

    @staticmethod
    def fact_palm4msa(M, p):
        return FaustFactCplx.fact_palm4msa_gen(M,p)

    @staticmethod
    def fact_palm4msa_gen(M, p, init_D=None):

        # M is supposed to be already verified (np.float or np.complex, fortran
        # contiguous, 1 or 2 dimensions)
        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] Mview_cplx

        cdef complex[:,:] tmp_mat_cplx

        # views for lambda and optionally D out buffer (FGFT)
        cdef complex[:] outbufview_cplx

        # only for FGFT
        cdef complex[:] init_D_view_cplx

        cdef FaustCoreCy.PyxParamsFactPalm4MSA[complex,double]* cpp_params_cplx

        # template parameter is always double (never complex) because no need
        # a threshold error of complex type
        cdef PyxStoppingCriterion[double] cpp_stop_crit
        cdef PyxConstraintGeneric** cpp_constraints


        cpp_stop_crit.is_criterion_error = p.stop_crit._is_criterion_error
        cpp_stop_crit.error_threshold = p.stop_crit.tol
        cpp_stop_crit.num_its = p.stop_crit.num_its
        cpp_stop_crit.max_num_its = p.stop_crit.maxiter

        calling_fft_algo = isinstance(init_D, np.ndarray)

        if(not p.init_facts):
            p.init_facts = [ None for i in range(p.num_facts) ]
            if(p.is_update_way_R2L):
                zeros_id = p.num_facts-1
            else:
                zeros_id = 0
            p.init_facts[zeros_id] = \
                np.zeros([p.constraints[zeros_id]._num_rows,p.constraints[zeros_id]._num_cols],
                        order='F', dtype=M.dtype)
            for i in [i for i in range(0, p.num_facts) if i != zeros_id]:
                p.init_facts[i] = np.eye(p.constraints[i]._num_rows,
                                        p.constraints[i]._num_cols, order='F',
                                        dtype=M.dtype)

        if(calling_fft_algo):
            # FFT/FGFT case, we store lambda in first position and the diagonal
            # of D in the next
            _out_buf = np.empty(init_D.shape[0]+1, dtype=M.dtype)
        else:
            # store only lambda as a return from Palm4MSA algo
            _out_buf = np.array([0], dtype=M.dtype)

        if(calling_fft_algo):
            cpp_params_cplx = new \
            FaustCoreCy.PyxParamsFactPalm4MSAFFT[complex,double]()
            init_D_view_cplx = init_D
            (<FaustCoreCy.PyxParamsFactPalm4MSAFFT[complex,double]*>cpp_params_cplx).init_D = &init_D_view_cplx[0]
        else:
            cpp_params_cplx = new \
            FaustCoreCy.PyxParamsFactPalm4MSA[complex,double]()

        Mview_cplx=M
        cpp_params_cplx.num_facts = p.num_facts
        cpp_params_cplx.num_facts = p.num_facts
        cpp_params_cplx.is_update_way_R2L = p.is_update_way_R2L
        cpp_params_cplx.init_lambda = p.init_lambda
        cpp_params_cplx.step_size = p.step_size
        cpp_params_cplx.grad_calc_opt_mode = p.grad_calc_opt_mode
        cpp_params_cplx.norm2_max_iter = int(p.norm2_max_iter)
        cpp_params_cplx.norm2_threshold = p.norm2_threshold
        cpp_params_cplx.stop_crit = cpp_stop_crit
        cpp_params_cplx.init_facts = <complex**> \
                PyMem_Malloc(sizeof(complex*)*p.num_facts)
        cpp_params_cplx.init_fact_sizes = <unsigned long*> \
        PyMem_Malloc(sizeof(unsigned long)*2*p.num_facts)
        cpp_params_cplx.is_verbose = p.is_verbose
        cpp_params_cplx.constant_step_size = p.constant_step_size
        outbufview_cplx = _out_buf

        cpp_constraints = \
        <PyxConstraintGeneric**> \
        PyMem_Malloc(sizeof(PyxConstraintGeneric*)*len(p.constraints))

        for i in range(0,len(p.constraints)):
            cons = p.constraints[i]
            #print("FaustFact.fact_palm4MSA() cons.name =", cons.name)
            if(cons.is_int_constraint()):
                #print("FaustFact.fact_palm4MSA() Int Constraint")
                cpp_constraints[i] = <PyxConstraintInt*> PyMem_Malloc(sizeof(PyxConstraintInt))
                (<PyxConstraintInt*>cpp_constraints[i]).parameter = cons._cons_value
            elif(cons.is_real_constraint()):
                #print("FaustFact.fact_palm4MSA() Real Constraint")
                cpp_constraints[i] = <PyxConstraintScalar[double]*> \
                PyMem_Malloc(sizeof(PyxConstraintScalar[double]))
                (<PyxConstraintScalar[double]*>cpp_constraints[i]).parameter =\
                        cons._cons_value
            elif(cons.is_mat_constraint()):
                #print("FaustFact.fact_palm4MSA() Matrix Constraint")
                cpp_constraints[i] = <PyxConstraintMat[complex]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[complex]))
                tmp_mat_cplx = cons._cons_value
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter =\
                        &tmp_mat_cplx[0,0]
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz
            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols


        cpp_params_cplx.constraints = cpp_constraints
        cpp_params_cplx.num_constraints = len(p.constraints)

        for i in range(0,p.num_facts):
            check_matrix(False, p.init_facts[i], message="while checking"
                         " palm4msa init facts: ")
            tmp_mat_cplx = p.init_facts[i].astype('complex')
            cpp_params_cplx.init_facts[i] = &tmp_mat_cplx[0,0]
            cpp_params_cplx.init_fact_sizes[i*2+0] = p.init_facts[i].shape[0]
            cpp_params_cplx.init_fact_sizes[i*2+1] = p.init_facts[i].shape[1]

        core = FaustCoreCplx(core=True)

        if(calling_fft_algo):
            core.core_faust_cplx = \
                    FaustCoreCy.fact_palm4MSAFFT[complex,double](&Mview_cplx[0,0],
                                                                M_num_rows,
                                                                M_num_cols,
                                                                <FaustCoreCy.PyxParamsFactPalm4MSAFFT[complex,double]*>cpp_params_cplx,
                                                                &outbufview_cplx[0])
        else:
            core.core_faust_cplx = FaustCoreCy.fact_palm4MSA[complex,double](&Mview_cplx[0,0], M_num_rows, M_num_cols,
                                 cpp_params_cplx, &outbufview_cplx[0])

        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        PyMem_Free(cpp_params_cplx.init_facts)
        PyMem_Free(cpp_params_cplx.init_fact_sizes)


        del cpp_params_cplx

        if(calling_fft_algo):
            return core, np.real(_out_buf[0]), _out_buf[1:]
        else:
            return core, np.real(_out_buf[0])

    @staticmethod
    def fact_hierarchical_fft(U, Lap, p, init_D):
        return FaustFactCplx.fact_hierarchical_gen(U, p, init_D, Lap)

    @staticmethod
    def fact_hierarchical(M, p):
        return FaustFactCplx.fact_hierarchical_gen(M, p)

    @staticmethod
    def fact_hierarchical_gen(M, p, init_D=None, Lap=None):

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] Mview_cplx

        cdef complex[:,:] Lapview_cplx

        # view for lambda and optionally D out buffer (FGFT)
        cdef complex[:] outbufview_cplx

        # only for FGFT
        cdef complex[:] init_D_view_cplx

        cdef complex[:,:] tmp_mat_cplx

        cdef FaustCoreCy.PyxParamsHierarchicalFact[complex,double]* cpp_params_cplx
        cdef PyxStoppingCriterion[double]* cpp_stop_crits
        # template parameter is always double (never complex) because no need
        # to have a threshold error of complex type (it wouldn't make sense)
        cdef PyxConstraintGeneric** cpp_constraints

        cpp_stop_crits = <PyxStoppingCriterion[double]*>\
        PyMem_Malloc(sizeof(PyxStoppingCriterion[double])*2)

        cpp_stop_crits[0].is_criterion_error = p.stop_crits[0]._is_criterion_error
        cpp_stop_crits[0].error_threshold = p.stop_crits[0].tol
        cpp_stop_crits[0].num_its = p.stop_crits[0].num_its
        cpp_stop_crits[0].max_num_its = p.stop_crits[0].maxiter
        cpp_stop_crits[1].is_criterion_error = p.stop_crits[1]._is_criterion_error
        cpp_stop_crits[1].error_threshold = p.stop_crits[1].tol
        cpp_stop_crits[1].num_its = p.stop_crits[1].num_its
        cpp_stop_crits[1].max_num_its = p.stop_crits[1].maxiter


        calling_fft_algo = isinstance(init_D, np.ndarray)

        if(calling_fft_algo):
            # FFT/FGFT case, we store lambda in first position and the diagonal
            # of D in the next
            _out_buf = np.empty(init_D.shape[0]+1, dtype=M.dtype)
        else:
            # store only lambda as a return from Palm4MSA algo
            _out_buf = np.array([0], dtype=M.dtype)

        if(calling_fft_algo):
            cpp_params_cplx = new \
            FaustCoreCy.PyxParamsHierarchicalFactFFT[complex,double]()
            init_D_view_cplx = init_D
            Lapview_cplx = Lap
            (<FaustCoreCy.PyxParamsHierarchicalFactFFT[complex,double]*>cpp_params_cplx).init_D = &init_D_view_cplx[0]
        else:
            cpp_params_cplx = new \
            FaustCoreCy.PyxParamsHierarchicalFact[complex,double]()
        Mview_cplx=M
        cpp_params_cplx.num_facts = p.num_facts
        cpp_params_cplx.num_facts = p.num_facts
        cpp_params_cplx.is_update_way_R2L = p.is_update_way_R2L
        cpp_params_cplx.init_lambda = p.init_lambda
        cpp_params_cplx.step_size = p.step_size
        cpp_params_cplx.grad_calc_opt_mode = p.grad_calc_opt_mode
        cpp_params_cplx.norm2_max_iter = int(p.norm2_max_iter)
        cpp_params_cplx.norm2_threshold =  p.norm2_threshold
        cpp_params_cplx.stop_crits = cpp_stop_crits
        cpp_params_cplx.is_verbose = p.is_verbose
        cpp_params_cplx.is_fact_side_left = p.is_fact_side_left
        cpp_params_cplx.constant_step_size = p.constant_step_size
        outbufview_cplx = _out_buf

        cpp_constraints = \
        <PyxConstraintGeneric**> \
        PyMem_Malloc(sizeof(PyxConstraintGeneric*)*len(p.constraints))

        for i in range(0,len(p.constraints)):
            cons = p.constraints[i]
            #print("FaustFact.fact_hierarchical() cons.name =", cons.name)
            if(cons.is_int_constraint()):
                #print("FaustFact.fact_hierarchical() Int Constraint")
                cpp_constraints[i] = <PyxConstraintInt*> PyMem_Malloc(sizeof(PyxConstraintInt))
                (<PyxConstraintInt*>cpp_constraints[i]).parameter = cons._cons_value
            elif(cons.is_real_constraint()):
                #print("FaustFact.fact_hierarchical() Real Constraint")
                cpp_constraints[i] = <PyxConstraintScalar[double]*> \
                PyMem_Malloc(sizeof(PyxConstraintScalar[double]))
                (<PyxConstraintScalar[double]*>cpp_constraints[i]).parameter =\
                        cons._cons_value
            elif(cons.is_mat_constraint()):
                #print("FaustFact.fact_hierarchical() Matrix Constraint")
                cpp_constraints[i] = <PyxConstraintMat[complex]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[complex]))
                tmp_mat_cplx = cons._cons_value
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter =\
                        &tmp_mat_cplx[0,0]
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz
            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols

        cpp_params_cplx.constraints = cpp_constraints
        cpp_params_cplx.num_rows = p.data_num_rows
        cpp_params_cplx.num_cols = p.data_num_cols
        cpp_params_cplx.num_constraints = len(p.constraints)

        core = FaustCoreCplx(core=True)
        if(calling_fft_algo):
            core.core_faust_cplx = \
                    FaustCoreCy.fact_hierarchical_fft[complex,
                                                      double](&Mview_cplx[0,0],
                                                              &Lapview_cplx[0,0], M_num_rows, M_num_cols,
                                                              <FaustCoreCy.PyxParamsHierarchicalFactFFT[complex,double]*>cpp_params_cplx,
                                                              &outbufview_cplx[0])
        else:
            core.core_faust_cplx = \
                    FaustCoreCy.fact_hierarchical[complex,
                                                  double](&Mview_cplx[0,0], M_num_rows, M_num_cols,
                                                          cpp_params_cplx, &outbufview_cplx[0])
        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        PyMem_Free(cpp_stop_crits)

        del cpp_params_cplx
        if(core.core_faust_cplx == NULL): raise Exception("fact_hierarchical"
                                                          " has failed.");

        if(calling_fft_algo):
            return core, np.real(_out_buf[0]), _out_buf[1:]
        else:
            return core, np.real(_out_buf[0])

    @staticmethod
    def palm4msa2020(M, p, full_gpu=True):
        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] Mview
        cdef complex[:,:] tmp_mat
        # view for lambda
        cdef double[:] outbufview

        cdef PyxStoppingCriterion[double] cpp_stop_crit
        cdef PyxMHTPParams[double] cpp_MHTPParams
        # template parameter is always double (never complex) because no need
        # to have a threshold error of complex type (it wouldn't make sense)
        cdef PyxConstraintGeneric** cpp_constraints

        cdef FaustCoreCy.FaustCoreCpp[complex]* core_faust_dbl_init_facts


        Mview = M
        _out_buf = np.array([0], dtype=M.dtype)
        _out_buf[0] = p.init_lambda;
        outbufview = _out_buf

        cpp_stop_crit.is_criterion_error = p.stop_crit._is_criterion_error
        cpp_stop_crit.error_threshold = p.stop_crit.tol
        cpp_stop_crit.num_its = p.stop_crit.num_its
        cpp_stop_crit.max_num_its = p.stop_crit.maxiter

        # use_MHTP is either False or a MHTPParams instance
        if p.use_MHTP != False:
            mhtpp = p.use_MHTP
            cpp_MHTPParams.used = True
            cpp_MHTPParams.stop_crit.is_criterion_error = mhtpp.stop_crit._is_criterion_error
            cpp_MHTPParams.stop_crit.error_threshold = mhtpp.stop_crit.tol
            cpp_MHTPParams.stop_crit.num_its = mhtpp.stop_crit.num_its
            cpp_MHTPParams.stop_crit.max_num_its = mhtpp.stop_crit.maxiter
            cpp_MHTPParams.constant_step_size = mhtpp.constant_step_size
            cpp_MHTPParams.step_size = mhtpp.step_size
            cpp_MHTPParams.updating_lambda = mhtpp.updating_lambda
            cpp_MHTPParams.palm4msa_period = mhtpp.palm4msa_period
        else:
            cpp_MHTPParams.used = False

        cpp_constraints = \
        <PyxConstraintGeneric**> \
        PyMem_Malloc(sizeof(PyxConstraintGeneric*)*len(p.constraints))

        p.factor_format = \
        pyfaust.factparams.ParamsFact.factor_format_str2int(p.factor_format)

        for i in range(0,len(p.constraints)):
            cons = p.constraints[i]
            #print("FaustFact.fact_palm4MSA() cons.name =", cons.name)
            if cons.is_int_constraint():
                #print("FaustFact.fact_palm4MSA() Int Constraint")
                cpp_constraints[i] = <PyxConstraintInt*> PyMem_Malloc(sizeof(PyxConstraintInt))
                (<PyxConstraintInt*>cpp_constraints[i]).parameter = cons._cons_value
            elif cons.is_real_constraint():
                #print("FaustFact.fact_palm4MSA() Real Constraint")
                cpp_constraints[i] = <PyxConstraintScalar[double]*> \
                PyMem_Malloc(sizeof(PyxConstraintScalar[double]))
                (<PyxConstraintScalar[double]*>cpp_constraints[i]).parameter =\
                        cons._cons_value
            elif cons.is_mat_constraint():
                #print("FaustFact.fact_palm4MSA() Matrix Constraint")
                cpp_constraints[i] = <PyxConstraintMat[complex]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[complex]))
                tmp_mat = cons._cons_value
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter =\
                        &tmp_mat[0,0]
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz
            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols

        if p.init_facts:
            # facts have been initialized from the wrapper
            # create a Faust
            F_facts = FaustCore(p.init_facts)
            # palm4msa2020_gen in FaustFact.hpp
            # is responsible to delete the object in case the
            # algorithm runs on GPU (hence the transform objects F_facts and
            # core are not the same)


        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = \
            FaustCoreCy.palm4msa2020[complex](&Mview[0,0], M_num_rows,
                                             M_num_cols,
                                             cpp_constraints,
                                             len(p.constraints),
                                             &outbufview[0],
                                             cpp_stop_crit,
                                             p.is_update_way_R2L,
                                             p.factor_format, p.packing_RL,
                                             cpp_MHTPParams,
                                             p.norm2_max_iter,
                                             p.norm2_threshold,
                                             p.is_verbose,
                                             p.constant_step_size,
                                             p.step_size,
                                             full_gpu,
                                             <FaustCoreCy.FaustCoreCpp[complex]*>NULL
                                              if not p.init_facts else F_facts.core_faust_cplx)

        if p.init_facts:
            F_facts.core_faust_cplx = NULL
            # needed to avoid double-free (because F_facts has the same
            # TransformHelper behind as core)

        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        if(core.core_faust_cplx == NULL): raise Exception("palm4msa2020"
                                                          " has failed.");

        return core, np.real(_out_buf[0])

    @staticmethod
    def fact_givens_fgft(Lap, J, t, verbosity=0, stoppingError = 0.0,
                         errIsRel=True, order=1, enable_large_Faust=False):
        if(order == 'ascend'): order = 1
        elif(order == 'descend'): order = -1
        elif(order == 'undef'): order = 0
        else: raise ValueError('order argument must be something among'
                               '\'ascend\', \'descend\' or \'undef\'')

        cdef unsigned int Lap_num_rows=Lap.shape[0]
        cdef unsigned int Lap_num_cols=Lap.shape[1]

        cdef complex[:,:] Lap_view
        cdef double[:] D_view

        Lap_view = Lap
        D = np.empty(Lap.shape[0], dtype='double')
        D_view = D

        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = FaustCoreCy.fact_givens_fgft_cplx[complex,double](&Lap_view[0,0],
                                                                                 Lap_num_rows,
                                                                                 Lap_num_cols, J, t,
                                                                                 &D_view[0], verbosity,
                                                                                 stoppingError,
                                                                                 errIsRel,
                                                                                 int(order),
                                                                                 enable_large_Faust)

        #from scipy.sparse import spdiags
        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if core.core_faust_cplx == NULL:
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return core, D

    @staticmethod
    def fact_givens_fgft_sparse(Lap, J, t, verbosity=0, stoppingError = 0.0,
                                     errIsRel=True, order='ascend',
                                     enable_large_Faust=False):
        from scipy.sparse import spdiags
        cdef complex[:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef unsigned int Lap_num_rows=Lap.shape[0]
        cdef unsigned int Lap_num_cols=Lap.shape[1]
        cdef double[:] D_view

        if(order == 'ascend'): order = 1
        elif(order == 'descend'): order = -1
        elif(order == 'undef'): order = 0

        data1d=Lap.data.astype(complex,'F')
        indices=Lap.indices.astype(np.int32, 'F')
        indptr=Lap.indptr.astype(np.int32, 'F')

        D = np.empty(Lap.shape[0], dtype='double')
        D_view = D

        core = FaustCoreCplx(core=True)
        core.core_faust_cplx = \
        FaustCoreCy.fact_givens_fgft_sparse_cplx[complex, double](&data1d[0],
                                                                  &indices[0],
                                                                  &indptr[0],
                                                                  Lap.nnz,
                                                                  Lap_num_rows,
                                                                  Lap_num_cols, J, t,
                                                                  &D_view[0],
                                                                  verbosity,
                                                                  stoppingError,
                                                                  errIsRel,
                                                                  int(order),
                                                                  enable_large_Faust)

        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if core.core_faust_cplx == NULL:
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return core, D

    @staticmethod
    def eigtj(M, nGivens=None, tol=0, relerr=True,  nGivens_per_fac=None, verbosity=0,
          order='ascend', enable_large_Faust=False):
        if nGivens == None:
            if tol == 0:
                raise Exception("You must specify nGivens or tol argument"
                        " (to define a stopping  criterion)")
            nGivens = 0
        if nGivens_per_fac == None: nGivens_per_fac = int(M.shape[0]/2)
        if isinstance(M, np.ndarray) and \
            not np.allclose(np.matrix(M, copy=False).H, M) or M.shape[0] != M.shape[1]:
            raise ValueError(" the matrix/array must be symmetric or hermitian.")
        if not isinstance(nGivens, int): raise TypeError("nGivens must be a int")
        if not isinstance(nGivens_per_fac, int): raise TypeError("nGivens_per_fac must be a int")
        nGivens_per_fac = max(nGivens_per_fac, 1)
        if nGivens > 0: nGivens_per_fac = min(nGivens_per_fac, nGivens)
        tol *= tol # the C++ impl. works on squared norms to measure errors

        if isinstance(M, np.ndarray):
            core_obj,D = FaustFactCplx.fact_givens_fgft(M, nGivens, nGivens_per_fac,
                    verbosity, tol,
                    relerr, order, enable_large_Faust)
        elif isinstance(M, csr_matrix):
            core_obj,D = FaustFactCplx.fact_givens_fgft_sparse(M, nGivens, nGivens_per_fac,
                    verbosity, tol,
                    relerr, order, enable_large_Faust)
        else:
            raise TypeError("The matrix to diagonalize must be a"
                            " scipy.sparse.csr_matrix or a numpy array.")
        return D, core_obj

    @staticmethod
    def svdtj_cplx(M, J, t, verbosity=0, stoppingError = 0.0,
              errIsRel=True, enable_large_Faust=False):

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] M_view
        cdef complex[:] S_view

        M_view = M
        S = np.empty(M.shape[0], dtype=M.dtype)
        S_view = S
        stoppingError *= stoppingError # the C++ impl. works on squared norms to measure errors

        coreU = FaustCoreCplx(core=True)
        coreV = FaustCoreCplx(core=True)
        FaustCoreCy.svdtj_cplx[complex,double](&(coreU.core_faust_cplx),
                                          &(coreV.core_faust_cplx),
                                          &S_view[0],
                                          &M_view[0,0],
                                          M_num_rows,
                                          M_num_cols, int(J), int(t),
                                          verbosity,
                                          stoppingError,
                                          errIsRel, enable_large_Faust)

        #from scipy.sparse import spdiags
        #S_spdiag = spdiags(S, [0], M.shape[0], M.shape[0])
        #return core, S_spdiag
        if coreU.core_faust_cplx == NULL:
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")
        return coreU, S, coreV

    @staticmethod
    def svdtj_sparse_cplx(M, J, t, verbosity=0, stoppingError = 0.0,
                     errIsRel=True, enable_large_Faust=False):
        from scipy.sparse import spdiags
        cdef complex [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]
        cdef complex[:] S_view


        data1d = M.data.astype(np.complex,'F')
        indices = M.indices.astype(np.int32, 'F')
        indptr = M.indptr.astype(np.int32, 'F')

        S = np.empty(M.shape[0], dtype=M.dtype)
        S_view = S
        stoppingError *= stoppingError # the C++ impl. works on squared norms to measure errors

        coreU = FaustCoreCplx(core=True)
        coreV = FaustCoreCplx(core=True)
        FaustCoreCy.svdtj_sparse_cplx[complex,double](
            &(coreU.core_faust_cplx),
            &(coreV.core_faust_cplx),
            &S_view[0],
            &data1d[0],
            &indices[0],
            &indptr[0],
            M.nnz,
            M_num_rows,
            M_num_cols, J, t,
            verbosity,
            stoppingError,
            errIsRel, enable_large_Faust)

        coreU._isReal = coreV._isReal = False
        #D_spdiag = spdiags(D, [0], M.shape[0], M.shape[0])
        #return core, D_spdiag

        if coreU.core_faust_cplx == NULL:
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")
        return coreU, S, coreV

    @staticmethod
    def butterfly_hierarchical(M, dir):

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] Mview_cplx

        if dir == "right":
            dir = 1
        elif dir == "left":
            dir = 0
        else:
            raise ValueError("dir argument must be 'right' or 'left'.")



        core = FaustCoreCplx(core=True)
        Mview_cplx = M
        core.core_faust_cplx = \
                FaustCoreCy.butterfly_hierarchical[complex](&Mview_cplx[0,0],
                                                            M_num_rows, M_num_cols,
                                                            dir)
        return core

def polyCoeffsCplx(d, basisX, coeffs, dev, out=None):
    K = coeffs.size-1
    # basisX must be a numpy array in fortran column major order and dtype complex128
    if basisX.ndim > 1:
        n = basisX.shape[1]
    else:
        n = 1

    ndim_M=basisX.ndim

    if (ndim_M > 2) | (ndim_M < 1):
        raise ValueError('input basisX invalid number of dimensions')

    cdef unsigned int nbrow_x=basisX.shape[0]
    cdef unsigned int nbcol_x #can't be assigned because we don't know yet if the input vector is 1D or 2D

    cdef unsigned int nbrow_y=d
    cdef unsigned int nbcol_y=n

    cdef complex[:] bxview_1D_cplx
    cdef complex[:,:] bxview_2D_cplx

    cdef complex[:] cview_cplx

    cdef y
    cdef complex[:,:] yview_cplx

    if ndim_M == 1:
        nbcol_x=1
        bxview_1D_cplx = basisX
    else:
        nbcol_x=basisX.shape[1]
        bxview_2D_cplx=basisX

    # TODO: make some sanity checks on argument sizes

    cview_cplx = coeffs

    if not isinstance(out, type(None)):
        y = out
        if y.ndim == 1:
            raise ValueError('out must have 2 dimensions.')
        if y.shape != (d,n):
            raise ValueError('out shape isn\'t valid.')
        dtype_err = ValueError('out dtype isn\'t valid.')
        if not y.flags['F_CONTIGUOUS']:
            raise ValueError('the array must be in fortran/column continous order.')
        if y.dtype != 'complex':
            raise dtype_err
        yview_cplx = y
    else:
        y = np.empty((nbrow_y, nbcol_y), dtype='complex', order='F')
        yview_cplx = y

    #void polyCoeffs(int d, int K, int n, const FPP* basisX, const FPP* coeffs, FPP* out) const;
    if ndim_M == 1:
        FaustCoreCy.polyCoeffs(d, K, n, &bxview_1D_cplx[0],
                   &cview_cplx[0], &yview_cplx[0,0], dev.startswith('gpu'))
        y = np.squeeze(y) # we want a single dim. (but we created two
        # above)
    else:
        FaustCoreCy.polyCoeffs(d, K, n, &bxview_2D_cplx[0,0],
                   &cview_cplx[0], &yview_cplx[0,0], dev.startswith('gpu'))
    return y

def polyGroupCoeffsCplx(d, K, basisX, coeffs, dev, out=None):
    # basisX must be a numpy array in fortran column major order and dtype complex128
    if basisX.ndim > 1:
        n = basisX.shape[1]
    else:
        n = 1

    ndim_M=basisX.ndim

    if (ndim_M > 2) | (ndim_M < 1):
        raise ValueError('input basisX invalid number of dimensions')

    cdef unsigned int nbrow_x=basisX.shape[0]
    cdef unsigned int nbcol_x #can't be assigned because we don't know yet if the input vector is 1D or 2D

    cdef unsigned int nbrow_y=d
    cdef unsigned int nbcol_y=n

    cdef complex[:] bxview_1D_cplx
    cdef complex[:,:] bxview_2D_cplx

    cdef complex[:, :] cview_cplx

    cdef y
    cdef size_t addr
    cdef complex** yview_cplx

    if ndim_M == 1:
        nbcol_x=1
        bxview_1D_cplx = basisX
    else:
        nbcol_x=basisX.shape[1]
        bxview_2D_cplx=basisX

    # TODO: make some sanity checks on argument sizes

    if coeffs.ndim != 2:
        raise ValueError('coeffs must have 2 dimensions.')

    if coeffs.shape[1] != K+1:
        raise ValueError('coeffs.shape[1] must be equal to K+1.')

    cview_cplx = coeffs

    ncoeffs = coeffs.shape[0]

    yview_cplx = <complex**> PyMem_Malloc(sizeof(complex**) * ncoeffs)

    if not isinstance(out, type(None)):
        if not isinstance(out, list):
            raise TypeError('out must be a list')
        if len(out) != ncoeffs:
            raise ValueError('out length must agree with coeffs.shape[0]')
        for i in range(ncoeffs):
            y = out[i]
            if y.ndim == 1:
                raise ValueError('out must have 2 dimensions.')
            if y.shape != (d,n):
                raise ValueError('out shape isn\'t valid.')
            dtype_err = ValueError('out dtype isn\'t valid.')
            if not y.flags['F_CONTIGUOUS']:
                raise ValueError('the array must be in fortran/column continous order.')
            if y.dtype != 'complex':
                raise dtype_err
            addr = y.__array_interface__['data'][0]
            yview_cplx[i] = <complex*> addr
    else:
        out = []
        for i in range(ncoeffs):
            out.append(np.empty((nbrow_y, nbcol_y), dtype='complex',
                                order='F'))
            addr = out[i].__array_interface__['data'][0]
            yview_cplx[i] = <complex*> addr

#void polyCoeffsSeq(int d, uint K, int n, const FPP* basisX, const FPP* coeffs, FPP** out, int n_out, bool on_gpu);
    if ndim_M == 1:
        FaustCoreCy.polyGroupCoeffs_(d, K, n, &bxview_1D_cplx[0],
                   &cview_cplx[0,0], yview_cplx, ncoeffs, dev.startswith('gpu'))
    else:
        FaustCoreCy.polyGroupCoeffs_(d, K, n, &bxview_2D_cplx[0,0],
                   &cview_cplx[0,0], yview_cplx, ncoeffs, dev.startswith('gpu'))

    PyMem_Free(yview_cplx)

    return out

cdef class ConstraintIntCoreCplx:

    # no need to create an object, because project() needs to create core object
    # for each call and for that purpose it needs to determine if mat is real or complex
    # so it can't create the object before the call
    # in that conditions a static method will suffice
    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, normalized=True,
                pos=False):
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')


        M_view_cplx = M
        M_out_view_cplx = M_out
        ret = FaustCoreCy.prox_int[complex](name, parameter, &M_view_cplx[0,0],
                                      num_rows, num_cols,
                                      &M_out_view_cplx[0,0], normalized,
                                      pos)

        if ret == -1:
            raise ZeroDivisionError("Can't normalize because norm is zero.");
        return M_out

cdef class ConstraintMatCoreCplx:

    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, parameter_sz, normalized=False,
                pos=False):
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx
        cdef complex[:,:] param_view_cplx

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')


        parameter = parameter.astype(M.dtype)

        M_view_cplx = M
        M_out_view_cplx = M_out
        param_view_cplx = parameter
        ret = FaustCoreCy.prox_mat[complex](name, &param_view_cplx[0,0],
                                      parameter_sz, &M_view_cplx[0,0],
                                      num_rows, num_cols,
                                      &M_out_view_cplx[0,0], normalized,
                                      pos)

        if ret == -1:
            raise ZeroDivisionError("Can't normalize because norm is zero.");
        return M_out

    @staticmethod
    def prox_blockdiag(M, block_shapes, normalized, pos):
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx
        cdef unsigned long int* m_ptr
        cdef unsigned long int* n_ptr

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')

        n_ptr = <unsigned long*>PyMem_Malloc(sizeof(unsigned long *)*len(block_shapes))
        m_ptr = <unsigned long*>PyMem_Malloc(sizeof(unsigned long *)*len(block_shapes))

        m_ptr[0] = block_shapes[0][0]
        n_ptr[0] = block_shapes[0][1]
        for i in range(1, len(block_shapes)):
            m_ptr[i] = block_shapes[i][0]+m_ptr[i-1]
            n_ptr[i] = block_shapes[i][1]+n_ptr[i-1]

        if(m_ptr[len(block_shapes)-1] != M.shape[0] or n_ptr[len(block_shapes)-1] != M.shape[1]):
            raise ValueError("The sum of block shapes is not equal to the matrix shapes.")

        M_view_cplx = M
        M_out_view_cplx = M_out
        ret = FaustCoreCy.prox_blockdiag[complex](&M_view_cplx[0,0], M.shape[0],
                                           M.shape[1], &m_ptr[0], &n_ptr[0],
                                          len(block_shapes),normalized, pos,
                                           &M_out_view_cplx[0,0])
        PyMem_Free(m_ptr)
        PyMem_Free(n_ptr)

        if ret == -1:
            raise ZeroDivisionError("Can't normalize because norm is zero.");

        return M_out

cdef class ConstraintRealCoreCplx:

    # no need to create an object, because project() needs to create core object
    # for each call and for that purpose it needs to determine if mat is real or complex
    # so it can't create the object before the call
    # in that conditions a static method will suffice
    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, normalized=False,
                pos=False):
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')

        M_view_cplx = M
        M_out_view_cplx = M_out
        ret = FaustCoreCy.prox_real[complex, double](name, parameter, &M_view_cplx[0,0],
                                      num_rows, num_cols,
                                               &M_out_view_cplx[0,0],
                                               normalized, pos)

        if ret == -1:
            raise ZeroDivisionError("Can't normalize because norm is zero.");

        return M_out
