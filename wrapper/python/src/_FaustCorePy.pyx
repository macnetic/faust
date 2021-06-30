##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython source file making the links between :                     ##
##        the Cython module FaustCoreCy and Python                          ##
##                                                                          ## 
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##                                                                          ##
##############################################################################

# importer le module Cython faisant le lien avec la class C++
cimport FaustCoreCy
from FaustCoreCy cimport PyxConstraintGeneric, PyxConstraintInt, \
PyxConstraintMat, PyxConstraintScalar, PyxStoppingCriterion, \
PyxParamsFactPalm4MSA, PyxMHTPParams

import numpy as np
cimport numpy as np
import copy

from libc.stdlib cimport malloc, free;
from libc.stdio cimport printf
from libc.string cimport memcpy, strlen;
from libcpp cimport bool, complex
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from scipy import sparse
from scipy.sparse import csr_matrix, csc_matrix
from re import match
import sys, os, pyfaust

cdef class FaustCore:

    #### ATTRIBUTE ########
    # classe Cython
    cdef FaustCoreCy.FaustCoreCpp[double]* core_faust_dbl
    cdef FaustCoreCy.FaustCoreCpp[complex]* core_faust_cplx
    #### CONSTRUCTOR ####
    #def __cinit__(self,np.ndarray[double, mode="fortran", ndim=2] mat):
    def  __cinit__(self,list_factors=None, alpha=1.0, core=False,
                   optimizedCopy=False):
        cdef double [:,:] data
        cdef double [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
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
            self.core_faust_dbl = new FaustCoreCy.FaustCoreCpp[double]()
            for factor in list_factors:
                nbrow=factor.shape[0];
                nbcol=factor.shape[1];
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
        elif(core): # trick to initialize a new FaustCoreCpp from C++ (see
        # transpose, conj and adjoint)
            pass
        #else:
        #TODO: raise error for undefined object here

    def clone(self, *args, **kwargs):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.clone()
        return core

    @staticmethod
    def _is_gpu_mod_enabled():
        return FaustCoreCy._is_gpu_mod_enabled()

    @staticmethod
    def enable_gpu_mod(libpaths=None, backend='cuda', silent=False, fatal=False):
        #TODO: extract out of Faust class
        cdef char * c_libpath
        cdef void* gm_handle
        if libpaths == None:
            # use default paths
            pyfaust_path = os.path.dirname(pyfaust.__file__)
            # don't use os.path.join or os.path.sep because anyway
            # lib name suffix or prefix depend on OS
            if sys.platform == 'linux':
                libpaths = [pyfaust_path+"/lib/libgm.so"]
            elif sys.platform == 'darwin':
                libpaths = [pyfaust_path+"/lib/libgm.dylib"]
            elif sys.platform == 'win32':
                libpaths = [pyfaust_path+"\lib\gm.dll"]
        for libpath in libpaths:
            blibpath = libpath.encode('utf-8')
            c_libpath = blibpath
            #c_libpath = libpath
            gm_handle = FaustCoreCy._enable_gpu_mod(c_libpath, silent);
            if gm_handle == NULL:
                if fatal and libpaths[-1] == libpath:
                    raise Exception("Can't load gpu_mod library, maybe the path ("
                                    +libpath+") is not"
                                    " correct or the backend (cuda) is not installed or"
                                    " configured properly so"
                                    " the libraries are not found.")

    @staticmethod
    def randFaust(faust_nrows, faust_ncols, t, field, min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density=0.1, per_row=True):
        if(field == 3): # real
            core = FaustCore(core=True)
            core.core_faust_dbl = \
            FaustCoreCy.FaustCoreCpp[double].randFaust(faust_nrows,
                                                       faust_ncols,
                                                       t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density, per_row)
            if(core.core_faust_dbl == NULL): raise MemoryError()
        elif(field == 4): # complex
            core = @FAUST_CORE_CPLX_TYPE@(core=True)
            core.core_faust_cplx = \
            FaustCoreCy.FaustCoreCpp[complex].randFaust(faust_nrows,
                                                        faust_ncols,
                                                        t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density, per_row)
            if(core.core_faust_cplx == NULL): raise MemoryError()
        else:
            raise ValueError("FaustCorePy.randFaust(): field must be 3 for real or"
                             " 4 for complex")
        return core

    def shape(self):
        cdef unsigned int nbrow = 0
        cdef unsigned int nbcol = 0
        nbrow = self.core_faust_dbl.getNbRow();
        nbcol = self.core_faust_dbl.getNbCol();
        return (nbrow,nbcol)

    def nbytes(self):
        nbytes = self.core_faust_dbl.getNBytes();
        return nbytes

    def multiply_csr_mat(self, X):
        cdef double [:] x_data1d
        cdef int [:] x_indices
        cdef int [:] x_indptr
        cdef double [:,:] y_data
        # X is supposed to be a csr_matrix
        x_indices = X.indices
        x_indptr = X.indptr
        x_nnz = X.nnz
        nbcol = X.shape[1]
        e = Exception("Dimensions must agree")
        x_data1d = X.data
        nbrow = self.core_faust_dbl.getNbRow()
        if(self.core_faust_dbl.getNbCol() != X.shape[0]): raise e
        y_data_arr = np.empty((nbrow,nbcol), dtype=np.double, order='F') # we don't know beforehand Y nnz
        y_data = y_data_arr
        self.core_faust_dbl.multiply(&y_data[0,0], nbrow, nbcol,
                                     &x_data1d[0], &x_indptr[0],
                                     &x_indices[0],
                                     x_nnz, X.shape[0], X.shape[1])
        return y_data_arr

    def get_product(self):
        cdef double [:,:] y_data

        y_arr = np.empty((self.core_faust_dbl.getNbRow(), self.core_faust_dbl.getNbCol()), dtype='d',
                         order='F')
        y_data = y_arr
        self.core_faust_dbl.get_product(&y_data[0,0], y_arr.shape[0],
                                        y_arr.shape[1])
        return y_arr

    def multiply_faust(self, F):
        if(isinstance(F, FaustCore)):
            if(F.isReal()):
                core = FaustCore(core=True)
                core.core_faust_dbl = \
                        self.core_faust_dbl.mul_faust((<FaustCore?>F).core_faust_dbl)
            else:
                core = @FAUST_CORE_CPLX_TYPE@(core=True)
                core = \
               <@FAUST_CORE_CPLX_TYPE@?> ((<FaustCore?> self)._ascomplex()).multiply_faust(F)
            return core
        raise ValueError("F must be a Faust object")

    cdef _ascomplex(self, scalar=1.0):
        cplx_facs = [self.get_fact_opt(i).astype(np.complex) for i in \
                              range(0,self.get_nb_factors())]
        cplx_facs[0] *= scalar # instead of passing the scal to the
        # construc. It avoids disp of
        # deprecation warning
        core = @FAUST_CORE_CPLX_TYPE@(cplx_facs)#, alpha=scalar)
        return core

    def device(self):
        # always returns cpu but calls the cpp code just in case
        cdef char c_str[256]
        self.core_faust_dbl.device(c_str)
        cdef length = strlen(c_str)
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        return py_str

    def multiply_scal(self, scalar):
        core = FaustCore(core=True)
        if(isinstance(scalar, int)):
            scalar = float(scalar)
        scalar_type_err = TypeError("The mul. scalar must be a real or a"
                                    " complex number")
        if(isinstance(scalar, float)):
            core.core_faust_dbl = \
                    self.core_faust_dbl.mul_scal(scalar)
        elif(isinstance(scalar, np.complex)):
            core = <@FAUST_CORE_CPLX_TYPE@?>self._ascomplex(scalar)
        else:
            raise scalar_type_err
        return core

    cdef _vertcat(self, F):
         if(not F.isReal()):
             core = @FAUST_CORE_CPLX_TYPE@(core=True)
             self = self._ascomplex()
             core.core_faust_cplx = self.core_faust_cplx.vertcat((<FaustCore?>F).core_faust_cplx)
         else:
             core = FaustCore(core=True)
             core.core_faust_dbl = self.core_faust_dbl.vertcat((<FaustCore?>F).core_faust_dbl)

         return core

    def vertcat(self,F):
        # TODO: merge with _vertcat 
        return self._vertcat(F)


    cdef _horzcat(self, F):
         #TODO: refactor with _vertcat(), maybe by passing func ref to a cat func
         #TODO/ F must be a FaustCore
         if(not F.isReal()):
             core = @FAUST_CORE_CPLX_TYPE@(core=True)
             self = self._ascomplex()
             core.core_faust_cplx = self.core_faust_cplx.horzcat((<FaustCore?>F).core_faust_cplx)
         else:
             core = FaustCore(core=True)
             core.core_faust_dbl = self.core_faust_dbl.horzcat((<FaustCore?>F).core_faust_dbl)
        #TODO:PyMem_Free
         return core

    def horzcat(self,F):
        # TODO: merge with _horzcat
        return self._horzcat(F)

    cdef _vertcatn(self, Fs):
        cdef FaustCoreCy.FaustCoreCpp[double]** _Fs

        _Fs = <FaustCoreCy.FaustCoreCpp[double]**> PyMem_Malloc(sizeof(void*) *
                                                                len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                _Fs[i] = (<FaustCore?>F).core_faust_dbl

        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.vertcatn(_Fs, len(Fs))
        #TODO:PyMem_Free
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
        cdef FaustCoreCy.FaustCoreCpp[double]** _Fs

        _Fs = <FaustCoreCy.FaustCoreCpp[double]**> PyMem_Malloc(sizeof(void*) *
                                                                len(Fs))
        for i, F in enumerate(Fs):
            if F.isReal():
                _Fs[i] = (<FaustCore?>F).core_faust_dbl
            else:
                # can't happen if vertcatn is the caller
                raise TypeError("cat Faust with different dtype")

        core = FaustCore(core=True)
        # print("self.isReal():", self.isReal())
        core.core_faust_dbl = self.core_faust_dbl.horzcatn(_Fs, len(Fs))

        #TODO:PyMem_Free
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
        core = FaustCore(core=True)
        cdef double[:] cview
        cview = coeffs
        core.core_faust_dbl = self.core_faust_dbl.polyCoeffs(&cview[0])
        return core

    def mulPolyCoeffs(self, coeffs, X, out=None):
        cdef double[:] cview
        cdef double[:,:] Xview
        cdef double[:,:] Yview

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
            if Y.dtype != 'd':
                raise dtype_err
            Yview = Y
            Xview = X
            cview = coeffs
        else:
            cview = coeffs
            Y = np.empty((X.shape[0], n), dtype='double', order='F')
            Yview = Y
            Xview = X

        self.core_faust_dbl.mulPolyCoeffs(&Xview[0,0],
                                          n,
                                          &Yview[0,0],
                                          &cview[0])
        if X_1dim:
            # X is a vector, Y must be one too
            Y = np.squeeze(Y)

        return Y

    def polyNext(self):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.polyNext()
        return core

    # Left-Multiplication by a Faust F
    # y=multiply(F,M) is equivalent to y=F*M
    def multiply(self,M):
        if isinstance(M, FaustCore):
            return self.multiply_faust(M)
        if not isinstance(M, np.ndarray):
            raise ValueError('input M must be a numpy.ndarray')
        M=M.astype(float,'F')
        if not M.dtype=='float':
           raise ValueError('input M must be a double array')
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

        cdef double[:] xview_1D
        cdef double[:,:] xview_2D

        if ndim_M == 1:
            nbcol_x=1
            xview_1D=M
        else:
            nbcol_x=M.shape[1]
            xview_2D=M

            if (nbrow_x != nbColThis):
                raise ValueError('y=F*M multiplication with Faust: invalid dimension of the input matrix M');

        #void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose);
        nbcol_y = nbcol_x;

        cdef y
        cdef double[:,:] yview
        y = np.empty([nbrow_y,nbcol_y], dtype='d',order='F')
        yview = y

        if ndim_M == 1:
            self.core_faust_dbl.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_1D[0],nbrow_x,nbcol_x)
            y = np.squeeze(y) # we want a single dim. (but we created two
            # above)
        else:
            self.core_faust_dbl.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_2D[0,0],nbrow_x,nbcol_x)

        return y

    def set_FM_mul_mode(self, mode):
        self.core_faust_dbl.set_FM_mul_mode(mode)

    def set_Fv_mul_mode(self, mode):
        self.core_faust_dbl.set_Fv_mul_mode(mode)


    # print information about the faust (size, number of factor, type of factor (dense/sparse) ...)
    def display(self):
        self.core_faust_dbl.Display()

    def to_string(self):
        cdef const char* c_str
        c_str = self.core_faust_dbl.to_string()
        cdef length = strlen(c_str)
        #printf("%s", c_str[:length])
        #py_str = str(c_str[:length], 'UTF-8')
        py_str = c_str[:length].decode('UTF-8', 'ignore')
        free(<void*>c_str)
        return py_str

    def nnz(self):
        cdef unsigned long long nnz = 0
        nnz = self.core_faust_dbl.nnz()
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
            norm = self.core_faust_dbl.norm(ord, threshold, max_num_its)
        elif(ord == np.inf):
            norm = self.core_faust_dbl.normInf()
        else:
            norm = self.core_faust_dbl.normFro()
        return norm


    def power_iteration(self, threshold, max_num_its):
        cdef double[:] _lambda_dbl_view
        _lambda = np.empty((1,), dtype=np.double)
        _lambda_dbl_view = _lambda
        self.core_faust_dbl.power_iteration(&_lambda_dbl_view[0], threshold, max_num_its)
        return _lambda[0]

    def normalize(self, ord):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.normalize(ord)
        return core

    def get_nb_factors(self):
        cdef int nb_factors
        nb_factors = int(self.core_faust_dbl.get_nb_factors())
        return nb_factors

    def get_fact(self,i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef double[:,:] fact_dbl_view
        fact = np.empty([self.core_faust_dbl.get_fact_nb_rows(i),
                         self.core_faust_dbl.get_fact_nb_cols(i)], dtype='d',
                        order='F')
        fact_dbl_view = fact
        self.core_faust_dbl.get_fact(i, &fact_dbl_view[0, 0])
        return fact

    def get_fact_opt(self, i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef double[:,:] fact_dbl_view
        cdef rowptr, col_ids, elts, nnz
        cdef int[:] rowptr_view, col_ids_view
        dtype = 'd'
        is_fact_sparse = self.core_faust_dbl.is_fact_sparse(i)
        is_transposed = self.core_faust_dbl.isTransposed()
        if(is_fact_sparse):
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
            fact = np.ndarray([self.core_faust_dbl.get_fact_nb_rows(i),
                         self.core_faust_dbl.get_fact_nb_cols(i)], dtype=dtype,
                        order=order)
            fact_dbl_view = fact
            self.core_faust_dbl.get_fact_dense(i, &fact_dbl_view[0, 0],
                                           <unsigned int*>NULL,
                                           <unsigned int*>NULL,
                                           is_transposed)
        return fact

    def slice(self, indices):
        # TODO: rename this function or cut in two: slice and fancy indexing
        core = FaustCore(core=True)
        start_row_id, end_row_id, start_col_id, end_col_id = (indices[0].start,
                                                              indices[0].stop,
                                                              indices[1].start,
                                                              indices[1].stop)
        core.core_faust_dbl = self.core_faust_dbl.slice(start_row_id,
                                                        end_row_id,
                                                        start_col_id,
                                                        end_col_id)
        return core

    def left(self, id):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.left(id)
        return core

    def right(self, id):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.right(id)
        return core

    def fancy_idx(self, indices):
        cdef unsigned long int[:] row_indices_view
        cdef unsigned long int[:] col_indices_view
        core = FaustCore(core=True)
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
        core.core_faust_dbl = \
        self.core_faust_dbl.fancy_idx(&row_indices_view[0], row_indices.size,
                                      &col_indices_view[0], col_indices.size)
        return core


    def save_mat_file(self,filepath):
        cdef char * cfilepath = <char*> PyMem_Malloc(sizeof(char) *
                                                     (len(filepath)+1))
        fparr = bytearray(filepath, "UTF-8");
        for i in range(0,len(filepath)):
            cfilepath[i] = fparr[i]
        cfilepath[i+1] = 0
        ret = self.core_faust_dbl.save_mat_file(cfilepath)
        if(not ret):
            raise Exception("Failed to save the file: "+filepath)
        PyMem_Free(cfilepath)

    def transpose(self):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.transpose()
        return core

    def swap_cols(self, id1, id2, permutation, inplace):
        if(inplace):
            self.core_faust_dbl.swap_cols(id1, id2,
                                          permutation, inplace)
            return self
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.swap_cols(id1, id2,
                                                            permutation, inplace)
        return core

    def swap_rows(self, id1, id2, permutation, inplace):
        if(inplace):
            self.core_faust_dbl.swap_rows(id1, id2,
                                          permutation, inplace)
            return self
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.swap_rows(id1, id2,
                                                            permutation, inplace)
        return core

    def optimize_storage(self, time=False):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.optimize_storage(time)
        return core

    def optimize(self, transp=False):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.optimize(transp)
        return core

    def optimize_time(self, transp=False, inplace=False, nsamples=1):
        if(inplace):
            self.core_faust_dbl.optimize_time(transp, inplace, nsamples)
        else:
            core = FaustCore(core=True)
            core.core_faust_dbl = self.core_faust_dbl.optimize_time(transp,
                                                                  inplace,
                                                                  nsamples)
            return core

    def conj(self):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.conjugate()
        return core

    def getH(self):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.adjoint()
        return core

    def zpruneout(self, nnz_tres, npasses, only_forward):
        core = FaustCore(core=True)
        core.core_faust_dbl = self.core_faust_dbl.zpruneout(nnz_tres,
                                                            npasses,
                                                            only_forward)
        return core

    def isReal(self):
        return True

    def __dealloc__(self):
        del self.core_faust_dbl

cdef check_matrix(isReal, M, message=""):
        if not isinstance(M, (np.ndarray) ):
            raise ValueError(message+'input must be a numpy ndarray')
        if(isReal):
#            M=M.astype(float,'F')
            if not M.dtype=='float':
                raise ValueError(message+'input numpy array dtype must be double (not'
                                 ' float)')
        else:
#            M=M.astype(complex,'F')
            if(M.dtype not in ['complex', 'complex128'] ): #could fail if complex128 etc.
                raise ValueError('input array must be complex(128) array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError(message+'input array must be Fortran contiguous (Colmajor)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError(message+'input matrix/array invalid number of dimensions')

import sys, os, pyfaust
# tries to load the libgm library silently,
# if not enabled at build time it will do nothing
FaustCore.enable_gpu_mod(silent=True)

cdef class FaustCoreCplxDummy(FaustCore):
    # this type is used only if complex support is not enabled
    def __cinit__(self, *args, **kwargs):
        pass
