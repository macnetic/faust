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
PyxConstraintMat, PyxConstraintScalar, PyxStoppingCriterion, PyxParamsFactPalm4MSA

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
                self.core_faust_dbl = new FaustCoreCy.FaustCoreCpp[double]()
            else:
                self.core_faust_cplx = new FaustCoreCy.FaustCoreCpp[complex]()
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

    @staticmethod
    def enable_gpu_mod(libpaths=None, backend='cuda', silent=False):
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
            if(gm_handle == NULL and not silent and libpaths[-1] == libpath):
                raise Exception("Can't load gpu_mod library, maybe the path ("
                                +libpath+") is not"
                                " correct or the backend (cuda) is not installed or"
                                " configured properly so"
                                " the libraries are not found.")

    @staticmethod
    def randFaust(t,field,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density=0.1, per_row=True):
        core = FaustCore(core=True)
        if(field == 3):
            core.core_faust_dbl = FaustCoreCy.FaustCoreCpp[double].randFaust(t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density, per_row)
            core._isReal = True
            if(core.core_faust_dbl == NULL): raise MemoryError()
        elif(field == 4):
            core.core_faust_cplx = FaustCoreCy.FaustCoreCpp[complex].randFaust(t,min_num_factors, max_num_factors, min_dim_size,
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
        core = FaustCore(core=True)
        core.core_faust_dbl = FaustCoreCy.FaustCoreCpp[double].hadamardFaust(n,
                                                                            norma)
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
        core = FaustCore(core=True)
        core.core_faust_cplx = \
                FaustCoreCy.FaustCoreCpp[complex].fourierFaust(n, norma)
        if(core.core_faust_cplx == NULL):
            raise MemoryError()
        # fourier is always a complex Faust
        core._isReal = False
        return core

    @staticmethod
    def eyeFaust(n, m, t='real'):
        core = FaustCore(core=True)
        if(t == 'real'):
            core.core_faust_dbl = FaustCoreCy.FaustCoreCpp[double].eyeFaust(n,
                                                                             m)
            core._isReal = True
        elif(t == 'complex'):
            core.core_faust_cplx = FaustCoreCy.FaustCoreCpp[complex].eyeFaust(n,
                                                                             m)
            core._isReal = False
        return core

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

    def nbytes(self):
        if(self._isReal):
            nbytes = self.core_faust_dbl.getNBytes();
        else:
            nbytes = self.core_faust_cplx.getNBytes();
        return nbytes

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
        if(X.dtype in [ 'float', 'float128',
                       'float16', 'float32',
                       'float64', 'double']):
            x_data1d = X.data
            nbrow = self.core_faust_dbl.getNbRow()
            if(self.core_faust_dbl.getNbCol() != X.shape[0]): raise e
            y_data_arr = np.empty((nbrow,nbcol), dtype=np.double, order='F') # we don't know beforehand Y nnz
            y_data = y_data_arr
            if(self._isReal):
                self.core_faust_dbl.multiply(&y_data[0,0], nbrow, nbcol,
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
                self.core_faust_cplx.multiply(&y_data_cplx[0,0], nbrow, nbcol,
                                          &x_data_cplx[0], &x_indptr[0],
                                          &x_indices[0],
                                          x_nnz, X.shape[0], X.shape[1])
            else:
                raise("x must be complex if Faust is")
                # shouldn't happen normally (avoided by calling function)
        return y_data_arr

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
        if(isinstance(F, FaustCore)):
            core = FaustCore(core=True)
            core._isReal = self._isReal
            if(self._isReal):
                if(F.isReal()):
                    core.core_faust_dbl = \
                            self.core_faust_dbl.mul_faust((<FaustCore?>F).core_faust_dbl)
                else:
                    core = \
                   <FaustCore?> (<FaustCore?> self._ascomplex()).multiply_faust(F)
            else:
                if(F.isReal()):
                    core.core_faust_cplx = \
                             self.core_faust_cplx.mul_faust((<FaustCore?>((<FaustCore?>((<FaustCore?>F)._ascomplex())))).core_faust_cplx)
                else:
                    core.core_faust_cplx = \
                            self.core_faust_cplx.mul_faust((<FaustCore?>F).core_faust_cplx)
            return core
        raise ValueError("F must be a Faust object")

    cdef _ascomplex(self, scalar=1.0):
        cplx_facs = [self.get_fact_opt(i).astype(np.complex) for i in \
                              range(0,self.get_nb_factors())]
        cplx_facs[0] *= scalar # instead of passing the scal to the
        # construc. It avoids disp of
        # deprecation warning
        core = FaustCore(cplx_facs)#, alpha=scalar)
        core._isReal = False
        return core

    def multiply_scal(self, scalar):
        core = FaustCore(core=True)
        core._isReal = self._isReal
        if(isinstance(scalar, int)):
            scalar = float(scalar)
        scalar_type_err = TypeError("The mul. scalar must be a real or a"
                                    " complex number")
        if(self._isReal):
            if(isinstance(scalar, float)):
                core.core_faust_dbl = \
                        self.core_faust_dbl.mul_scal(scalar)
            elif(isinstance(scalar, np.complex)):
                core = <FaustCore?>self._ascomplex(scalar)
                core._isReal = False
            else:
                raise scalar_type_err
        else:
            if(isinstance(scalar, np.complex) or isinstance(scalar,
                                                            float)):
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_scal(scalar)
            else:
                raise scalar_type_err
        return core

    cdef _vertcat(self, F):
         core = FaustCore(core=True)
         #TODO/ F must be a FaustCore
         if(self._isReal):
             if(not F.isReal()):
                 self = self._ascomplex()
                 core.core_faust_cplx = self.core_faust_cplx.vertcat((<FaustCore?>F).core_faust_cplx)
             else:
                 core.core_faust_dbl = self.core_faust_dbl.vertcat((<FaustCore?>F).core_faust_dbl)
         else:
             core.core_faust_cplx = \
                     self.core_faust_cplx.vertcat((<FaustCore?>F).core_faust_cplx)
         core._isReal = core.core_faust_dbl != NULL
         return core

    def vertcat(self,F):
        if(F.isReal() and not self.isReal()):
            return self._vertcat((<FaustCore?>F)._ascomplex())
        return self._vertcat(F)

    cdef _horzcat(self, F):
         #TODO: refactor with _vertcat(), maybe by passing func ref to a cat func
         core = FaustCore(core=True)
         #TODO/ F must be a FaustCore
         if(self._isReal):
             if(not F.isReal()):
                 self = self._ascomplex()
                 core.core_faust_cplx = self.core_faust_cplx.horzcat((<FaustCore?>F).core_faust_cplx)
             else:
                 core.core_faust_dbl = self.core_faust_dbl.horzcat((<FaustCore?>F).core_faust_dbl)
         else:
             core.core_faust_cplx = \
                     self.core_faust_cplx.horzcat((<FaustCore?>F).core_faust_cplx)
         core._isReal = core.core_faust_dbl != NULL
         return core

    def horzcat(self,F):
        if(F.isReal() and not self.isReal()):
            return self._horzcat((<FaustCore?>F)._ascomplex())
        return self._horzcat(F)

    # Left-Multiplication by a Faust F
    # y=multiply(F,M) is equivalent to y=F*M
    def multiply(self,M):
        if(isinstance(M, FaustCore)):
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

    def set_mul_order_opt_mode(self, mode):
        if(self._isReal):
            self.core_faust_dbl.set_mul_order_opt_mode(mode)
        else:
            self.core_faust_cplx.set_mul_order_opt_mode(mode)

    def set_Fv_mul_mode(self, mode):
        if(self._isReal):
            self.core_faust_dbl.set_Fv_mul_mode(mode)
        else:
            self.core_faust_cplx.set_Fv_mul_mode(mode)


    # print information about the faust (size, number of factor, type of factor (dense/sparse) ...)
    def display(self):
        if(self._isReal):
            self.core_faust_dbl.Display()
        else:
            self.core_faust_cplx.Display()

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

    def normalize(self, ord):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.normalize(ord)
        else:
            core.core_faust_cplx = self.core_faust_cplx.normalize(ord)
        core._isReal = self._isReal
        return core

    def get_nb_factors(self):
        cdef int nb_factors
        if(self._isReal):
            nb_factors = int(self.core_faust_dbl.get_nb_factors())
        else:
            nb_factors = int(self.core_faust_cplx.get_nb_factors())
        return nb_factors

    def get_fact(self,i):
        if(i >=  self.get_nb_factors() or i < 0):
            raise ValueError("factor index must be greater or equal 0 and "
                             "lower than "+str(self.get_nb_factors())+".")
        cdef fact
        cdef double[:,:] fact_dbl_view
        cdef complex[:,:] fact_cplx_view
        if(self._isReal):
            fact = np.zeros([self.core_faust_dbl.get_fact_nb_rows(i),
                             self.core_faust_dbl.get_fact_nb_cols(i)], dtype='d',
                            order='F')
            fact_dbl_view = fact
            self.core_faust_dbl.get_fact(i, &fact_dbl_view[0, 0])
        else:
            fact = np.zeros([self.core_faust_cplx.get_fact_nb_rows(i),
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

    def slice(self, indices):
        # TODO: rename this function or cut in two: slice and fancy indexing
        core = FaustCore(core=True)
        start_row_id, end_row_id, start_col_id, end_col_id = (indices[0].start,
                                                              indices[0].stop,
                                                              indices[1].start,
                                                              indices[1].stop)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.slice(start_row_id,
                                                            end_row_id,
                                                            start_col_id,
                                                            end_col_id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.slice(start_row_id,
                                                              end_row_id,
                                                              start_col_id,
                                                              end_col_id)

        core._isReal = self._isReal
        return core

    def left(self, id):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.left(id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.left(id)
        core._isReal = self._isReal
        return core

    def right(self, id):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.right(id)
        else:
            core.core_faust_cplx = self.core_faust_cplx.right(id)
        core._isReal = self._isReal
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
        if(self._isReal):
            core.core_faust_dbl = \
            self.core_faust_dbl.fancy_idx(&row_indices_view[0], row_indices.size,
                                          &col_indices_view[0], col_indices.size)
        else:
            core.core_faust_cplx = \
            self.core_faust_cplx.fancy_idx(&row_indices_view[0],
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

    def transpose(self):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.transpose()
        else:
            core.core_faust_cplx = self.core_faust_cplx.transpose()
        core._isReal = self._isReal
        return core

    def optimize_storage(self, time=False):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.optimize_storage(time)
        else:
            core.core_faust_cplx = self.core_faust_cplx.optimize_storage(time)
        core._isReal = self._isReal
        return core

    def optimize(self, transp=False):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.optimize(transp)
        else:
            core.core_faust_cplx = self.core_faust_cplx.optimize(transp)
        core._isReal = self._isReal
        return core

    def optimize_time(self, transp=False, inplace=False, nsamples=1):
        if(inplace):
            if(self._isReal):
                self.core_faust_dbl.optimize_time(transp, inplace, nsamples)
            else:
                self.core_faust_cplx.optimize_time(transp, inplace, nsamples)
        else:
            core = FaustCore(core=True)
            if(self._isReal):
                core.core_faust_dbl = self.core_faust_dbl.optimize_time(transp,
                                                                      inplace,
                                                                      nsamples)
            else:
                core.core_faust_cplx = self.core_faust_cplx.optimize_time(transp,
                                                                         inplace,
                                                                         nsamples)
            core._isReal = self._isReal
            return core


    def conj(self):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.conjugate()
        else:
            core.core_faust_cplx = self.core_faust_cplx.conjugate()
        core._isReal = self._isReal
        return core

    def getH(self):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.adjoint()
        else:
            core.core_faust_cplx = self.core_faust_cplx.adjoint()
        core._isReal = self._isReal
        return core

    def zpruneout(self, nnz_tres, npasses, only_forward):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.zpruneout(nnz_tres,
                                                                npasses,
                                                                only_forward)
        else:
            core.core_faust_cplx = self.core_faust_cplx.zpruneout(nnz_tres,
                                                                  npasses,
                                                                  only_forward)
        core._isReal = self._isReal
        return core

    def isReal(self):
        if(self._isReal): return True
        return False

    def __dealloc__(self):
        if(self._isReal):
            del self.core_faust_dbl
        else:
            del self.core_faust_cplx

cdef check_matrix(isReal, M):
        if not isinstance(M, (np.ndarray) ):
            raise ValueError('input must be a numpy ndarray')
        if(isReal):
            M=M.astype(float,'F')
            if not M.dtype=='float':
                raise ValueError('input numpy array dtype must be double (not'
                                 ' float)')
        else:
            M=M.astype(complex,'F')
            if(M.dtype not in ['complex', 'complex128', 'complex64'] ): #could fail if complex128 etc.
                raise ValueError('input array must be complex array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError('input array must be Fortran contiguous (Colmajor)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError('input matrix/array invalid number of dimensions')

cdef class ConstraintIntCore:

    # no need to create an object, because project() needs to create core object
    # for each call and for that purpose it needs to determine if mat is real or complex
    # so it can't create the object before the call
    # in that conditions a static method will suffice
    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, normalized=True,
                pos=False):
        cdef double[:,:] M_view_dbl
        cdef double[:,:] M_out_view_dbl
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx

        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        #order = 'C'
        #if(np.isfortran(M)):
        #    order = 'FORTRAN'
        M = np.asfortranarray(M)

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')


        check_matrix(isReal, M)
        if(isReal):
            M_view_dbl = M
            M_out_view_dbl = M_out
            FaustCoreCy.prox_int[double](name, parameter, &M_view_dbl[0,0], num_rows,
                                         num_cols,&M_out_view_dbl[0,0],
                                         normalized, pos)
        else:
            M_view_cplx = M
            M_out_view_cplx = M_out
            FaustCoreCy.prox_int[complex](name, parameter, &M_view_cplx[0,0],
                                          num_rows, num_cols,
                                          &M_out_view_cplx[0,0], normalized,
                                          pos)

        return M_out

cdef class ConstraintMatCore:

    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, parameter_sz, normalized=False,
                pos=False):
        cdef double[:,:] M_view_dbl
        cdef double[:,:] M_out_view_dbl
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx
        cdef double[:,:] param_view_dbl
        cdef complex[:,:] param_view_cplx

        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']

        M = np.asfortranarray(M)

        M_out = np.empty(M.shape, dtype=M.dtype, order='F')


        check_matrix(isReal, M)
        parameter = parameter.astype(M.dtype)
        check_matrix(isReal, parameter)
        if(isReal):
            M_view_dbl = M
            M_out_view_dbl = M_out
            param_view_dbl = parameter
            FaustCoreCy.prox_mat[double](name, &param_view_dbl[0,0],
                                         parameter_sz, &M_view_dbl[0,0], num_rows,
                                         num_cols,&M_out_view_dbl[0,0],
                                         normalized, pos)
        else:
            M_view_cplx = M
            M_out_view_cplx = M_out
            param_view_cplx = parameter
            FaustCoreCy.prox_mat[complex](name, &param_view_cplx[0,0],
                                          parameter_sz, &M_view_cplx[0,0],
                                          num_rows, num_cols,
                                          &M_out_view_cplx[0,0], normalized,
                                          pos)

        return M_out

def prox_blockdiag(M, block_shapes, normalized, pos):
    cdef double[:,:] M_view_dbl
    cdef double[:,:] M_out_view_dbl
    cdef complex[:,:] M_view_cplx
    cdef complex[:,:] M_out_view_cplx
    cdef unsigned long int* m_ptr
    cdef unsigned long int* n_ptr

    isReal = M.dtype in [ 'float', 'float128',
                         'float16', 'float32',
                         'float64', 'double']

    M = np.asfortranarray(M)

    M_out = np.empty(M.shape, dtype=M.dtype, order='F')

    n_ptr = <unsigned long*>PyMem_Malloc(sizeof(unsigned long *)*len(block_shapes))
    m_ptr = <unsigned long*>PyMem_Malloc(sizeof(unsigned long *)*len(block_shapes))

    for i in range(len(block_shapes)):
        m_ptr[i] = block_shapes[i][0]
        n_ptr[i] = block_shapes[i][1]

    check_matrix(isReal, M)
    # check_matrix(isReal, parameter)
    if(isReal):
        M_view_dbl = M
        M_out_view_dbl = M_out
        FaustCoreCy.prox_blockdiag[double](&M_view_dbl[0,0], M.shape[0],
                                           M.shape[1], &m_ptr[0], &n_ptr[0],
                                          len(block_shapes),normalized, pos,
                                           &M_out_view_dbl[0,0])
    else:
        M_view_cplx = M
        M_out_view_cplx = M_out
        FaustCoreCy.prox_blockdiag[complex](&M_view_cplx[0,0], M.shape[0],
                                           M.shape[1], &m_ptr[0], &n_ptr[0],
                                          len(block_shapes),normalized, pos,
                                           &M_out_view_cplx[0,0])
    PyMem_Free(m_ptr)
    PyMem_Free(n_ptr)

    return M_out

cdef class ConstraintRealCore:

    # no need to create an object, because project() needs to create core object
    # for each call and for that purpose it needs to determine if mat is real or complex
    # so it can't create the object before the call
    # in that conditions a static method will suffice
    @staticmethod
    def project(M, name, num_rows, num_cols, parameter, normalized=False,
                pos=False):
        cdef double[:,:] M_view_dbl
        cdef double[:,:] M_out_view_dbl
        cdef complex[:,:] M_view_cplx
        cdef complex[:,:] M_out_view_cplx

        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']

#        order = 'C'
#        if(np.isfortran(M)):
#            order = 'FORTRAN'
        M = np.asfortranarray(M)
        M_out = np.empty(M.shape, dtype=M.dtype, order='F')

        check_matrix(isReal, M)
        if(isReal):
            M_view_dbl = M
            M_out_view_dbl = M_out
            FaustCoreCy.prox_real[double, double](name, parameter, &M_view_dbl[0,0], num_rows,
                                         num_cols,&M_out_view_dbl[0,0],
                                                  normalized, pos)
        else:
            M_view_cplx = M
            M_out_view_cplx = M_out
            FaustCoreCy.prox_real[complex, double](name, parameter, &M_view_cplx[0,0],
                                          num_rows, num_cols,
                                                   &M_out_view_cplx[0,0],
                                                   normalized, pos)

        return M_out



cdef class FaustFact:

    @staticmethod
    def fact_palm4msa_fft(Lap, p):
        return FaustFact.fact_palm4msa_gen(Lap, p, p.init_D)

    @staticmethod
    def fact_palm4msa(M, p):
        return FaustFact.fact_palm4msa_gen(M,p)

    @staticmethod
    def fact_palm4msa_gen(M, p, init_D=None):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        # if not float nor complex, raise exception
        check_matrix(isReal, M)

        M = np.asfortranarray(M)

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview
        cdef complex[:,:] Mview_cplx

        cdef double[:,:] tmp_mat
        cdef complex[:,:] tmp_mat_cplx

        # views for lambda and optionally D out buffer (FGFT)
        cdef double[:] outbufview
        cdef complex[:] outbufview_cplx
        # only for FGFT
        cdef double[:] init_D_view
        cdef complex[:] init_D_view_cplx

        cdef FaustCoreCy.PyxParamsFactPalm4MSA[double,double]* cpp_params
        cdef FaustCoreCy.PyxParamsFactPalm4MSA[complex,double]* cpp_params_cplx
        cdef PyxStoppingCriterion[double] cpp_stop_crit
        # template parameter is always double (never complex) because no need
        # to have a threshold error of complex type (it wouldn't make sense)
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
                np.zeros([p.constraints[zeros_id]._num_rows,p.constraints[zeros_id]._num_cols])
            if(not isReal):
                p.init_facts[zeros_id] = p.init_facts[zeros_id].astype(np.complex)
            for i in [i for i in range(0, p.num_facts) if i != zeros_id]:
                p.init_facts[i] = np.eye(p.constraints[i]._num_rows,
                                        p.constraints[i]._num_cols)
                if(not isReal):
                    p.init_facts[i] = p.init_facts[i].astype(np.complex)

        if(calling_fft_algo):
            # FFT/FGFT case, we store lambda in first position and the diagonal
            # of D in the next
            _out_buf = np.empty(init_D.shape[0]+1, dtype=M.dtype)
        else:
            # store only lambda as a return from Palm4MSA algo
            _out_buf = np.array([0], dtype=M.dtype)

        if(isReal):
            if(calling_fft_algo):
                cpp_params = new \
                FaustCoreCy.PyxParamsFactPalm4MSAFFT[double,double]()
                init_D_view = init_D
                (<FaustCoreCy.PyxParamsFactPalm4MSAFFT[double,double]*>cpp_params).init_D = &init_D_view[0]
            else:
                cpp_params = new \
                FaustCoreCy.PyxParamsFactPalm4MSA[double,double]()
            Mview=M
            cpp_params.num_facts = p.num_facts
            cpp_params.is_update_way_R2L = p.is_update_way_R2L
            cpp_params.init_lambda = p.init_lambda
            cpp_params.step_size = p.step_size
            cpp_params.grad_calc_opt_mode = p.grad_calc_opt_mode
            cpp_params.norm2_max_iter = int(p.norm2_max_iter)
            cpp_params.norm2_threshold = p.norm2_threshold
            cpp_params.stop_crit = cpp_stop_crit
            cpp_params.init_facts = <double**> \
                    PyMem_Malloc(sizeof(double*)*p.num_facts)
            cpp_params.init_fact_sizes = <unsigned long*> \
            PyMem_Malloc(sizeof(unsigned long)*2*p.num_facts)
            cpp_params.is_verbose = p.is_verbose
            cpp_params.constant_step_size = p.constant_step_size
            outbufview = _out_buf
        else:
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
                if(isReal):
                    cpp_constraints[i] = <PyxConstraintMat[double]*> \
                            PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                    tmp_mat = cons._cons_value
                    (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                            &tmp_mat[0,0]
                    (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter_sz =\
                            cons._cons_value_sz
                else:
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


        if(isReal):
            cpp_params.constraints = cpp_constraints
            cpp_params.num_constraints = len(p.constraints)
        else:
            cpp_params_cplx.constraints = cpp_constraints
            cpp_params_cplx.num_constraints = len(p.constraints)

        for i in range(0,p.num_facts):
            if(isReal):
                tmp_mat = p.init_facts[i]
                cpp_params.init_facts[i] = &tmp_mat[0,0]
                cpp_params.init_fact_sizes[i*2+0] = p.init_facts[i].shape[0]
                cpp_params.init_fact_sizes[i*2+1] = p.init_facts[i].shape[1]
            else:
                tmp_mat_cplx = p.init_facts[i].astype('complex')
                cpp_params_cplx.init_facts[i] = &tmp_mat_cplx[0,0]
                cpp_params_cplx.init_fact_sizes[i*2+0] = p.init_facts[i].shape[0]
                cpp_params_cplx.init_fact_sizes[i*2+1] = p.init_facts[i].shape[1]

        core = FaustCore(core=True)
        if(isReal):
            if(calling_fft_algo):
                core.core_faust_dbl = \
                FaustCoreCy.fact_palm4MSAFFT[double,double](&Mview[0,0],
                                                            M_num_rows,
                                                            M_num_cols,
                                                            <FaustCoreCy.PyxParamsFactPalm4MSAFFT[double,double]*>cpp_params,
                                                           &outbufview[0])
            else:
                core.core_faust_dbl = FaustCoreCy.fact_palm4MSA[double,double](&Mview[0,0], M_num_rows, M_num_cols,
 #           FaustCoreCy.fact_palm4MSA(&Mview[0,0], M_num_rows, M_num_cols,
                                      cpp_params, &outbufview[0])
            core._isReal = True
        else:
            if(calling_fft_algo):
                core.core_faust_cplx = \
                        FaustCoreCy.fact_palm4MSAFFT[complex,double](&Mview_cplx[0,0],
                                                                    M_num_rows,
                                                                    M_num_cols,
                                                                    <FaustCoreCy.PyxParamsFactPalm4MSAFFT[complex,double]*>cpp_params,
                                                                    &outbufview_cplx[0])
            else:
                core.core_faust_cplx = FaustCoreCy.fact_palm4MSA[complex,double](&Mview_cplx[0,0], M_num_rows, M_num_cols,
                                     cpp_params_cplx, &outbufview_cplx[0])
            core._isReal = False
        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)
        if(isReal):
            PyMem_Free(cpp_params.init_facts)
            PyMem_Free(cpp_params.init_fact_sizes)
        else:
            PyMem_Free(cpp_params_cplx.init_facts)
            PyMem_Free(cpp_params_cplx.init_fact_sizes)
#

        if(isReal):
            del cpp_params
        else:
            del cpp_params_cplx

        if(calling_fft_algo):
            return core, np.real(_out_buf[0]), _out_buf[1:]
        else:
            return core, np.real(_out_buf[0])

    @staticmethod
    def fact_hierarchical_fft(U, Lap, p, init_D):
        return FaustFact.fact_hierarchical_gen(U, p, init_D, Lap)

    @staticmethod
    def fact_hierarchical(M, p):
        return FaustFact.fact_hierarchical_gen(M, p)

    @staticmethod
    def fact_hierarchical_gen(M, p, init_D=None, Lap=None):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        M = np.asfortranarray(M)
        check_matrix(isReal, M)
        if(isReal):
            return FaustFact.fact_hierarchical_gen_real(M, p, init_D, Lap)
        else:
            return FaustFact.fact_hierarchical_gen_cplx(M, p, init_D, Lap)

    @staticmethod
    def fact_hierarchical_gen_real(M, p, init_D=None, Lap=None):

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview

        cdef double[:,:] Lapview

        # views for lambda and optionally D out buffer (FGFT)
        cdef double[:] outbufview

        # only for FGFT
        cdef double[:] init_D_view

        cdef double[:,:] tmp_mat

        cdef FaustCoreCy.PyxParamsHierarchicalFact[double,double]* cpp_params
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
            cpp_params = new \
            FaustCoreCy.PyxParamsHierarchicalFactFFT[double,double]()
            Lapview = Lap
            init_D_view = init_D
            (<FaustCoreCy.PyxParamsHierarchicalFactFFT[double,double]*>cpp_params).init_D = &init_D_view[0]
        else:
            cpp_params = new \
            FaustCoreCy.PyxParamsHierarchicalFact[double,double]()
        Mview=M
        cpp_params.num_facts = p.num_facts
        cpp_params.is_update_way_R2L = p.is_update_way_R2L
        cpp_params.init_lambda = p.init_lambda
        cpp_params.step_size = p.step_size
        cpp_params.grad_calc_opt_mode = p.grad_calc_opt_mode
        cpp_params.norm2_max_iter = int(p.norm2_max_iter)
        cpp_params.norm2_threshold = p.norm2_threshold
        cpp_params.stop_crits = cpp_stop_crits
        cpp_params.is_verbose = p.is_verbose
        cpp_params.is_fact_side_left = p.is_fact_side_left
        cpp_params.constant_step_size = p.constant_step_size
        outbufview = _out_buf

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
                cpp_constraints[i] = <PyxConstraintMat[double]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                tmp_mat = cons._cons_value
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                        &tmp_mat[0,0]
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz

            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols

        cpp_params.constraints = cpp_constraints
        cpp_params.num_rows = p.data_num_rows
        cpp_params.num_cols = p.data_num_cols
        cpp_params.num_constraints = len(p.constraints)

        core = FaustCore(core=True)
        if(calling_fft_algo):
            core.core_faust_dbl = \
            FaustCoreCy.fact_hierarchical_fft[double,double](&Mview[0,0],
                                                             &Lapview[0,0], M_num_rows, M_num_cols,
                                                             <FaustCoreCy.PyxParamsHierarchicalFactFFT[double,double]*>cpp_params,
                                                             &outbufview[0])
        else:
            core.core_faust_dbl = \
            FaustCoreCy.fact_hierarchical[double,double](&Mview[0,0], M_num_rows, M_num_cols,
#           FaustCoreCy.fact_hierarchical(&Mview[0,0], M_num_rows, M_num_cols,
                                  cpp_params, &outbufview[0])
        core._isReal = True

        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        PyMem_Free(cpp_stop_crits)

        del cpp_params
        if(core.core_faust_dbl == NULL): raise Exception("fact_hierarchical"
                                                          " has failed.");

        if(calling_fft_algo):
            return core, np.real(_out_buf[0]), _out_buf[1:]
        else:
            return core, np.real(_out_buf[0])

    @staticmethod
    def fact_hierarchical_gen_cplx(M, p, init_D=None, Lap=None):

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
            elif(cons._is_mat_constraint()):
                #print("FaustFact.fact_hierarchical() Matrix Constraint")
                cpp_constraints[i] = <PyxConstraintMat[complex]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[complex]))
                tmp_mat_cplx = cons._cons_value
                (<PyxConstraintMat[complex]*>cpp_constraints[i]).parameter =\
                        &tmp_mat_cplx[0,0]

            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols

        cpp_params_cplx.constraints = cpp_constraints
        cpp_params_cplx.num_rows = p.data_num_rows
        cpp_params_cplx.num_cols = p.data_num_cols
        cpp_params_cplx.num_constraints = len(p.constraints)

        core = FaustCore(core=True)
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
        core._isReal = False
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
    def fact_givens_fgft_sparse(Lap, J, t, verbosity=0, stoppingError = 0.0,
                                errIsRel=True, order='ascend', enable_large_Faust=False):
        from scipy.sparse import spdiags
        cdef double [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef unsigned int Lap_num_rows=Lap.shape[0]
        cdef unsigned int Lap_num_cols=Lap.shape[1]
        cdef double[:] D_view

        if(order == 'ascend'): order = 1
        elif(order == 'descend'): order = -1
        elif(order == 'undef'): order = 0

        data1d=Lap.data.astype(float,'F')
        indices=Lap.indices.astype(np.int32, 'F')
        indptr=Lap.indptr.astype(np.int32, 'F')

        D = np.empty(Lap.shape[0], dtype=Lap.dtype)
        D_view = D

        core = FaustCore(core=True)
        core.core_faust_dbl = FaustCoreCy.fact_givens_fgft_sparse[double,
                                                           double](&data1d[0],
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

        core._isReal = True
        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if(core._isReal and core.core_faust_dbl == NULL or \
           not core._isReal and core.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")


        return core, D

    @staticmethod
    def fact_givens_fgft(Lap, J, t, verbosity=0, stoppingError = 0.0,
                         errIsRel=True, order=1, enable_large_Faust=False):
        isReal = Lap.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        check_matrix(isReal, Lap)

        if(order == 'ascend'): order = 1
        elif(order == 'descend'): order = -1
        elif(order == 'undef'): order = 0
        else: raise ValueError('order argument must be something among'
                               '\'ascend\', \'descend\' or \'undef\'')

        cdef unsigned int Lap_num_rows=Lap.shape[0]
        cdef unsigned int Lap_num_cols=Lap.shape[1]

        cdef double[:,:] Lap_view
        cdef double[:] D_view

        Lap_view = Lap
        D = np.empty(Lap.shape[0], dtype=Lap.dtype)
        D_view = D

        core = FaustCore(core=True)
        core.core_faust_dbl = FaustCoreCy.fact_givens_fgft[double, double](&Lap_view[0,0],
                                                                           Lap_num_rows,
                                                                           Lap_num_cols, J, t,
                                                                           &D_view[0], verbosity,
                                                                           stoppingError,
                                                                           errIsRel,
                                                                           int(order),
                                                                           enable_large_Faust)

        core._isReal = True
        #from scipy.sparse import spdiags
        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if(core._isReal and core.core_faust_dbl == NULL or \
           not core._isReal and core.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")


        return core, D

    @staticmethod
    def fact_givens_fgft_cplx(Lap, J, t, verbosity=0, stoppingError = 0.0,
                         errIsRel=True, order=1, enable_large_Faust=False):
        isReal = Lap.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        check_matrix(isReal, Lap)

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

        core = FaustCore(core=True)
        core.core_faust_cplx = FaustCoreCy.fact_givens_fgft_cplx[complex,double](&Lap_view[0,0],
                                                                                 Lap_num_rows,
                                                                                 Lap_num_cols, J, t,
                                                                                 &D_view[0], verbosity,
                                                                                 stoppingError,
                                                                                 errIsRel,
                                                                                 int(order),
                                                                                 enable_large_Faust)

        core._isReal = False
        #from scipy.sparse import spdiags
        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if(core._isReal and core.core_faust_dbl == NULL or \
           not core._isReal and core.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return core, D

    @staticmethod
    def fact_givens_fgft_sparse_cplx(Lap, J, t, verbosity=0, stoppingError = 0.0,
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

        core = FaustCore(core=True)
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

        core._isReal = False
        #D_spdiag = spdiags(D, [0], Lap.shape[0], Lap.shape[0])
        #return core, D_spdiag
        if(core._isReal and core.core_faust_dbl == NULL or \
           not core._isReal and core.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return core, D

    @staticmethod
    def eigtj(M, nGivens=None, tol=0, relerr=True,  nGivens_per_fac=None, verbosity=0,
          order='ascend', enable_large_Faust=False):
        if(nGivens == None):
            if(tol == 0):
                raise Exception("You must specify nGivens or tol argument"
                        " (to define a stopping  criterion)")
            nGivens = 0
        if(nGivens_per_fac == None): nGivens_per_fac = int(M.shape[0]/2)
        if(isinstance(M, np.ndarray) and \
            not np.allclose(np.matrix(M, copy=False).H, M) or M.shape[0] != M.shape[1]):
            raise ValueError(" the matrix/array must be symmetric or hermitian.")
        if(not isinstance(nGivens, int)): raise TypeError("nGivens must be a int")
        if(not isinstance(nGivens_per_fac, int)): raise TypeError("nGivens_per_fac must be a int")
        nGivens_per_fac = max(nGivens_per_fac, 1)
        if(nGivens > 0): nGivens_per_fac = min(nGivens_per_fac, nGivens)
        tol *= tol # the C++ impl. works on squared norms to measure errors

        M_is_real = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        if(isinstance(M, np.ndarray)):
            M = np.asfortranarray(M)
            if(M_is_real):
                core_obj,D = FaustFact.fact_givens_fgft(M, nGivens, nGivens_per_fac,
                        verbosity, tol,
                        relerr, order, enable_large_Faust)
            else: #complex
                core_obj,D = FaustFact.fact_givens_fgft_cplx(M, nGivens, nGivens_per_fac,
                        verbosity, tol,
                        relerr, order, enable_large_Faust)

        elif(isinstance(M, csr_matrix)):
            if(M_is_real):
                core_obj,D = FaustFact.fact_givens_fgft_sparse(M, nGivens,
                        nGivens_per_fac,
                        verbosity,
                        tol,
                        relerr,
                        order, enable_large_Faust)
            else: #complex
                core_obj,D = FaustFact.fact_givens_fgft_sparse_cplx(M, nGivens, nGivens_per_fac,
                        verbosity, tol,
                        relerr, order, enable_large_Faust)
        else:
            raise TypeError("The matrix to diagonalize must be a"
                            " scipy.sparse.csr_matrix or a numpy array.")
        return D, core_obj

    @staticmethod
    def svdtj(M, J, t, verbosity=0, stoppingError = 0.0,
              errIsRel=True, enable_large_Faust=False):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        check_matrix(isReal, M)

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] M_view
        cdef double[:] S_view

        M_view = np.asfortranarray(M)
        S = np.empty(M.shape[0], dtype=M.dtype)
        S_view = S
        stoppingError *= stoppingError # the C++ impl. works on squared norms to measure errors

        coreU = FaustCore(core=True)
        coreV = FaustCore(core=True)
        FaustCoreCy.svdtj[double, double](&(coreU.core_faust_dbl),
                                          &(coreV.core_faust_dbl),
                                          &S_view[0],
                                          &M_view[0,0],
                                          M_num_rows,
                                          M_num_cols, int(J), int(t),
                                          verbosity,
                                          stoppingError,
                                          errIsRel,
                                          enable_large_Faust)

        coreU._isReal = coreV._isReal = True
        #from scipy.sparse import spdiags
        #S_spdiag = spdiags(S, [0], M.shape[0], M.shape[0])
        #return core, S_spdiag
        if(coreU._isReal and coreU.core_faust_dbl == NULL or \
           not coreU._isReal and coreU.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return coreU, S, coreV

    @staticmethod
    def svdtj_sparse(M, J, t, verbosity=0, stoppingError = 0.0,
                     errIsRel=True, enable_large_Faust=False):
        from scipy.sparse import spdiags
        cdef double [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]
        cdef double[:] S_view


        data1d = M.data.astype(float,'F')
        indices = M.indices.astype(np.int32, 'F')
        indptr = M.indptr.astype(np.int32, 'F')

        S = np.empty(M.shape[0], dtype=M.dtype)
        S_view = S
        stoppingError *= stoppingError # the C++ impl. works on squared norms to measure errors

        coreU = FaustCore(core=True)
        coreV = FaustCore(core=True)
        FaustCoreCy.svdtj_sparse[double,double](
            &(coreU.core_faust_dbl),
            &(coreV.core_faust_dbl),
            &S_view[0],
            &data1d[0],
            &indices[0],
            &indptr[0],
            M.nnz,
            M_num_rows,
            M_num_cols, J, t,
            verbosity,
            stoppingError,
            errIsRel,
            enable_large_Faust)

        coreU._isReal = coreV._isReal = True
        #D_spdiag = spdiags(D, [0], M.shape[0], M.shape[0])
        #return core, D_spdiag
        if(coreU._isReal and coreU.core_faust_dbl == NULL or \
           not coreU._isReal and coreU.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")

        return coreU, S, coreV

    @staticmethod
    def svdtj_cplx(M, J, t, verbosity=0, stoppingError = 0.0,
              errIsRel=True, enable_large_Faust=False):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        check_matrix(isReal, M)

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef complex[:,:] M_view
        cdef complex[:] S_view

        M_view = np.asfortranarray(M)
        S = np.empty(M.shape[0], dtype=M.dtype)
        S_view = S
        stoppingError *= stoppingError # the C++ impl. works on squared norms to measure errors

        coreU = FaustCore(core=True)
        coreV = FaustCore(core=True)
        FaustCoreCy.svdtj_cplx[complex,double](&(coreU.core_faust_cplx),
                                          &(coreV.core_faust_cplx),
                                          &S_view[0],
                                          &M_view[0,0],
                                          M_num_rows,
                                          M_num_cols, int(J), int(t),
                                          verbosity,
                                          stoppingError,
                                          errIsRel, enable_large_Faust)

        coreU._isReal = coreV._isReal = False
        #from scipy.sparse import spdiags
        #S_spdiag = spdiags(S, [0], M.shape[0], M.shape[0])
        #return core, S_spdiag
        if(coreU._isReal and coreU.core_faust_dbl == NULL or \
           not coreU._isReal and coreU.core_faust_cplx == NULL):
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

        coreU = FaustCore(core=True)
        coreV = FaustCore(core=True)
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

        if(coreU._isReal and coreU.core_faust_dbl == NULL or \
           not coreU._isReal and coreU.core_faust_cplx == NULL):
            raise Exception("Empty transform (nGivens is too big ? Set"
                            " enable_large_Faust to True to force the computation).")
        return coreU, S, coreV

    @staticmethod
    def palm4msa2020(M, p, on_gpu=False):
        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview
        cdef double[:,:] tmp_mat
        # view for lambda
        cdef double[:] outbufview

        cdef PyxStoppingCriterion[double] cpp_stop_crit
        # template parameter is always double (never complex) because no need
        # to have a threshold error of complex type (it wouldn't make sense)
        cdef PyxConstraintGeneric** cpp_constraints

        if(M.dtype == np.complex):
            raise TypeError("2020 palm4msa implementation doesn't support"
                            " complex matrices.")

        M = np.asfortranarray(M)

        Mview=M
        _out_buf = np.array([0], dtype=M.dtype)
        _out_buf[0] = p.init_lambda;
        outbufview = _out_buf

        cpp_stop_crit.is_criterion_error = p.stop_crit._is_criterion_error
        cpp_stop_crit.error_threshold = p.stop_crit.tol
        cpp_stop_crit.num_its = p.stop_crit.num_its
        cpp_stop_crit.max_num_its = p.stop_crit.maxiter

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
                cpp_constraints[i] = <PyxConstraintMat[double]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                tmp_mat = cons._cons_value
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                        &tmp_mat[0,0]
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz
            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols
        core = FaustCore(core=True)
        core.core_faust_dbl = \
            FaustCoreCy.palm4msa2020[double](&Mview[0,0], M_num_rows,
                                             M_num_cols,
                                             cpp_constraints,
                                             len(p.constraints),
                                             &outbufview[0],
                                             cpp_stop_crit,
                                             p.is_update_way_R2L,
                                             p.use_csr, p.packing_RL,
                                             p.norm2_max_iter,
                                             p.norm2_threshold,
                                             p.is_verbose,
                                             p.constant_step_size,
                                             p.step_size, on_gpu)
        core._isReal = True

        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        if(core.core_faust_dbl == NULL): raise Exception("palm4msa2020"
                                                          " has failed.");

        return core, np.real(_out_buf[0])

    @staticmethod
    def hierarchical2020(M, p, on_gpu=False):

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview

        # view for lambda
        cdef double[:] outbufview

        cdef double[:,:] tmp_mat

        cdef PyxConstraintGeneric** cpp_constraints

        M = np.asfortranarray(M)
        is_update_way_R2L = p.is_update_way_R2L
        is_fact_side_left = p.is_fact_side_left
        use_csr = p.use_csr
        packing_RL = p.packing_RL
        norm2_max_iter = p.norm2_max_iter
        norm2_threshold = p.norm2_threshold
        cdef PyxStoppingCriterion[double]* cpp_stop_crits
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

        constraints = p.constraints

        if(M.dtype == np.complex):
            raise TypeError("2020 Hierarchical implementation doesn't support"
                            " complex matrices.")
        # store only lambda as a return from Palm4MSA algo
        _out_buf = np.array([0], dtype=M.dtype)
        _out_buf[0] = p.init_lambda;

        Mview=M
        outbufview = _out_buf

        num_constraints = len(constraints)
        num_facts = int(num_constraints/2+1)

        cpp_constraints = \
        <PyxConstraintGeneric**> \
        PyMem_Malloc(sizeof(PyxConstraintGeneric*)*num_constraints)


        for i in range(0,num_constraints):
            cons = constraints[i]
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
                cpp_constraints[i] = <PyxConstraintMat[double]*> \
                        PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                tmp_mat = cons._cons_value
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                        &tmp_mat[0,0]
                (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter_sz =\
                        cons._cons_value_sz

            else:
                raise ValueError("Constraint type/name is not recognized.")
            cpp_constraints[i].name = cons.name
            cpp_constraints[i].num_rows = cons._num_rows
            cpp_constraints[i].num_cols = cons._num_cols

        core = FaustCore(core=True)
        core.core_faust_dbl = \
                FaustCoreCy.hierarchical2020[double](&Mview[0,0], M_num_rows,
                                                     M_num_cols,
                                                     cpp_stop_crits,
                                                     cpp_constraints,
                                                     num_constraints,
                                                     num_facts,
                                                     &outbufview[0],
                                                     is_update_way_R2L,
                                                     is_fact_side_left,
                                                     use_csr, packing_RL,
                                                     norm2_max_iter,
                                                     norm2_threshold,
                                                     p.is_verbose,
                                                     p.constant_step_size,
                                                     p.step_size, on_gpu)
        core._isReal = True

        for i in range(0,num_constraints):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)
        PyMem_Free(cpp_stop_crits)
        if(core.core_faust_dbl == NULL): raise Exception("hierarchical2020"
                                                          " has failed.");

        return core, np.real(_out_buf[0])

import sys, os, pyfaust
# tries to load the libgm library silently,
# if not enabled at build time it will do nothing
FaustCore.enable_gpu_mod(silent=True)
