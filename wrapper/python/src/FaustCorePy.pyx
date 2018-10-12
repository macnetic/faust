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

cdef class FaustCore:

    #### ATTRIBUTE ########
    # classe Cython
    cdef FaustCoreCy.FaustCoreCpp[double]* core_faust_dbl
    cdef FaustCoreCy.FaustCoreCpp[complex]* core_faust_cplx
    cdef bool _isReal
    #### CONSTRUCTOR ####
    #def __cinit__(self,np.ndarray[double, mode="fortran", ndim=2] mat):
    def  __cinit__(self,list_factors=None, alpha=1.0, core=False):
        cdef double [:,:] data
        cdef double [:] data1d #only for csr mat factor
        cdef int [:] indices # only for csr mat
        cdef int [:] indptr # only for csr mat
        cdef complex [:,:] data_cplx
        cdef complex [:] data1d_cplx #only for csr mat factor
        cdef unsigned int nbrow
        cdef unsigned int nbcol
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
                   raise ValueError("Faust factors must be a numpy.ndarray or "
                                    "a scipy.sparse.csr.csr_matrix")
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
                      self.core_faust_dbl.push_back(&data[0,0], nbrow, nbcol)
                   else: #csr.csr_matrix
                      data1d=factor.data.astype(float,'F')
                      indices=factor.indices.astype(np.int32, 'F')
                      indptr=factor.indptr.astype(np.int32, 'F')
                      self.core_faust_dbl.push_back(&data1d[0], &indptr[0],
                                                    &indices[0], factor.nnz, nbrow, nbcol)
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
                        self.core_faust_cplx.push_back(&data_cplx[0,0], nbrow, nbcol)
                    else:
                        #print("FaustCore, factor dims:", nbrow, nbcol)
                        data1d_cplx = factor.data.astype(np.complex128, 'F')
                        indices=factor.indices.astype(np.int32, 'F')
                        indptr=factor.indptr.astype(np.int32, 'F')
                        self.core_faust_cplx.push_back(&data1d_cplx[0], &indptr[0],
                                                    &indices[0], factor.nnz, nbrow, nbcol)
        elif(core): # trick to initialize a new FaustCoreCpp from C++ (see
        # transpose, conj and adjoint)
            pass
        #else:
        #TODO: raise error for undefined object here


    @staticmethod
    def randFaust(t,field,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density=0.1):
        core = FaustCore(core=True)
        if(field == 3):
            core.core_faust_dbl = FaustCoreCy.FaustCoreCpp[double].randFaust(t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density)
            core._isReal = True
            if(core.core_faust_dbl == NULL): raise MemoryError()
        elif(field == 4):
            core.core_faust_cplx = FaustCoreCy.FaustCoreCpp[complex].randFaust(t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density)
            if(core.core_faust_cplx == NULL): raise MemoryError()
            core._isReal = False
        else:
            raise ValueError("FaustCorePy.randFaust(): field must be 3 for real or"
                             " 4 for complex")
        return core

    #### METHOD ####
    #~ 	def getNbRow(self):
        #~ 		#return self.core_faust_dbl.getNbRow();
        #~ 		(dim1,dim2)=self.shape();
        #~ 		return dim1

#~ 	def getNbCol(self):
    #~ 		(dim1,dim2)=self.shape();
    #~ 		return dim2


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


    cdef multiply_faust(self, F):
        if(isinstance(F, FaustCore)):
            core = FaustCore(core=True)
            core2 = FaustCore(core=True)
            if(self._isReal):
                core.core_faust_dbl = \
                self.core_faust_dbl.mul_faust((<FaustCore?>F).core_faust_dbl)
            else:
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_faust((<FaustCore?>F).core_faust_cplx)
            core._isReal = self._isReal
            return core
        raise ValueError("F must be a Faust object")


    def multiply_scal(self, scalar):
        core = FaustCore(core=True)
        if(isinstance(scalar, int)):
            scalar = float(scalar)
        if(self._isReal):
            if(isinstance(scalar, float)):
                core.core_faust_dbl = \
                        self.core_faust_dbl.mul_scal(scalar)
            else:
                raise Exception("You cannot multiply a real Faust by a"
                                " complex scalar (not yet implemented).")
        else:
            if(isinstance(scalar, np.complex) or isinstance(scalar,
                                                            float)):
                core.core_faust_cplx = \
                        self.core_faust_cplx.mul_scal(scalar)
            else:
                raise Exception("The multiplicative scalar must be a real or "
                                "a complex number.")
        core._isReal = self._isReal
        return core




    # numpy.matrix is always 2-dimensional but C-style (RowMajor)  stored
    # not Fortran-style (ColMajor) as Faust lib use Fortran-style storage

    # Left-Multiplication by a Faust F
    # y=multiply(F,M) is equivalent to y=F*M
    def multiply(self,M):
        if(isinstance(M, FaustCore)):
            return self.multiply_faust(M)
        if not isinstance(M, (np.ndarray) ):
            raise ValueError('input M must a numpy ndarray')
        if(self._isReal):
           M=M.astype(float,'F')
           if not M.dtype=='float':
               raise ValueError('input M must be double array')
        else:
           M=M.astype(complex,'F')
           if(M.dtype not in ['complex', 'complex128', 'complex64'] ): #could fail if complex128 etc.
               raise ValueError('input M must be complex array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError('input M must be Fortran contiguous (Colmajor)')

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
                                              nbcol_y, &xview_1D_cplx[0], nbrow_x,nbcol_x)
        else:
            if(self._isReal):
                self.core_faust_dbl.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_2D[0,0],nbrow_x,nbcol_x)
            else:
                self.core_faust_cplx.multiply(&yview_cplx[0,0],nbrow_y,nbcol_y,&xview_2D_cplx[0,0],nbrow_x,nbcol_x)

        return y

    

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

    def norm(self, ord):
        cdef double norm
        if(str(ord).lower() not in ["1","2","fro"]):
            raise ValueError("FaustCorePy.norm() invalid type of norm asked.")
        if(self._isReal):
            if(isinstance(ord,int)):
                norm = self.core_faust_dbl.norm(ord)
            else:
                norm = self.core_faust_dbl.normFro()
        else:
            if(isinstance(ord,int)):
                norm = self.core_faust_cplx.norm(ord)
            else:
                norm = self.core_faust_cplx.normFro()
        return norm

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

    def slice(self, start_row_id, end_row_id, start_col_id, end_col_id):
        core = FaustCore(core=True)
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


    def save_mat_file(self,filepath):
        cdef char * cfilepath = <char*> PyMem_Malloc(sizeof(char) *
                                                     (len(filepath)+1))
        fparr = bytearray(filepath, "UTF-8");
        for i in range(0,len(filepath)):
            cfilepath[i] = fparr[i]
        cfilepath[i+1] = 0
        if(self._isReal):
            self.core_faust_dbl.save_mat_file(cfilepath)
        else:
            self.core_faust_cplx.save_mat_file(cfilepath)
        PyMem_Free(cfilepath)

    def transpose(self):
        core = FaustCore(core=True)
        if(self._isReal):
            core.core_faust_dbl = self.core_faust_dbl.transpose()
        else:
            core.core_faust_cplx = self.core_faust_cplx.transpose()
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
            raise ValueError('input M must a numpy ndarray')
        if(isReal):
            M=M.astype(float,'F')
            if not M.dtype=='float':
                raise ValueError('input M must be double array')
        else:
            M=M.astype(complex,'F')
            if(M.dtype not in ['complex', 'complex128', 'complex64'] ): #could fail if complex128 etc.
                raise ValueError('input M must be complex array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError('input M must be Fortran contiguous (Colmajor)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError('input M invalid number of dimensions')


cdef class FaustFact:

    @staticmethod
    def fact_palm4MSA(M, p):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        #TODO: it not float nor complex, raise exception
        check_matrix(isReal, M)

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview
        cdef complex[:,:] Mview_cplx

        cdef double[:,:] tmp_mat
        cdef complex[:,:] tmp_mat_cplx

        cdef FaustCoreCy.PyxParamsFactPalm4MSA[double,double] cpp_params
        cdef FaustCoreCy.PyxParamsFactPalm4MSA[complex,double] cpp_params_cplx
        cdef PyxStoppingCriterion[double] cpp_stop_crit
        # template parameter is always double (never complex) because no need
        # to have a treshold error of complex type (it wouldn't make sense)
        cdef PyxConstraintGeneric** cpp_constraints

        cpp_stop_crit.is_criterion_error = p.stop_crit.is_criterion_error
        cpp_stop_crit.error_treshold = p.stop_crit.error_treshold
        cpp_stop_crit.num_its = p.stop_crit.num_its
        cpp_stop_crit.max_num_its = p.stop_crit.max_num_its

        if(not p.init_facts):
            if(p.is_update_way_R2L):
                zeros_id = p.num_facts-1
            else:
                zeros_id = 0
            p.init_facts = [
                np.zeros([p.constraints[zeros_id]._num_rows,p.constraints[zeros_id]._num_cols]) ]
            if(not isReal):
                p.init_facts[0] = p.init_facts[0].astype(np.complex)
            for i in [i for i in range(0,p.num_facts) if i != zeros_id]:
                p.init_facts += [np.eye(p.constraints[i]._num_rows,
                                        p.constraints[i]._num_cols)]
                if(not isReal):
                    p.init_facts[i] = p.init_facts[i].astype(np.complex)


        if(isReal):
            Mview=M
            cpp_params.num_facts = p.num_facts
            cpp_params.is_update_way_R2L = p.is_update_way_R2L
            cpp_params.init_lambda = p.init_lambda
            cpp_params.step_size = p.step_size
            cpp_params.stop_crit = cpp_stop_crit
            cpp_params.init_facts = <double**> \
                    PyMem_Malloc(sizeof(double*)*p.num_facts)
            cpp_params.init_fact_sizes = <unsigned long*> \
            PyMem_Malloc(sizeof(unsigned long)*2*p.num_facts)
            cpp_params.is_verbose = p.is_verbose
            cpp_params.constant_step_size = p.constant_step_size
        else:
            Mview_cplx=M
            cpp_params_cplx.num_facts = p.num_facts
            cpp_params_cplx.num_facts = p.num_facts
            cpp_params_cplx.is_update_way_R2L = p.is_update_way_R2L
            cpp_params_cplx.init_lambda = p.init_lambda
            cpp_params_cplx.step_size = p.step_size
            cpp_params_cplx.stop_crit = cpp_stop_crit
            cpp_params_cplx.init_facts = <complex**> \
                    PyMem_Malloc(sizeof(complex*)*p.num_facts)
            cpp_params_cplx.init_fact_sizes = <unsigned long*> \
            PyMem_Malloc(sizeof(unsigned long)*2*p.num_facts)
            cpp_params_cplx.is_verbose = p.is_verbose
            cpp_params_cplx.constant_step_size = p.constant_step_size


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
            elif(cons._is_mat_constraint()):
                #print("FaustFact.fact_palm4MSA() Matrix Constraint")
                if(isReal):
                    cpp_constraints[i] = <PyxConstraintMat[double]*> \
                            PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                    tmp_mat = cons._cons_value
                    (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                            &tmp_mat[0,0]
                else:
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
                tmp_mat_cplx = p.init_facts[i]
                cpp_params_cplx.init_facts[i] = &tmp_mat_cplx[0,0]
                cpp_params_cplx.init_fact_sizes[i*2+0] = p.init_facts[i].shape[0]
                cpp_params_cplx.init_fact_sizes[i*2+1] = p.init_facts[i].shape[1]

        core = FaustCore(core=True)
        if(isReal):
            core.core_faust_dbl = FaustCoreCy.fact_palm4MSA[double,double](&Mview[0,0], M_num_rows, M_num_cols,
 #           FaustCoreCy.fact_palm4MSA(&Mview[0,0], M_num_rows, M_num_cols,
                                      &cpp_params)
            core._isReal = True
            #TODO: FPP == complex not yet supported by C++ code
        else:
            core.core_faust_cplx = FaustCoreCy.fact_palm4MSA[complex,double](&Mview_cplx[0,0], M_num_rows, M_num_cols,
                                     &cpp_params_cplx)
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
        return core

    @staticmethod
    def fact_hierarchical(M, p):
        isReal = M.dtype in [ 'float', 'float128',
                             'float16', 'float32',
                             'float64', 'double']
        # double == float64
        #TODO: it not float nor complex, raise exception
        check_matrix(isReal, M)

        cdef unsigned int M_num_rows=M.shape[0]
        cdef unsigned int M_num_cols=M.shape[1]

        cdef double[:,:] Mview
        cdef complex[:,:] Mview_cplx

        cdef double[:,:] tmp_mat
        cdef complex[:,:] tmp_mat_cplx

        cdef FaustCoreCy.PyxParamsHierarchicalFact[double,double] cpp_params
        cdef FaustCoreCy.PyxParamsHierarchicalFact[complex,double] cpp_params_cplx
        cdef PyxStoppingCriterion[double]* cpp_stop_crits
        # template parameter is always double (never complex) because no need
        # to have a treshold error of complex type (it wouldn't make sense)
        cdef PyxConstraintGeneric** cpp_constraints

        cpp_stop_crits = <PyxStoppingCriterion[double]*>\
        PyMem_Malloc(sizeof(PyxStoppingCriterion[double])*2)

        cpp_stop_crits[0].is_criterion_error = p.stop_crits[0].is_criterion_error
        cpp_stop_crits[0].error_treshold = p.stop_crits[0].error_treshold
        cpp_stop_crits[0].num_its = p.stop_crits[0].num_its
        cpp_stop_crits[0].max_num_its = p.stop_crits[0].max_num_its
        cpp_stop_crits[1].is_criterion_error = p.stop_crits[1].is_criterion_error
        cpp_stop_crits[1].error_treshold = p.stop_crits[1].error_treshold
        cpp_stop_crits[1].num_its = p.stop_crits[1].num_its
        cpp_stop_crits[1].max_num_its = p.stop_crits[1].max_num_its



        if(isReal):
            Mview=M
            cpp_params.num_facts = p.num_facts
            cpp_params.is_update_way_R2L = p.is_update_way_R2L
            cpp_params.init_lambda = p.init_lambda
            cpp_params.step_size = p.step_size
            cpp_params.stop_crits = cpp_stop_crits
            cpp_params.is_verbose = p.is_verbose
            cpp_params.is_fact_side_left = p.is_fact_side_left
            cpp_params.constant_step_size = p.constant_step_size
        else:
            Mview_cplx=M
            cpp_params_cplx.num_facts = p.num_facts
            cpp_params_cplx.num_facts = p.num_facts
            cpp_params_cplx.is_update_way_R2L = p.is_update_way_R2L
            cpp_params_cplx.init_lambda = p.init_lambda
            cpp_params_cplx.step_size = p.step_size
            cpp_params_cplx.stop_crits = cpp_stop_crits
            cpp_params_cplx.is_verbose = p.is_verbose
            cpp_params_cplx.is_fact_side_left = p.is_fact_side_left
            cpp_params_cplx.constant_step_size = p.constant_step_size

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
                if(isReal):
                    cpp_constraints[i] = <PyxConstraintMat[double]*> \
                            PyMem_Malloc(sizeof(PyxConstraintMat[double]))
                    tmp_mat = cons._cons_value
                    (<PyxConstraintMat[double]*>cpp_constraints[i]).parameter =\
                            &tmp_mat[0,0]
                else:
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

        if(isReal):
            cpp_params.constraints = cpp_constraints
            cpp_params.num_rows = p.data_num_rows
            cpp_params.num_cols = p.data_num_cols
            cpp_params.num_constraints = len(p.constraints)
        else:
            cpp_params_cplx.constraints = cpp_constraints
            cpp_params_cplx.num_rows = p.data_num_rows
            cpp_params_cplx.num_cols = p.data_num_cols
            cpp_params_cplx.num_constraints = len(p.constraints)

        core = FaustCore(core=True)
        if(isReal):
            core.core_faust_dbl = FaustCoreCy.fact_hierarchical[double,double](&Mview[0,0], M_num_rows, M_num_cols,
 #           FaustCoreCy.fact_hierarchical(&Mview[0,0], M_num_rows, M_num_cols,
                                      &cpp_params)
            core._isReal = True
            #TODO: FPP == complex not yet supported by C++ code
        else:
            core.core_faust_cplx = \
            FaustCoreCy.fact_hierarchical[complex, double](&Mview_cplx[0,0], M_num_rows, M_num_cols,
                                     &cpp_params_cplx)
            core._isReal = False
        for i in range(0,len(p.constraints)):
            PyMem_Free(cpp_constraints[i])
        PyMem_Free(cpp_constraints)

        PyMem_Free(cpp_stop_crits)

        return core
