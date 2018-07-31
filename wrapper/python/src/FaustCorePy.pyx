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

import numpy as np
cimport numpy as np
import copy

from libc.stdlib cimport malloc, free;
from libc.stdio cimport printf
from libc.string cimport memcpy, strlen;
from libcpp cimport bool
from libcpp cimport complex
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from scipy import sparse
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
        elif(field == 4):
            core.core_faust_cplx = FaustCoreCy.FaustCoreCpp[complex].randFaust(t,min_num_factors, max_num_factors, min_dim_size,
                   max_dim_size, density)
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
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError('input M must be Fortran contiguous (Colmajor)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError('input M invalid number of dimensions')

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
