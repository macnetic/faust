##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython source file making the links between :                     ##
##        the Cython module FaustCoreCy and Python                          ##
##                                                                          ## 
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.gforge.inria.fr>                         ##
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
from libc.string cimport memcpy;
from libcpp cimport bool
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

cdef class FaustCore:

    #### ATTRIBUTE ########
    # classe Cython
    cdef FaustCoreCy.FaustCoreCpp[double] m_faust


    #### CONSTRUCTOR ####
    #def __cinit__(self,np.ndarray[double, mode="fortran", ndim=2] mat):
    def  __cinit__(self,list_factors):
        #print 'inside cinit'
        self.m_faust = FaustCoreCy.FaustCoreCpp[double]();
        cdef double [:,:] data
        cdef unsigned int nbrow
        cdef unsigned int nbcol
        #print 'longueur list'
        #print len(list_factors)
        #print 'avant boucle'
        for factor in list_factors:
            data=factor.astype(float,'F');
            nbrow=factor.shape[0];
            nbcol=factor.shape[1];
            #	print nbrow
            #	print nbcol
            self.m_faust.push_back(&data[0,0],nbrow,nbcol);
            #print(self.__dict__)
            #print 'apres boucle'



    #### METHOD ####
    #~ 	def getNbRow(self):
        #~ 		#return self.m_faust.getNbRow();
        #~ 		(dim1,dim2)=self.shape();
        #~ 		return dim1

#~ 	def getNbCol(self):
    #~ 		(dim1,dim2)=self.shape();
    #~ 		return dim2


    def shape(self,transpose_flag=False):
        cdef unsigned int nbrow = 0
        cdef unsigned int nbcol = 0
        self.m_faust.setOp(transpose_flag,nbrow,nbcol)
        return (nbrow,nbcol)







    #avoid using np.ndarray as input type
    #because there is no distinction between row vector and column vector
    #WARNING : a vector is always a Column Vector
    #
    # numpy.matrix is always 2-dimensional but C-style (RowMajor)  stored
    # not Fortran-style (ColMajor) as Faust lib use Fortran-style storage

    # Left-Multiplication by a Faust F
    # y=multiply(F,x) is equivalent to y=F*x 
    def multiply(self,x,transpose_flag):
        if not isinstance(x, (np.ndarray) ):
            raise NameError('input x must a numpy ndarray')
        #transform into float F continous  matrix
        x=x.astype(float,'F')
        if not x.dtype=='float':
            raise NameError('input x must be double array')
        if not x.flags['F_CONTIGUOUS']:
            raise NameError('input x must be Fortran contiguous (Colmajor)')

        ndim_x=x.ndim;

        if (ndim_x > 2) | (ndim_x < 1):
            raise NameError('input x invalid number of dimensions')

        cdef unsigned int nbrow_x=x.shape[0]
        cdef unsigned int nbcol_x #can't be assigned because we don't know yet if the input vector is 1D or 2D

        dimThis=self.shape(transpose_flag)
        cdef unsigned int nbRowThis=dimThis[0];
        cdef unsigned int nbColThis=dimThis[1];


        cdef unsigned int nbrow_y=nbRowThis
        cdef unsigned int nbcol_y

        cdef double[:] xview_1D
        cdef double[:,:] xview_2D

        if ndim_x == 1:
            nbcol_x=1
            xview_1D=x;
        else:
            nbcol_x=x.shape[1]
            xview_2D=x;

            if (nbrow_x != nbColThis):
                raise NameError('y=F*x multiplication with Faust : invalid dimension of the input matrix x');

        #void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose);
        nbcol_y = nbcol_x;

        cdef y = np.zeros([nbrow_y,nbcol_y], dtype='d',order='F')
        cdef double[:,:] yview=y
        if ndim_x == 1:
            self.m_faust.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_1D[0],nbrow_x,nbcol_x,transpose_flag);
        else:
            self.m_faust.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_2D[0,0],nbrow_x,nbcol_x,transpose_flag);

        return y



    # print information about the faust (size, number of factor, type of factor (dense/sparse) ...)	
    def display(self,transpose_flag):
        print("Faust transposition " + str(transpose_flag))
        self.m_faust.Display();

    def nnz(self):
        cdef unsigned long long nnz = 0
        nnz = self.m_faust.nnz()
        return nnz

    def norm(self):
        cdef double norm
        norm = self.m_faust.norm()
        return norm

    def get_nb_factors(self):
        cdef int nb_factors
        nb_factors = int(self.m_faust.get_nb_factors())
        return nb_factors

    def get_fact(self,i):
        cdef fact = np.zeros([self.m_faust.get_fact_nb_rows(i),
                              self.m_faust.get_fact_nb_cols(i)], dtype='d',
                             order='F')
        cdef double[:,:] fact_view = fact
        self.m_faust.get_fact(i, &fact_view[0, 0])
        return fact

    def save_mat_file(self,filepath):
        cdef char * cfilepath = <char*> PyMem_Malloc(sizeof(char) *
                                                     (len(filepath)+1))
        fparr = bytearray(filepath, "UTF-8");
        for i in range(0,len(filepath)):
            cfilepath[i] = fparr[i]
        cfilepath[i+1] = 0
        self.m_faust.save_mat_file(cfilepath)
        PyMem_Free(cfilepath)
