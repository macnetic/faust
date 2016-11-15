##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython source file making the links between :                     ##
##        the Cython module CyFaust and Python                              ##
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
cimport CyFaust 

import numpy as np 
cimport numpy as np

from libc.stdlib cimport malloc, free;
from libc.string cimport memcpy;
from libcpp cimport bool

cdef class Faust:
	
	#### ATTRIBUTE ########
	# classe Cython
	cdef CyFaust.FaustCpp[double] m_faust
	
	#### CONSTRUCTOR ####
	#def __cinit__(self,np.ndarray[double, mode="fortran", ndim=2] mat):
	def  __cinit__(self,list_factors):
		#print 'inside cinit'
		self.m_faust = CyFaust.FaustCpp[double]();
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
		#print 'apres boucle'
		
	
	
	#### METHOD ####
	def getNbRow(self):
		return self.m_faust.getNbRow();
		
	def getNbCol(self):
		return self.m_faust.getNbCol();
		
		
	def shape(self):
		return (self.m_faust.getNbRow(),self.m_faust.getNbRow())
	





	#avoid using np.ndarray as input type
	#because there is no distinction between row vector and column vector
	#WARNING : a vector is always a Column Vector
	#
	# numpy.matrix is always 2-dimensional but C-style (RowMajor)  stored
	# not Fortran-style (ColMajor) as Faust lib use Fortran-style storage

	# Left-Multiplication by a Faust F
	# y=multiply(F,x) is equivalent to y=F*x 
	def multiply(self,x):
		if not isinstance(x, (np.ndarray) ):
			raise NameError('input x must a numpy ndarray')
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
		
		cdef unsigned int nbRowThis=self.getNbRow();
		cdef unsigned int nbColThis=self.getNbCol();
		
		
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
			
		
		#void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,int nbrow_x,int nbcol_x,bool isTranspose);
		nbcol_y = nbcol_x;
		
		cdef y = np.zeros([nbrow_y,nbcol_y], dtype='d',order='F')
		cdef double[:,:] yview=y
		if ndim_x == 1:
			self.m_faust.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_1D[0],nbrow_x,nbcol_x,False);
		else:
			self.m_faust.multiply(&yview[0,0],nbrow_y,nbcol_y,&xview_2D[0,0],nbrow_x,nbcol_x,False);
		
		return y
		
		
	#overloading of multiplication operator *,
	# y = F * x ,with F a Faust
	def __mul__(self, x):
		return self.multiply(x)
