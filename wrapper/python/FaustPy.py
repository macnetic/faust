##############################################################################
##                              Description:                                ##
##                                                                          ##
##          FaustPy is a python module  which delivered                     ##
##          a class named Faust which represents a dense matrix             ##
##          by a product of 'sparse' factors (i.e Faust)                    ##
##          Python wrapper class implemented in C++                         ##
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


import copy

import numpy as np
import FaustCorePy


class Faust:
	""" This class represents a dense matrix by a product of 'sparse' factors (i.e Faust) 
		The aim of the Faust representatio is to speed-up multiplication by this matrix
	"""

	def  __init__(F,list_factors):
		""" create a Faust from a list of factor.
		
		Parameter
		---------
		list_factors : list/tuple of numpy matrices
		"""
		F.m_faust = FaustCorePy.FaustCore(list_factors);
		F.m_transpose_flag=0;
		F.shape=F.m_faust.shape(F.m_transpose_flag)


	def getNbRow(F):
		""" return the number of row of the current Faust. """
		return F.shape[0]
		
	def getNbCol(F):
		""" return the number of column of the current Faust. """
		return F.shape[1]
		
		
	def transpose(F):
		""" transpose the current Faust. """
		F_trans=copy.copy(F)
		F_trans.m_transpose_flag=not (F.m_transpose_flag)
		F_trans.shape=(F.shape[1],F.shape[0])
		
		return F_trans;
		
	def display(F):
		""" display information of the current Faust. """
		print("Struct : ")
		F.m_faust.display(F.m_transpose_flag);
		
		
	def __mul__(F,x):
		""" multiplication between the current Faust and a numpy matrix.
		(overload the python operator *, return F*x)
		
		Parameter
		---------
		x : 2D numpy ndarray of double scalar, must be FORTRAN contiguous 
		"""
		return F.m_faust.multiply(x,F.m_transpose_flag)
		
	def todense(F):
		""" convert the current Faust into a numpy matrix. 
			return a numpy matrix """
		
		
		identity=np.eye(F.getNbCol(),F.getNbCol());
		F_dense=F*identity
		return F_dense


	def __getitem__(F,list_index):
		""" Slicing : Return the value of the current Faust at index list_index.
		(overload of python built-in operator F(indexRow,indexCol)
		
		Parameter
		---------
		list_index : tab of length 2 with its elements must be slice, integer or Ellipsis(...) 
		
		Example of use :
		F[2,3], F[0:dim1,...], F[::-1,::-1] 
		
		"""
		#check if list_index has a 2 index (row and column) 
		if (len(list_index) != 2):
			raise ValueError('list_index must contains 2 elements, the row index and the col index')
		
		#check if the index are slice or integer or Ellipsis
		for id in list_index:
			if (not isinstance(id, slice)) and (not isinstance(id,int) and (id !=Ellipsis)):
				raise ValueError('list_index must contains slice (1:n) or Ellipsis (...) or int (2)')
		
		keyCol=list_index[1]
		keyRow=list_index[0]
		identity=np.eye(F.getNbCol(),F.getNbCol());
		identity=identity[...,keyCol]
		submatrix=F*identity
		submatrix=submatrix[keyRow,:]
		return submatrix




