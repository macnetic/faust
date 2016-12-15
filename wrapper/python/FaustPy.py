##############################################################################
##                              Description:                                ##
##                                                                          ##
##          FaustPy is class wrapping a C++ class                           ##
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

#import sys
#FaustPath='/home/nbellot/Documents/faust_root/faust/trunk/devcpp/build/wrapper/python'
#sys.path.append(FaustPath)
import copy

import numpy as np
import FaustCorePy


class Faust:

	#### CONSTRUCTOR ####
	def  __init__(self,list_factors):
		#print 'inside cinit'
		self.m_faust = FaustCorePy.FaustCore(list_factors);
		self.m_transpose_flag=0;
		self.shape=self.m_faust.shape(self.m_transpose_flag)

	#### METHOD ####
	def getNbRow(self):
		#return self.m_faust.getNbRow();
		return self.shape[0]
		
	def getNbCol(self):
		return self.shape[1]
		
		
	def transpose(self):
		F_trans=copy.copy(self)
		#F_trans=FaustPy([]);
		#F_trans.m_faust=self.m_faust
		F_trans.m_transpose_flag=not (self.m_transpose_flag)
		F_trans.shape=(self.shape[1],self.shape[0])
		
		return F_trans;
		
	def afficher(self):
		print("Struct : ")
		print(self.m_faust)
		print(self.m_transpose_flag)
		print(self.shape)
		
	def __mul__(self,x):

		return self.m_faust.multiply(x,self.m_transpose_flag)
		
	def todense(self):
		identity=np.eye(self.getNbCol(),self.getNbCol());
		self_dense=self*identity
		return self_dense


	def __getitem__(self,list_index):
		#check if list_index has a 2 index (row and column) 
		if (len(list_index) != 2):
			raise ValueError('list_index must contains 2 elements, the row index and the col index')
		
		#check if the index are slice or integer or Ellipsis
		for id in list_index:
			if (not isinstance(id, slice)) and (not isinstance(id,int) and (id !=Ellipsis)):
				raise ValueError('list_index must contains slice (1:n) or Ellipsis (...) or int (2)')
		
		keyCol=list_index[1]
		keyRow=list_index[0]
		identity=np.eye(self.getNbCol(),self.getNbCol());
		identity=identity[...,keyCol]
		submatrix=self*identity
		submatrix=submatrix[keyRow,:]
		return submatrix




