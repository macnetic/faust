##############################################################################
##                              Description:                                ##
##  Test the Python wrapper                                                 ##
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
##############################################################################
import sys
if len(sys.argv) != 2 :
	raise ValueError('test_pyFaust.py : invalid number of input arguments')

FaustPath=sys.argv[1]
sys.path.append(FaustPath)

import numpy as np
#import PyFaust 
import pyfaust
dim1 = 100
dim2 = 200
dim3 = 7
nb_factor = 3
int_max= 100


# function : multiplication par un faust representer par une liste de facteur sous Python
def faust_multiply(list_factor,x):
	y=x
	#for fact in list_factor:
	for fact in reversed(list_factor):
		y=fact.dot(y)
	
	return y


print('**** CONFIG FAUST F ****');
print('dim1 : '+str(dim1) )
print('dim2 : '+str(dim2))
print('nb_factor : '+str(nb_factor))


#initialisation de la liste des facteurs
list_factor=[0]*nb_factor
for i in range(nb_factor):
	list_factor[i]=np.random.randint(int_max, size=(dim1,dim1))

list_factor[nb_factor-1]=np.random.randint(int_max, size=(dim1,dim2))



######################################
print("*** CONTRUCTOR ***")
F = pyfaust.Faust(list_factor)
print("lister les attribut de F")

print("*** display ***")
F.display()
print("Ok")


######################################
print("*** DIMENSION ***")
dim1F=F.shape[0];
dim2F=F.shape[1];

print("F dim1 : ",dim1F)
print("F dim2 : ",dim2F)

if (dim1F != dim1) | (dim2F != dim2):
	raise ValueError('invalid dimension')
print("Ok")



####################################
print("*** TRANSPOSE ***")
F_trans=F.transpose()
if not (F_trans.shape == (dim2, dim1)) :
	print("expected size : "+str([dim2, dim1]))
	print("got : "+str(F_trans.shape))
	raise ValueError('F_trans : invalid  size of the dense matrix')


if not (F.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F.shape))
	raise ValueError('F_trans : modify size invalid  size of the dense matrix')


F_trans_trans=F_trans.transpose()
if not (F_trans_trans.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_trans_trans.shape))
	raise ValueError('F_trans_trans invalid  size of the dense matrix')
	
if not (F_trans.shape == (dim2, dim1)) :
	print("expected size : "+str([dim2, dim1]))
	print("got : "+str(F_trans.shape))
	raise ValueError('F_trans_trans modify size of the dense matrix')

print("Ok")

######################################
print("*** CONVERSION DENSE MATRIX ***")
F_dense_expected=faust_multiply(list_factor,np.eye(dim2))

F_dense=F.toarray()

if not (F_dense.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_dense.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (F_dense==F_dense_expected).all():
	raise ValueError('invalid value')
	


F_trans_dense=F_trans.toarray()
if not (F_trans_dense.shape == (dim2, dim1)) :
	print("expected size : "+str([dim2, dim1]))
	print("got : "+str(F_dense_trans.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (np.transpose(F_dense)==F_trans_dense).all():
	raise ValueError('invalid value')
	
	
F_trans_trans_dense=F_trans_trans.toarray()

if not (F_trans_trans_dense.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_dense.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (F_dense==F_dense_expected).all():
	raise ValueError('invalid value')
print("Ok")






#####################################
print("*** MULTIPLICATION VECTOR ***")
x=np.random.randint(int_max, size=(dim2,1))
x_old=x.copy()

expected_y=F_dense.dot(x)


y=F*x

if (y.shape[0] != dim1) | (y.shape[1] != 1):
	print('expected size :  (' + str(dim1) + ',1)')
	print('got  ' + str(y.shape[0])+','+str(y.shape[1]))
	raise ValueError('multiplication : invalid size of ouput vector')

if not (y==expected_y).all():
	raise ValueError('multiplication : invalid ouput vector y')
	
if not (x_old==x).all():
	raise ValueError('multiplication : input vector x has changed')
print("Ok")




######################################
print("*** MULTIPLICATION MATRIX ***")
X=np.random.randint(int_max, size=(dim2,dim3))
X_old=X.copy()

expected_Y=F_dense.dot(X)


Y=F*X

if (Y.shape[0] != dim1) | (Y.shape[1] != dim3):
	print('expected size :  (' + str(dim1) + ','+str(dim3)+')')
	print('got  ' + str(Y.shape[0])+','+str(Y.shape[1]))
	raise ValueError('multiplication : invalid size of ouput matrix')

if not (Y==expected_Y).all():
	raise ValueError('multiplication : invalid ouput matrix Y')


if not (X_old==X).all():
	raise ValueError('multiplication : input matrix X has changed')
print("Ok")



######################################
print("*** SLICING ***")
for i in range(dim1):
	for j in range(dim2):
		F_i_j=F[i,j]
		F_trans_j_i=F_trans[j,i]
		
		if F_i_j != F_dense[i,j]:
			raise ValueError('invalid value')
			
		if F_i_j != F_trans_dense[j,i]:
			raise ValueError('invalid value')
		
F_dense_slice = F[...,...]
if not (F_dense_slice.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_dense_slice.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (F_dense==F_dense_slice).all():
	raise ValueError('invalid value')


F_trans_dense_slice = F_trans[...,...]
if not (F_trans_dense_slice.shape == (dim2, dim1)) :
	print("expected size : "+str([dim2, dim1]))
	print("got : "+str(F_trans_dense_slice.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (F_trans_dense==F_trans_dense_slice).all():
	raise ValueError('invalid value')
	


F_dense_slice_1 = F[0:dim1,0:dim2]
if not (F_dense_slice_1.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_dense_slice_1.shape))
	raise ValueError('invalid  size of the dense matrix')
if not (F_dense==F_dense_slice_1).all():
	raise ValueError('invalid value')
	
	
#~ slice_row=dim1:0:-1
#~ slice_col=dim2:0:-1

#~ F_dense_slice_2 = F[slice_row,slice_col]
F_dense_slice_2 = F[::-1,::-1]
if not (F_dense_slice_2.shape == (dim1, dim2)) :
	print("expected size : "+str([dim1, dim2]))
	print("got : "+str(F_dense_slice_2.shape))
	print(F_dense_slice_2)
	print("")
	print(F_dense)
	raise ValueError('invalid  size of the dense matrix')
if not (F_dense[::-1,::-1]==F_dense_slice_2).all():
	raise ValueError('invalid value')
print("Ok")
