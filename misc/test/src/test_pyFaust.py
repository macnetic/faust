##############################################################################
##                              Description:                                ##
##  Test the Python wrapper                                                 ##
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
##############################################################################
import sys
if len(sys.argv) != 2 :
	raise ValueError('test_pyFaust.py : invalid number of input arguments')

FaustPath=sys.argv[1]
sys.path.append(FaustPath)

import numpy as np
import PyFaust 

dim1 = 5
dim2 = 10
dim3 = 7
nb_factor = 5
int_max= 100

print('**** CONFIG FAUST F ****');
print 'dim1 : '+str(dim1) 
print 'dim2 : '+str(dim2)
print 'nb_factor : '+str(nb_factor)


#initialisation de la liste des facteurs
list_factor=[0]*nb_factor
for i in range(nb_factor):
	list_factor[i]=np.random.randint(int_max, size=(dim1,dim1))

list_factor[nb_factor-1]=np.random.randint(int_max, size=(dim1,dim2))



######################################
print "*** CONTRUCTOR ***"
F = PyFaust.Faust(list_factor)
print "Ok"




######################################
print "*** DIMENSION ***"
dim1F=F.getNbRow();
dim2F=F.getNbCol();

print("F dim1 : ",dim1F)
print("F dim2 : ",dim2F)

if (dim1F != dim1) | (dim2F != dim2):
	raise ValueError('invalid dimension')
print "Ok"




#####################################
print "*** MULTIPLICATION VECTOR ***"
x=np.random.randint(int_max, size=(dim2,1))
x_old=x.copy()

expected_y=x
for fact in reversed(list_factor):
	#print i
	fact_mat=np.matrix(fact);
	#print "shape factor "+str(fact_mat.shape)
	#print "shape y "+str(expected_Y_vec.shape)
	expected_y=fact_mat*expected_y;

y=F*x

if (y.shape[0] != dim1) | (y.shape[1] != 1):
	print 'expected size :  (' + str(dim1) + ',1)'
	print 'got  ' + str(y.shape[0])+','+str(y.shape[1])
	raise ValueError('multiplication : invalid size of ouput vector')

if not (y==expected_y).all():
	raise ValueError('multiplication : invalid ouput vector y')
	
if not (x_old==x).all():
	raise ValueError('multiplication : input vector x has changed')
print "Ok"




######################################
print "*** MULTIPLICATION MATRIX ***"
X=np.random.randint(int_max, size=(dim2,dim3))
X_old=X.copy()

expected_Y=X
for fact in reversed(list_factor):
	#print i
	fact_mat=np.matrix(fact);
	#print "shape factor "+str(fact_mat.shape)
	#print "shape y "+str(expected_Y_vec.shape)
	expected_Y=fact_mat*expected_Y;

Y=F*X

if (Y.shape[0] != dim1) | (Y.shape[1] != dim3):
	print 'expected size :  (' + str(dim1) + ','+str(dim3)+')'
	print 'got  ' + str(Y.shape[0])+','+str(Y.shape[1])
	raise ValueError('multiplication : invalid size of ouput matrix')

if not (Y==expected_Y).all():
	raise ValueError('multiplication : invalid ouput matrix Y')


if not (X_old==X).all():
	raise ValueError('multiplication : input matrix X has changed')
print "Ok"





