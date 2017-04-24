##############################################################################
##                              Description:                                ##
##    This demo shows that a Faust is handled as matlab builtin matrix,     ##
## presenting  all the function that are overloaded for Faust class         ## 
##                        (size,mtimes,transpose...)                        ##
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
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##############################################################################



#import the FaustPy package
import FaustPy;

# import module to generate data
import scipy
from scipy import sparse as sp
import numpy as np


# generate the factors of the Faust
dim1 = 1000;
dim2 = 2000;
nb_factor = 2;
list_factor_sparse=[0]*nb_factor
list_factor=[0]*nb_factor
int_max = 120
density_per_fact=0.1;
list_factor_sparse[0]=int_max*sp.random(dim1,dim1,density=density_per_fact,format='csr',dtype=np.float64);
list_factor[0]=list_factor_sparse[0].todense();
list_factor_sparse[1]=int_max*sp.random(dim1,dim2,density=density_per_fact,format='csr',dtype=np.float64);
list_factor[1]=list_factor_sparse[1].todense();




# create a Faust named F from its factors
A = FaustPy.Faust(list_factor)

# get the size of the Faust
print "dimension of the Faust : ", A.shape

# transpose a Faust
A_trans = A.transpose()

# multiplication a Numpy matrix by a Faust
x = np.random.randint(int_max, size=(dim2,1))
y = A * x

# convert a faust to numpy matrix
A_numpy = A.todense()

# slicing
coeff = A[0,0]
col_2nd = A[:,1]; 
submatrix_A = A[3:5,2:3]



# speed-up multiplication
import time  
nb_mult = 100
t_dense = 0.0;
t_faust = 0.0;

for i in range(nb_mult):
	x=np.random.randint(int_max, size=(dim2,1))
	
	t_begin = time.time()
	y_dense = A_numpy.dot(x)
	t_elapsed = time.time() - t_begin
	t_dense += t_elapsed
	
	t_begin = time.time()
	y_faust=A*x;
	t_elapsed = time.time()-t_begin 
	t_faust += t_elapsed

print "multiplication SPEED-UP using Faust"
print "Faust is "+str(t_dense/t_faust)+" faster than a full matrix"
