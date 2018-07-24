##############################################################################
##                              Description:                                ##
##  Test the Python wrapper (time comparison for multiplication             ##
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
# test multiplication time

import sys
import os
if len(sys.argv) != 3 :
	raise ValueError('test_pyFaust.py : invalid number of input arguments')

FaustPath=sys.argv[1]
sys.path.append(FaustPath)

FigPath=sys.argv[2]


import scipy 
from scipy import sparse as sp
import numpy as np
plotlib_present=1
try:
   import matplotlib.pyplot as plt
except ImportError:
    print('Python module matplotlib is not present, no figure will be made')
    plotlib_present=0

import time
import pyfaust 


# function : multiplication par un faust sous Scipy
def faust_multiply(list_factor,x):
	y=x
	#for fact in list_factor:
	for fact in reversed(list_factor):
		y=fact*y
	
	return y

nb_mult=100
nb_factor=4
#dim=4096
RCG=10
density=1./RCG
density_per_fact=density/nb_factor
# comparaison dans la justesse des resultats (obliger de travailler avec des entiers)
result_comparison=1;


#dim_2pow_list=[8,9,10,11,12,13]
dim_2pow_list=[8,9,10]
nb_dim=len(dim_2pow_list)
dim_list=[0]*nb_dim

for i in range(nb_dim):
	dim_list[i]=pow(2,dim_2pow_list[i])
	

print("******* CONFIGURATION ***********")
print("dimension matrix : "+str(dim_list))
print("nb factor : "+str(nb_factor))
print("RCG : "+str(RCG))
print("nbr multiplication : "+str(nb_mult))
print("*********************************\n")




t_list=np.zeros((nb_mult,nb_dim,5))#dense,scipy,faust,dense transposee,faust transposee



print("** FAUST VECTOR MULTIPLICATON ***")
for j in range(nb_dim):
	dim=dim_list[j]
	print("dim : "+str(dim))
	list_factor_sp=[0]*nb_factor
	list_factor_dense=[0]*nb_factor
	int_max=100

	#initialisation 
	for i in range(nb_factor):
		if (scipy.__version__ < '0.18.1'):		
			list_factor_sp[i]=int_max*sp.rand(dim,dim,density=density_per_fact,format='csr',dtype=np.float64)
		else:
			list_factor_sp[i]=int_max*sp.random(dim,dim,density=density_per_fact,format='csr',dtype=np.float64)
		#transforme en entier
		if result_comparison:
			list_factor_sp[i]=np.floor(list_factor_sp[i])
		
		list_factor_dense[i]=list_factor_sp[i].toarray()
		#print(list_factor_dense[i])
		
		
		#print("factor"+str(i)+":")
		#print(list_factor_sp[i])

	#print("list_factor_dense")
	#for fact in list_factor_dense:
	#	print("factor ")
	#	print(fact)


	F=pyfaust.Faust(list_factor_dense)
	F_dense=F.toarray()

	if(F.shape[0] != dim) or (F.shape[1] != dim):
		raise ValueError('invalid Faust size')





	for i in range(nb_mult):
		#print(str(i)+"/"+str(nb_mult))
		x=np.random.randint(int_max, size=(dim,1))
		
		t=time.time()
		y_dense=F_dense.dot(x)
		t_list[i,j,0]=time.time()-t
		
		t=time.time()
		y_scipy=faust_multiply(list_factor_sp,x)
		t_list[i,j,1]=time.time()-t


		t=time.time()
		y_Faust=F*x
		t_list[i,j,2]=time.time()-t
		assert((y_dense == y_Faust).all())
		
		
		t=time.time()
		y_dense_trans=np.transpose(F_dense).dot(x)
		t_list[i,j,3]=time.time()-t

		
		t=time.time()
		y_Faust_trans=F.transpose()*x
		t_list[i,j,4]=time.time()-t
		assert((y_dense_trans == y_Faust_trans).all())



		#print(F_dense)
		#print("t_Faust "+str(t_Faust))
		#print("t_scipy "+str(t_scipy))
		#print("t_dense "+str(t_dense))
		if result_comparison:
			if not (y_Faust==y_scipy).all():
				raise ValueError('multiplication : invalid ouput vector Y')
			
			if not (y_Faust==y_dense).all():
				print("Error")
				print(y_Faust)
				print(y_dense)
				raise ValueError('multiplication : dense_multiplication different from Faust one')
			
			if not (y_Faust_trans==y_dense_trans).all():
				print("Error")
				print(y_Faust)
				print(y_dense)
				raise ValueError('multiplication :transpose  dense_multiplication different from Faust one')
			

#print(t_list)
t_list_mean=np.mean(t_list,0)
t_dense=t_list_mean[:,0]
t_scipy=t_list_mean[:,1]
t_faust=t_list_mean[:,2]
t_dense_trans=t_list_mean[:,3]
t_faust_trans=t_list_mean[:,4]


print("************ RESULT *************")
print("dim "+str(dim_list))
print("t_Faust "+str(t_faust))
print("t_scipy "+str(t_scipy))
print("t_dense "+str(t_dense))
print("t_dense_trans "+str(t_dense_trans))
print("t_faust_trans "+str(t_faust_trans))

speed_up_Faust=t_dense/t_faust
speed_up_Scipy=t_dense/t_scipy
speed_up_Faust_trans=t_dense_trans/t_faust_trans
print("speed-up Faust : "+str(speed_up_Faust))
print("speed-up scipy : "+str(speed_up_Scipy))
print("speed-up Faust (transpose) : "+str(speed_up_Faust_trans))


#print("size dim :"+str(dim_list.shape))
print("t_dense :"+str(t_dense.shape))






if(plotlib_present and not "DONT_PYPLOT_FAUST_TIME" in os.environ):
	print("saving figure in "+FigPath)
	tickness=4
	size_font=18
	font_style='bold'
	plt.loglog(dim_list,t_dense, 'r', basex=2,label='dense (numpy)',linewidth=tickness)
	plt.loglog(dim_list,t_faust, 'b', basex=2,label='Faust (pyfaust)',linewidth=tickness)
	plt.loglog(dim_list,t_scipy, 'g', basex=2,label='Faust (scipy)',linewidth=tickness)
	plt.title('Time Comparison : Faust-Vector multiplication',fontweight=font_style)
	plt.legend(loc=2)
	plt.xlabel('dimension',fontsize=size_font,fontweight=font_style)
	plt.ylabel('time (sec)',fontsize=size_font,fontweight=font_style)
	plt.grid(linestyle="--",which="both")
	plt.grid(linestyle="-")
	plt.savefig(FigPath+'/time_comparison.png')
	plt.clf()

	plt.loglog(dim_list,speed_up_Faust, 'b', basex=2,label='Faust (pyfaust)',linewidth=tickness)
	plt.loglog(dim_list,speed_up_Scipy, 'g', basex=2,label='Faust (scipy)',linewidth=tickness)
	plt.title('Speed-up : Faust-Vector multiplication',fontweight=font_style)
	plt.legend(loc=2)
	plt.xlabel('dimension',fontsize=size_font,fontweight=font_style)
	plt.ylabel('speed-up',fontsize=size_font,fontweight=font_style)
	plt.grid(linestyle="--",which="both")
	plt.grid(linestyle="-")
	plt.savefig(FigPath+'/speed-up.png')


print("*********************************")

#
