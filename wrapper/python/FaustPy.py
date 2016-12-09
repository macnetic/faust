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



#~ dim1 = 5
#~ dim2 = 10
#~ dim3 = 7
#~ nb_factor = 2
#~ int_max= 100

#~ print('**** CONFIG FAUST F ****');
#~ print 'dim1 : '+str(dim1) 
#~ print 'dim2 : '+str(dim2)
#~ print 'nb_factor : '+str(nb_factor)


#~ #initialisation de la liste des facteurs
#~ list_factor=[0]*nb_factor
#~ for i in range(nb_factor):
	#~ list_factor[i]=np.random.randint(int_max, size=(dim1,dim1))

#~ list_factor[nb_factor-1]=np.random.randint(int_max, size=(dim1,dim2))


#~ print "*** CONTRUCTOR ***"
#~ F = FaustPy(list_factor)
#~ print "Ok"
#~ print("dim : "+str(F.shape))


#~ Fbis=F.transpose()
#~ print("**** F ****")
#~ F.afficher()
#~ print("**** F trans ****")
#~ Fbis.afficher()
#~ print("**** attributs F ****")
#~ print(F.__dict__)
#~ dir(Fbis)
