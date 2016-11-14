import numpy as np
import sys
import PyFaust 


fact2=np.array([[1,2,3],[4,5,6]],dtype='d',order='F');
fact1=np.array([[1,2],[3,4],[5,6],[7,8]],dtype='d',order='F');




print "fact1"
print fact1.shape

print "fact2"
print fact2.shape 

list_factor = [fact1,fact2];
print list_factor
F = PyFaust.Faust(list_factor)

NbColF=F.getNbCol();
NbRowF=F.getNbRow();

print "F shape"
print F.shape

print "F nbRow : "
print NbRowF

print "F nbCol : "
print NbColF



X_array_vec=np.array([1,2,3],dtype='d',order='F');
X_array_mat=np.array([[1,2],[3,4],[5,6]],dtype='d',order='F');
X_matrix_vec=np.matrix(X_array_vec);
X_matrix_mat=np.matrix(X_array_mat);

print "X_vec"
print X_array_vec.shape


print "X_mat"
print X_array_mat.shape






#Y_array_vec=F*X_array_vec;


Faust_Y_vec=F*X_array_vec;
Faust_Y_mat=F*X_array_mat;









numpy_Y_vec=X_matrix_vec.transpose();
numpy_Y_mat=X_matrix_mat;

for fact in reversed(list_factor):
	fact_mat=np.matrix(fact);
	numpy_Y_vec=fact_mat*numpy_Y_vec;
	numpy_Y_mat=fact_mat*numpy_Y_mat;
	
	

F.shape

print "Y_mat_numpy"
print numpy_Y_mat


print "Faust Y_vec"
print Faust_Y_mat

print "Numpy_Y_vec"
print numpy_Y_vec

print "Faust_Y_vec"
print Faust_Y_vec






