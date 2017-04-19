##########################################################################
########                 README for developer                  ###########
########                    FAuST Toolbox PYTHON wrapper        ##########
######## Flexible Approximate Multi-Layer Sparse Transform      ##########
##########################################################################


##########################################################################
General purpose :
The Python wrapper intefaces the Faust C++ class with Python.


##########################################################################
Required Component :

Python 2 : https://www.python.org/downloads/
The wrapper is not compatible with Python 3.x versions, onlu Python 2.x versions.

Numpy : (http://www.numpy.org/), it's a python package useful for numerical calulation.
Don't worry, very often it's a nativ python package.
However, you may experienced trouble  if your version of numpy is older than 1.9.2 


Cython : http://cython.org/
This module allows to write Python extension in C/C++.




##########################################################################
File's Structure :

*****************
Main Files :
*****************

src/FaustPy.py :
	--> this python module is directly utilised by the user.
	    It allows the handling of FaustPy.Faust class representing our type of matrices

setup.py.in : 
	--> file that will be configured into setup.py in the build 
	directory with CMake.
	A setup.py file is equivalent to a Makefile in C/C++ but for Python. 
	It's that the different paths such include directories, 
	library's flags are set.
	The behaviour of setup.py file is welldocumented on this link :
	https://docs.python.org/2/distutils/configfile.html

        In our case, it compiles the module FaustCorePy 
	(which consists in only one class FaustCorePy.FaustCore)
	This class is used by the FaustPy.Faust class



*****************
Sources Files :
*****************

src/FaustCore.h and src/FaustCore.hpp :
	--> C++ Header file wrapping the  Class Faust::Transform<FPP,Cpu>,          
        the methods are quite the same but this works with pointer             
	which are easier to handle from Cython/Python

src/FaustCoreCy.pxd :
	--> Cython Header making the links between C++ class FaustCpp       
            and Cython environment, works like C/C++ Header, 
	    creates the Cython module named FaustCoreCy.

src/FaustCorePy.pyx :
	--> creates the Python module FaustCorePy
	    Cython source file making the links between :                     
	    the Cython module FaustCoreCy and Python		






##########################################################################
License:

Copyright (2016):          Luc Le Magoarou, Remi Gribonval, 
                      Nicolas Bellot, Adrien Leman, Thomas Gautrais
                                INRIA Rennes, FRANCE
                                http://www.inria.fr/

The FAuST Toolbox is distributed under the terms of the GNU Affero General 
Public License.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public 
License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
##########################################################################


##########################################################################
Contacts:	

	Luc Le Magoarou: luc.le-magoarou@inria.fr
	Remi Gribonval : remi.gribonval@inria.fr
	Nicolas Bellot  : nicolas.bellot@inria.fr
	Adrien Leman    : adrien.leman@inria.fr	
	Thomas Gautrais : thomas.gautrais@inria.fr	
##########################################################################


##########################################################################
References:

[1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
	approximations of matrices and applications", Journal of Selected 
	Topics in Signal Processing, 2016.
	<https://hal.archives-ouvertes.fr/hal-01167948v1>
##########################################################################	

