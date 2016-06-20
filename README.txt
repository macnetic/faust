##########################################################################
####################		FAuST Toolbox			######################
##########################################################################


##########################################################################
General purpose:

The FAuST toolbox contains code implementing a general framework designed 
to factorize matrices of interest into multiple sparse factors. 
The algorithms implemented here are described in details in [1]- Le Magoarou

For more information on the FAuST Project, please visit the website of the 
project: <http://faust.gforge.inria.fr>
##########################################################################


##########################################################################
License:

Copyright (2016):	Luc Le Magoarou, Remi Gribonval
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
Install:

1-	Unpack the directory. 
2-	mkdir ./build 
3-  cd build
4-	cmake ..
5-	make install
##########################################################################


##########################################################################
Main tools :

The main functions offered by the FAuST toolbox are implemented in two files
	- hierarchical_fact.m
        [lambda, facts, errors] = hierarchical_fact(params);
        input:
            - params: Structure contained data matrice, desired number 
                      of factors, constraint sets...   
        output: 
            - lambda: multiplicative scalar
            - facts : estimated factorization (cell array of sparse matrices)	
            - errors
	- palm4MSA.m
        [lambda, facts] = palm4MSA(params)
        input:
            - params: Structure contained data matrice, desired number 
                      of factors, constraint sets...
        output:
            - lambda: multiplicative scalar
            - facts : estimated factorization (cell array of sparse matrices)
##########################################################################


##########################################################################
Demonstration:

This package contains the Matlab code which represented the results of [1].
One demonstration is available in "./misc/demo/" directory: 
-	Source_localization using MEG data (cf. Sec.V. of [1];)

See related ./misc/demo/README file for more detail. 
##########################################################################

##########################################################################
Naming conventions in the FAuST project :

	- namespace 		:	Faust::xxx
	- class 			: 	Faust::MyClass	/	Faust::Class
	- attributs			:	m_myAttribut	/	m_attribut
	- methods			:	m_my_method()	/	m_method()

	- function			:	Faust::my_function()
	- variable			:	myVariable
	- object			:	myObjet

	- files class  		: 	faust_MyClass.x	/	faust_Class.x
	- files function	:	faust_my_function.x

	- files class gpu	: 	faust_MyClass_gpu.x	/	faust_Class_gpu.x
	- files function gpu:	faust_my_function_gpu.x

 	- class template gpu	: 	Faust::MyClass<FPP,gpu>	/
##########################################################################


##########################################################################
Contacts:	
	Luc Le Magoarou: luc.le-magoarou@inria.fr
	Remi Gribonval : remi.gribonval@inria.fr
##########################################################################


##########################################################################
References:

[1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse 
	approximations of matrices and applications", Journal of Selected 
	Topics in Signal Processing, 2016.
	<https://hal.archives-ouvertes.fr/hal-01167948v1>
##########################################################################	


