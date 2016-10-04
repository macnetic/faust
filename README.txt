##########################################################################
####################		FAuST Toolbox					##########
######## Flexible Approximate Multi-Layer Sparse Transform  ############
##########################################################################

##########################################################################
General purpose:

The FAuST toolbox contains a C++ code implementing a general framework 
designed to factorize matrices of interest into multiple sparse factors. 
It contains a template CPU/GPU C++ code and a Matlab wrapper.
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
Installation:

Please refer to the document "./Faust.pdf" to install the FAUST toolbox.
The FAUST toolbox have been tested on multiplatform :
- LINUX (fedora 20, 21, 22, 23 / Ubuntu) 
- MAC
- WINDOWS (windows 7)

Quick install on UNIX :
1-	Unpack the directory. 
2-	mkdir ./build
3-	cd ./build
4-	cmake .. OR ccmake .. (with Graphical User Interface)
5-	make
6-	make install

Warning : 
mex function working with specific version of gcc depend to the platform used. 
For that, please refers to the matworks information: 
http://fr.mathworks.com/support/compilers/R2016a/index.html
##########################################################################


##########################################################################
Main tools :

The main tools offered by the FAuST toolbox are implemented in two mex function: 
[lambda, facts, errors] = mexHierarchical_fact(params)
[lambda, facts] = mexPalm4MSA(params)

These functions are available in the wrapper matlab.   
###### [lambda, facts, errors] = mexHierarchical_fact(params) #######################
#  This mex function runs the hierarchical
#  matrix factorization algorithm (Algorithm 3 of [1])on the specified
#  input matrix, returning the factors in "facts" (cell array of sparse matrices), 
#  the multiplicative scalar in "lambda" and the errors in "errors".
#
#  Required fields in PARAMS:
#  --------------------------
#    'data' - Training data.
#      A matrix to hierarchically factorize.
#    'nfacts' - Number of factors.
#      Specifies the desired number of factors.
#    'cons' - Constraint sets.
#      Specifies the constraint sets in which eafaust_params_palm.ch factor should lie. It
#      should be a cell-array of size 2*(nfacts-1), where the jth columns
#      sets the constraints for the jth factorization in two factors:
#      cons(1,j) specifies the constraints for the left factor and
#      cons(2,j) for the right factor. cons(i,j) should be itself a
#      cell-array of size 1*4 taking this form for a factor of size m*n:
#      {'constraint name', 'constraint parameter', m, n}
#
#  Optional fields in PARAMS:
#  --------------------------
#    'niter1' - Number of iterations for the 2-factorisation.
#      Specifies the desired number of iteration for the factorisations in
#      2 factors. The default value is 500.
#    'niter2' - Number of iterations for the global optimisation.
#      Specifies the desired number of iteration for the global
#      optimisation. The default value is 500.
#    'fact_side' - Side to be factorized iteratively: 1-the left or
#      0-the right. The default value is 0;
#    'verbose' - Verbosity of the function. if verbose=1, the function
#      outputs the error at each iteration. if verbose=0, the function runs
#      in silent mode. The default value is 0.
#    'update_way' - Way in which the factors are updated. If update_way = 1
#      ,the factors are updated from right to left, and if update_way = 0,
#      the factors are updated from left to right. The default value is 0.

###### [lambda, facts] = mexPalm4MSA(params)
#  Factorization of a data matrix into multiple factors using PALM.
#  This mex function runs the PALM algorithm on the
#  specified set of signals (Algorithm 2 of [1]), returning the factors in
#  "facts" and the multiplicative scalar in "lambda".
#
#  Required fields in PARAMS:
#  --------------------------
#    'data' - Training data.
#      A matrix containing the training signals as its columns.
#    'nfacts' - Number of factors.
#      Specifies the desired number of factors.
#    'cons' - Constraint sets.
#      Specifies the constraint sets in which each factor should lie. It
#      should be a cell-array of size 1*nfacts, where the jth column sets
#      the constraints for the jth factor (starting from the left). cons(j)
#      should be itself a cell-array of size 1*4 taking this form for a
#      factor of size m*n:
#      {'constraint name', 'constraint parameter', m, n}
#    'niter' - Number of iterations.
#      Specifies the number of iterations to run.
#    'init_facts' - Initialization of "facts".
#      Specifies a starting point for the algorithm.
#  Optional fields in PARAMS:
#  --------------------------
#    'init_lambda' - Initialization of "lambda".
#      Specifies a starting point for the algorithm. The default value is 1
#    'verbose' - Verbosity of the function. if verbose=1, the function
#      outputs the error at each iteration. if verbose=0, the function runs
#      in silent mode. The default value is 0.
#    'update_way' - Way in which the factors are updated. If update_way = 1
#      ,the factors are updated from right to left, and if update_way = 0,
#      the factors are updated from left to right. The default value is 0.
##########################################################################
class FAUST :
It represents a given dense matrix by a product of sparse matrix (i.e faust)
in order to speed-up multiplication by this matrix. 
##########################################################################

##########################################################################
Demonstration:

This package contains the Matlab code which represented the results of [1].
One demonstration is available in "./DIR_INSTALL/demo/" directory: 
-	Brain_source_localization using MEG data (cf. Sec.V. of [1];)
-	Hadamard_factorization
-	Runtime_comparison

This package contains a CTest tool to execute binary file automatically.
(cf. CMakeLists.txt)
##########################################################################

##########################################################################
Naming conventions in the C++ FAuST project :

	- namespace 		:	Faust::xxx
	- class 			: 	Faust::MyClass	/	Faust::Class
	- attributs			:	m_myAttribut	/	m_attribut
	#- methods			:	my_method()		/	method()
	- methods			:	myMethod()		/	method()

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


