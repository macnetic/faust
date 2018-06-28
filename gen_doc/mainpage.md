\mainpage FAµST Home

\section intro_sec Introduction
 FAµST is a C++ library that implements a way to decompose a given matrix into a product of sparse matrix in order to reduce multiplication costs by the matrix. <br>
FaµST can be widely used to speed up iterative algorithm commonly used for solving high dimensional linear inverse problems. <br>
It is delivered with Matlab wrapper. <br>
<br>
The description of Faust's project is available in articles written by Luc Le Magoarou and Rémi Gribonval at : <a href="https://hal.archives-ouvertes.fr/hal-01167948v1"> link1 </a> and <a href="https://hal.archives-ouvertes.fr/hal-01156478v1"> link2 </a>. <br>
FAµST is developped at <a href="http://www.inria.fr/en/centre/rennes"> Rennes INRIA</a> by <a href="https://team.inria.fr/panama/fr/">Panama team </a>. <br>
For more information on the FAµST Project, please visit the <a href="http://faust.inria.fr"> website</a>. <br> 

<HR>

\section install_sec Installation
See document ./gettingStartedFAuST-versionX_X.pdf to install the toolbox FAµST. 

\section NamingConvention Naming conventions in the C++ FAuST project
	- namespace 		:	Faust::xxx
	- class 			: 	Faust::MyClass	/	Faust::Class
	- attributs			:	m_myAttribut	/	m_attribut
	- methods			:	myMethod()		/	method()

	- function			:	Faust::my_function()
	- variable			:	myVariable
	- object			:	myObjet

	- files class  		: 	faust_MyClass.x	/	faust_Class.x
	- files function	:	faust_my_function.x

	- files class gpu	: 	faust_MyClass_gpu.x	/	faust_Class_gpu.x
	- files function gpu:	faust_my_function_gpu.x

 	- class template gpu	: 	Faust::MyClass<FPP,gpu>	/



<HR>
\authors Adrien Leman, Nicolas Bellot, Thomas Gautrais
\date @DOXYGEN_CURRENT_DATE@

<HR>

