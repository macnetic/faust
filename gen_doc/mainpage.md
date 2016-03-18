\mainpage FAµST

\section intro_sec Introduction
FAµST is a C++ library that implements a way to decompose a given matrix into a product of sparse matrix in order to reduce multiplication costs by the matrix. <br>
FaµST can be widely used to speed up iterative algorithm commonly used for solving high dimensional linear inverse problems. <br>
It is delivered with Matlab wrapper. <br>
<br>
The description of Faust's project is available in articles written by Luc Le Magoarou and Rémi Gribonval at : <a href="https://hal.archives-ouvertes.fr/hal-01167948v1"> link1 </a> and <a href="https://hal.archives-ouvertes.fr/hal-01156478v1"> link2 </a>. <br>
FAµST is developped at <a href="http://www.inria.fr/en/centre/rennes"> Rennes INRIA</a> by <a href="https://team.inria.fr/panama/fr/">Panama team </a>. <br>
For more information on the FAµST Project, please visit the <a href="http://faust.gforge.inria.fr"> website</a>. <br> 

\section install_sec Installation

\subsection step1 Step 1: Download FAµST
FAµST Project is available on subversion repository. 
	- svn co ./FAUST  <br>

\subsection step2 Step 2:  Dependency 
	- library  OpenBLAS http://www.openblas.net/
	- library  HDF5 https://www.hdfgroup.org/HDF5/release/obtainsrc.html#conf
	- library  MATIO is an C library for reading and writing MATLAB MAT files. https://sourceforge.net/projects/matio/ <br>
		   matio must be build with hdf5 : Support for MAT file version 7.3 requires the HDF5 library <br>
		   ./configure --enable-extended-sparse=yes --with-matlab=matlab --with-hdf5=/home/library/hdf5-1.8.16/src/.libs <br>
		   make <br>
		   make check <br>
	- MATLAB   Version 2014b <br>
	- GCC 4.7  Comptatible with GCC matlab for mex tools <br>


\subsection step3 Step 3: Configure
	- ./cmake ../..  <br>
	[option] : <br>
		-G "CodeBlocks - Unix Makefiles"	-->	indicates that create both makefile and codeBlock projets. <br>
		-DCMAKE_BUILD_TYPE="Debug" 		-->	indicates that compilation is in debug mode.  <br>

\subsection step4 Step 4: Make
	- make -j4 		--> 	build with 4 threads <br>
 	- build from codeBlock 	--> 	In debug mode for example <br>

\section Example_sec Target & Example

Here is some examples of algorithm of the project FAUST. <br>

\subsection Target1 Target 1 : faust_hier
hierarchical_fact_test.cpp <br>
Run the hierarchical factorization from a dense matrix contained in a file. <br>
target name : testing/bin/faust_hier <br>

\subsection Target2 Target 2 : launch_hierarchical_fact
configuration file : cmd_line/cmdline_function/launch_hierarchical_fact.cpp <br>
wrapper/cmd_line/src/launch_hierarchical_fact.cpp <br>
Compute the hierarchical factorization of a given data matrix A in cmdline mode. <br>
target name : wrapper/cmd_line/bin/launch_hierarchical_fact <br>

\subsection Target3 Target 3 : launch_palm4MSA
configuration file : cmd_line/cmdline_function/launch_palm4MSA.cpp <br>
wrapper/cmd_line/src/launch_palm4MSA.cpp <br>
Launch_palm4MSA command-line factorizes an input matrix corresponding to an operator of interest into J sparse factors and converges to a stationary point. <br>
target name : wrapper/cmd_line/bin/launch_palm4MSA <br>

\subsection Target4 Target 4 : faust_test
faust_test.cpp <br>
	 - Compute the multiplication between an faust-matrix and a vector.  <br>
	 - Evaluation of the time difference using faust-matrix and classical matrix (apply on MEG experiment or dictionnary learning).  <br>

target name : testing/bin/faust_test <br>

\subsection Target5 Target 5 : meg
MEG_fact.cpp <br>
Run the hierarchical factorization of the MEG data <br>
target name : meg <br>

\subsection Target6 Target 6 : palm
palm4MSA_test.cpp <br>
Run a test of palm4MSA <br>
target name : testing/bin/palm <br>




\section Demo Demonstration
\subsection Demo1 Demo 1 : Source localization in the context of functional brain imaging :
Faust was used in a experience of source localization in brain image.<br>
Different FAUST approximation of a Magnetoencephalography (MEG) gain matrix "Mi" was computed using hierarchical_fact.<br>
A SparseCoding algorithm was used to solve this source localization problem using the MEG gain matrix or its FAUST approximation.<br>
The following pictures illustrate the different trade-offs between speed-up and error of localization using MEG matrix or FAUST matrix.<br>

\image html MEG_computed_time.jpg "computing time" width=2cm
\image html MEG_distance.jpg "distance between estimated sources and true ones" width=2cm
in progress...<br>

\subsection Demo2 Demo 2 : Image denoising :
The goal of this demo is to enhancement a noisy image.<br>
in progress...<br>


\authors Nicolas Bellot, Thomas Gautrais, Adrien Leman
\date 03/2016



