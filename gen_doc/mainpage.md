*! \mainpage FAµST mainpage
\section intro_sec Introduction
FAµST is a C++ library that implements a way to decompose a given matrix into a product of sparse matrix in order to reduce multiplication costs by the matrix. <br>
FaµST can be widely used to speed up iterative algorithm commonly used for solving high dimensional linear inverse problems. <br>
<br>
FAµST is delivered with Matlab wrapper. <br>
<br>
FAµST is developped at <a href="http://www.inria.fr/en/centre/rennes"> Rennes INRIA</a> by <a href="https://team.inria.fr/panama/fr/">Panama team </a>. <br>


\section install_sec Installation

\subsection step1 Step 1: download package libFaust
\subsection step2 Step 2: ./configure
\subsection step1 Step 3: ./make


\section Demo
Faust was used in a experience of source localization in brain image.<br>
Different FAUST approximation of a Magnetoencephalography (MEG) gain matrix (M<i>) was computed using hierarchical_fact.<br>
A SparseCoding algorithm was used to solve this source localization problem using the MEG gain matrix or its FAUST approximation.<br>
The following pictures illustrate the different trade-offs between speed-up and error of localization using MEG matrix or FAUST matrix.<br>


\image html  MEG_computed_time.jpg "computing time" width=5cm
\image html MEG_distance.jpg "distance between estimated sources and true ones" width=2cm
