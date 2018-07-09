\mainpage FAµST Library API Documentation

\section intro_sec Introduction

The FAµST toolbox provides algorithms and data structures to decompose a given dense matrix into a product of sparse matrices in order to reduce its computational complexity (both for storage and manipulation).

FaµST can be used to:

- speedup / reduce the memory footprint of iterative algorithms commonly used for solving high dimensional linear inverse problems,
- learn dictionaries with an intrinsically efficient implementation,
- compute (approximate) fast Fourier transforms on graphs.

A general introduction to FAµST is available in the following paper:

[1] [Le Magoarou L. and Gribonval R., “Flexible multi-layer sparse approximations of matrices and applications”](https://hal.archives-ouvertes.fr/hal-01167948), Journal of Selected Topics in Signal Processing, 2016.

The next paper can come as complement:

[2][Le Magoarou L. and Gribonval R., Gramfort A., “FAµST: speeding up linear transforms for tractable inverse problems“](https://hal.archives-ouvertes.fr/hal-01156478v1)

The FAµST toolbox was initially released as a Matlab implementation (versions 1.x). A C++ implementation (versions 2.x), including Matlab wrappers, is now available and the object of further developments. A Python interface is being developed.

FAµST is developped at [Rennes INRIA](https://hal.archives-ouvertes.fr/hal-01156478v1). <br>
For more information on the FAµST Project, please visit the website [faust.inria.fr](http://faust.inria.fr). <br> 


<HR>


@INCLUDE_SPECIFIC_DOC@


<HR>
\authors Adrien Leman, Nicolas Bellot, Thomas Gautrais, Hakim H.
\date @DOXYGEN_CURRENT_DATE@

<HR>

