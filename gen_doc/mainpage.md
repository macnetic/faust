\mainpage FAµST Toolbox API Documentation

\section intro_sec What is the FAµST Toolbox ?

The FAµST toolbox provides algorithms and data structures to decompose a given dense matrix into a product of sparse matrices in order to reduce its computational complexity (both for storage and manipulation).

FaµST can be used to:

- speedup / reduce the memory footprint of iterative algorithms commonly used for solving high dimensional linear inverse problems,
- learn dictionaries with an intrinsically efficient implementation,
- compute (approximate) fast Fourier transforms on graphs.

The FAµST toolbox is organized in several parts:

- The FAµST core library, a C++ backend implementing FAµST-related data structures and algorithms,
- [matfaust](./namespacematfaust.html), a Matlab frontend to the FAµST library,
- [pyfaust](./namespacepyfaust.html), a Python frontend to the FAµST library.


@INCLUDE_SPECIFIC_DOC@

\section more_info Learning more about the FAµST Framework

A general introduction to the FAµST framework is available in the following paper:

[1] [Le Magoarou L. and Gribonval R., “Flexible multi-layer sparse approximations of matrices and applications”](https://hal.archives-ouvertes.fr/hal-01167948), Journal of Selected Topics in Signal Processing, 2016.

The following papers can come as complement:

[2] [Le Magoarou L. and Gribonval R., Gramfort A., “FAµST: speeding up linear transforms for tractable inverse problems“](https://hal.archives-ouvertes.fr/hal-01156478v1), European Signal Processing Conference (EUSIPCO), Aug 2015, Nice, France.

[3] [Le Magoarou L., Gribonval R., Tremblay N., “Approximate fast graph Fourier transforms via multi-layer sparse approximation“](https://hal.inria.fr/hal-01416110),Transactions on Signal and Information Processing over Networks

The FAµST toolbox was initially released as a Matlab implementation ([versions 1.x](http://faust.inria.fr/download/faust-1-x/)).
<br/>Since version 2.0, it has been implemented in C++. Besides, the development of wrappers has made this C++ core accessible from Matlab and Python programming languages.


\section Credits

FAµST has been developed in the [PANAMA team](https://team.inria.fr/panama/) and [DANTE team](https://team.inria.fr/dante/). <br>
For further information on the FAµST Project, please visit the website [faust.inria.fr](https://faust.inria.fr). <br>

Researchers: Luc Le Magoarou, Rémi Gribonval, Le Quoc Tung, Amélie Barbe, Léon Zheng, Elisa Riccietti, Mathurin Massias  
Software engineers: Adrien Leman, Nicolas Bellot, Thomas Gautrais, Hakim HADJ-DJILANI


<HR>
\authors  (Documentation) Rémi Gribonval, Hakim HADJ-DJILANI
\date (of doc generation) @DOXYGEN_CURRENT_DATE@

<HR>

