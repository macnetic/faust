# FAuST Toolbox -- Flexible Approximate Multi-Layer Sparse Transform


General purpose
===============

The FAuST toolbox contains a C++ code implementing a general framework
designed to factorize matrices of interest into multiple sparse factors.
It contains a template CPU/GPU C++ code and a Matlab wrapper.
The algorithms implemented here are described in details in [1]- Le Magoarou

For more information on the FAuST Project, please visit the website of the
project: [FAµST website](http://faust.inria.fr)

---

Installation
============

Please refer to the document "./gettingStartedFAuST-version2_0.pdf"
to install the FAUST toolbox.
The FAUST toolbox has been tested on the following environments:
- LINUX (fedora 20, 21, 22, 23 / Ubuntu)
- MACOS X
- WINDOWS (windows 7)

---

Quick install on UNIX
=====================

	Unpack the directory.
	mkdir ./build
	cd ./build
	cmake .. OR ccmake .. (with Graphical User Interface)
	make
	make install

**Warning**:
The Matlab interface of FAuST requires compiling mex files. The mex compiler
compatible with specific versions of gcc depending on the platform used.
For more information, please refer to the [Mathworks website](http://fr.mathworks.com/support/compilers/R2016a/index.html).

---

License
========

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
along with this program.  If not, see [licenses](http://www.gnu.org/licenses/).

---

Contacts
========
	Luc Le Magoarou: luc.le-magoarou@inria.fr
	Remi Gribonval : remi.gribonval@inria.fr
	Nicolas Bellot  : nicolas.bellot@inria.fr
	Adrien Leman    : adrien.leman@inria.fr
	Thomas Gautrais : thomas.gautrais@inria.fr
---

References
==========

[1]	[Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
	approximations of matrices and applications", Journal of Selected
	Topics in Signal Processing, 2016.](https://hal.archives-ouvertes.fr/hal-01167948v1)

