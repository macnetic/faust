[![pipeline status](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/pipeline.svg)](https://gitlab.inria.fr/faustgrp/faust/commits/hakim_branch)
![FAµST Logo](./gen_doc/images/logo.png)

# FAuST Toolbox -- Flexible Approximate Multi-Layer Sparse Transform


General purpose
===============

The FAuST toolbox contains a C++ code implementing a general framework
designed to factorize matrices of interest into multiple sparse factors.
It contains a template CPU/GPU C++ code and a Matlab wrapper.
A Python wrapper is also available.
The algorithms implemented here are described in details in [1]- Le Magoarou

For more information on the FAuST Project, please visit the website of the
project: [FAµST website](http://faust.inria.fr)

---

Dependencies
============
- cuda (preferably cuda 12 latest version but 9 and 11 are also supported).
There is a known bug on cuda 11.4 (issue #305). CUDA is optional, only used if cmake option ``USE_GPU_MOD`` is ON.
The ``gpu_mod`` submodule must be checked out in order to enable this function.
- Eigen 3.4.x.
- [matio](https://github.com/tbeu/matio) version >= 1.5.7 (current latest version 1.5.23 is supported and advised).
matio own dependencies, as hdf5 and zlib. matio dependency can be disabled through cmake option ``NO_MATIO``.
- Python3 (with numpy, scipy and cython) to build the python wrappers.
- Matlab (>= R2017b) to build the matlab wrappers
(there is a constraint on the gcc compiler version depending on the used Matlab version,
the CMake script indicate if the match is not appropriate).
- libxml2 (Optional, needed with CMake ``BUILD_READ_XML_FILE``).
- CMake >= 3.21.0
- TODO: other?


---

Build on UNIX
=====================

	Unpack the directory.
	mkdir ./build
	cd ./build
	cmake .. OR ccmake .. (with Graphical User Interface)
	make
	make install

**Warning 1**:
The Matlab interface of FAuST requires compiling mex files. The mex compiler
compatible with specific versions of gcc depending on the platform used.
For more information, please refer to the [Mathworks website](http://fr.mathworks.com/support/compilers/R2016a/index.html).

**Warning 2**:
Many CMake build options are available (cf. [CMakeLists.txt](./CMakeLists.txt)). It might be quite complicated to deal with them at start (refer to gitlab ci building scripts in [./misc/continuous\_integration/jobs/](./misc/continuous_integration/jobs/) or the [.gitlab-ci.yml](.gitlab-ci.yml) root script to get some insight).

---
Quick Build of the python wrappers (pyfaust) on UNIX (without MATLAB and MATIO)
=========================================================

(With Eigen, Python3 with Cython, numpy and scipy installed)

	git clone git@gitlab.inria.fr:faustgrp/faust.git --depth=1  faust_no_matio
	cd faust_no_matio/
	mkdir build
	cd build
	cmake -DBUILD_WRAPPER_PYTHON=ON -DNO_MATIO=ON -DNOCPPTESTS=ON ..
    # on Linux if clang compiler is not installed add the cmake option -DLINUX_DEFAULT_COMPILER_FOR_PYTHON=ON to defaulty use gcc
	make faust_python

---

Quickest Install on Linux, Windows and macOS (pre-built pakages)
============================================

Please refer to the document [Installation guides](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/md_README.html)
to install the FAUST toolbox.
The FAUST toolbox has been tested on the following environments:
- LINUX (fedora 35 - 37 / centos 7, 8 / Ubuntu)
- MACOS X
- WINDOWS (windows 10)

Pre-compiled packages from Gitlab Continuous Integration are also available. Except of course PIP packages, all packages include matlab and python wrappers, below are the latest release links.  
- [macOS (.pkg) installer](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_macos_release)  
- [Windows (.exe) NSI installer](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_win_release)  
- [Linux (.rpm, .deb) packages](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_linux_release)  
- [Linux (.rpm, .deb) packages with embedded static matio library](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_linux_release)  
- Python PIP (pre-compiled) packages: for [Linux](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_linux_purepy_release), [macOS](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_macos_purepy_release) and [Windows 10](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/master/download?job=package_win_purepy_release). Or preferably, refer to [pypi pyfaust](https://pypi.org/project/pyfaust) or [anaconda pyfaust](https://anaconda.org/pyfaust/pyfaust).

---
License
========

Cf. license.txt

---

Contacts
========

	Rémi Gribonval: remi.gribonval@inria.fr
	Hakim: hakim.hadj-djilani@inria.fr


Credits
========

Researchers:
Luc Le Magoarou
Remi Gribonval
TODO: add others

Software engineers:
Adrien Leman (2016), Nicolas Bellot(2015-2016), Thomas Gautrais (2015), Hakim Hadj-Djilani (2018-), Pascal Carrivain (2023-).

---

References
==========

[1]	[Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
	approximations of matrices and applications", Journal of Selected
	Topics in Signal Processing, 2016.](https://hal.archives-ouvertes.fr/hal-01167948v1)

