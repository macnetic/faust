# Readme for developers


[[_TOC_]]

### General purpose and features

The FAuST toolbox contains a C++ code implementing a general framework
designed to factorize matrices of interest into multiple sparse factors.
It contains a template CPU/GPU C++ code and a Matlab wrapper.
The algorithms implemented here are described in details in [1]- Le Magoarou
For more information on the FAuST Project, please visit the website of the
project: [https://faust.inria.fr](https://faust.inria.fr)


### Installation

[https://faust.inria.fr](https://faust.inria.fr) (for pre-built packages) and README.md (for building)


### FAuST structure directory

- ``./CMake/``
	 contains ".cmake" files used to execute some internal cmake
	command as to find the path of the externals library, Matlab Path,
	Python Path, ...
- ``.gen_doc/``
     contains both user and developer documentation:
		- Developer documentation with Doxygen tool
		- User tutorials: jupyter notebooks and matlab live scripts
- ``.gpu_mod``
     git submodule for GPU implementations.

- ``.misc/``
     contains the tests (Ctest tool), the demonstrations, the data (or a script to download them),
	the configuration files, and Continuous
	Integration job scripts (used in .gitlab-ci.yml)
- ``.src/``
     contains the C++ sources of the project
- ``/wrapper/``
     contains the wrapper of the project (Matlab, Python)




### API Documentation: Doxygen

The Doxygen documentation is available in the following directory:
``./gen_doc/``
To build Doxygen documentation, build the project with following option:
``cmake .. -DBUILD_DOCUMENTATION="ON"``
The configuration of Doxygen are defined in the file ``./gen_doc/Doxyfile.in``
Main page of doxygen is available in following link :
``./build/doc/html/index.html``

Specific wrapper API documentation has been added and is also generated through doxygen
when wrapper building cmake variables are set. And more precisely the C++ doxydoc can be excluded if the
cmake variable ``EXCLUDE_FAUST_LIB_INSTALL`` is set (this feature is used for release packages).
The documentation is automatically uploaded to gitlab each time a (release) tag is added to the gitlab project.

The online Doxygen doc is available
[here](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/index.html).




### Continuous Integration

The Continuous Integration for the project FAUST is based on the CDash tool
(see. http://www.cdash.org/). The building and test are available on the
public link:
[CDash FAµST project](https://cdash-ci.inria.fr/index.php?project=faust)
CDash is not so used anymore, the test reports are directly uploaded on gitlab (cf. README-tests-ci-release.md)

For more details concerning CI (Virtual Machines, etc.), refer to ci.inria.fr and cloud stack.
https://ci.inria.fr project faust2


### Contributing to FAµST

Please consult the guide [here](CONTRIBUTING.md)


### References

[1]	Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
	approximations of matrices and applications", Journal of Selected
Topics in Signal Processing, 2016.
<https://hal.archives-ouvertes.fr/hal-01167948v1>


