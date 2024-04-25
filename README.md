[![pipeline status](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/pipeline.svg?ignore_skipped=true)](https://gitlab.inria.fr/faustgrp/faust/commits/hakim_branch)
![pyfaust test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=pyfaust_test_code_coverage&key_text=pyfaustcov&key_width=90)
![matfaust test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=matfaust_test_code_coverage&key_text=matfaustcov&key_width=90)
![C++ test coverage](https://gitlab.inria.fr/faustgrp/faust/badges/hakim_branch/coverage.svg?job=ctest&key_text=C%2B%2Bcov)  

# ![FAµST](./gen_doc/images/logo.png) Toolbox -- Flexible Approximate Multi-Layer Sparse Transform


[[_TOC_]]


### General purpose

The FAuST toolbox contains a C++ code implementing a general framework
designed to factorize matrices of interest into multiple sparse factors.
It contains a template CPU/GPU C++ code and a Matlab wrapper.
A Python wrapper is also available.
The algorithms implemented here are described in details in [1]- Le Magoarou

For more information on the FAuST Project, please visit the website of the
project: [FAµST website](http://faust.inria.fr)

---


### Dependencies

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


### Build on UNIX

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

### Quick Build of the python wrappers (pyfaust) on UNIX (without MATLAB and MATIO)

(With Eigen, Python3 with Cython, numpy and scipy installed)

	git clone git@gitlab.inria.fr:faustgrp/faust.git --depth=1  faust_no_matio
	cd faust_no_matio/
	mkdir build
	cd build
	cmake -DBUILD_WRAPPER_PYTHON=ON -DNO_MATIO=ON -DNOCPPTESTS=ON ..
    # on Linux if clang compiler is not installed add the cmake option -DLINUX_DEFAULT_COMPILER_FOR_PYTHON=ON to defaulty use gcc
	make faust_python

---



### Using Docker for a quick build of the python wrappers (pyfaust) on Linux without any dependency burden

First you need to install [docker](https://docs.docker.com/engine/install/) or ``podman-docker``.
Then follow the next commands:

    git clone https://gitlab.inria.fr/faustgrp/faust.git faust
    cd faust/docker_linux
    # build the docker image, naming it faust_fedora
    docker build -t faust_fedora .
    # run a terminal in the docker container you just built
    docker run -v $PWD/../../faust:/faust:z -it faust_fedora bash
    # the :z in the mapping exp. is only needed if your system is selinux enabled

Then in the bash terminal launched above in the docker container, type the following commands to build and test pyfaust (with its wrappers).

    [root@faust_fedora]# mkdir build_docker && cd build_docker
    [root@faust_fedora]# cmake -DBUILD_WRAPPER_PYTHON=ON -DMATIO_INC_DIR=/usr/include -DMATIO_LIB_FILE=/usr/lib64/libmatio.so ..
    [root@faust_fedora]# make faust_python
    [root@faust_fedora]# cd wrapper/python/
    [root@faust_fedora]# python3 -c "import pyfaust as pf; print(pf.rand(10,10))"
    Faust size 10x10, density 2.5, nnz_sum 250, 5 factor(s):
    - FACTOR 0 (double) SPARSE, size 10x10, density 0.5, nnz 50
    - FACTOR 1 (double) SPARSE, size 10x10, density 0.5, nnz 50
    - FACTOR 2 (double) SPARSE, size 10x10, density 0.5, nnz 50
    - FACTOR 3 (double) SPARSE, size 10x10, density 0.5, nnz 50
    - FACTOR 4 (double) SPARSE, size 10x10, density 0.5, nnz 50

Now you can modify the C++ code of FAµST in ``src`` or the python wrappers in ``wrapper/python``.
Thanks to the docker directory mapping you can do it outside of the docker container
 in your own environment and easily build pyfaust again from the container as explained above.

---


### Quickest Install on Linux, Windows and macOS (pre-built packages)

Please refer to the document [Installation guides](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/md_README.html)
to install the FAUST toolbox.
The FAUST toolbox has been tested on the following environments:

- LINUX (fedora 35 - 38 / Centos 7, 8 / RHEL / Ubuntu)
- MACOS X
- WINDOWS (windows 10)


**Latest pre-compiled release packages** from Gitlab Continuous Integration are also available.
 The links for the latest release are available on the main website
 [install page](https://faust.inria.fr/download/faust-3-x/).  
(All system packages include Matlab and Python wrappers. Of course, PIP
 packages include only Python wrappers)  
You might also refer directly to [PyPI pyfaust](https://pypi.org/project/pyfaust)
 or [Anaconda pyfaust](https://anaconda.org/pyfaust/pyfaust) projects.

**Latest pre-built revision/development packages** (not release ones!) are also available as
 artifacts of the Gitlab pipelines:

- [Revision Linux packages](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build?job=pkg_linux)  
- [Revision Mac OS X](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build?job=pkg_macos)  
- [Revision Windows package](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build?job=pkg_win)  
- [Revision Linux PIP package](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build/wrapper/python/dist?job=pkg_linux_purepy_rev)  
- [Revision Mac OS X PIP package](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build/wrapper/python/dist?job=pkg_macos_purepy_rev)  
- [Revision Windows PIP package](https://gitlab.inria.fr/faustgrp/faust/-/jobs/artifacts/main/browse/build/wrapper/python/dist?job=pkg_win_purepy_rev)  

---


### Contributing to FAµST

Please consult the guide [here](CONTRIBUTING.md) and the [README for developers](README.developer.md).

---


### License -- Credits -- Contacts -- References

#### License

Cf. license.txt

---

#### Contacts

	Rémi Gribonval: remi.gribonval@inria.fr
	Hakim: hakim.hadj-djilani@inria.fr
    Pascal Carrivain: pascal.carrivain@inria.fr


#### Credits

Researchers:
Luc Le Magoarou
Remi Gribonval
TODO: add others

Software engineers:
Adrien Leman (2016), Nicolas Bellot(2015-2016), Thomas Gautrais (2015), Hakim Hadj-Djilani (2018-), Pascal Carrivain (2023-).

---

#### References

[1]	[Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse
	approximations of matrices and applications", Journal of Selected
	Topics in Signal Processing, 2016.](https://hal.archives-ouvertes.fr/hal-01167948v1)

