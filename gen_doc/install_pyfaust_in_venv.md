\page install_pyfaust_in_venv Installing pyfaust in Python virtual environments (Python venv & Anaconda)

## SUMMARY
[1. Python @PY3_VER@](#python_venv)  
[2. Anaconda](#anaconda)

\anchor python_venv
1. Python @PY3_VER@
================

First, you need to create your virtual environment into a directory reserved to pyfaust:

	$ python3 -m venv test_pyfaust-@CPACK_PACKAGE_VERSION@

Second, you need to enable the virtual environment:

	$ source ./test_pyfaust-@CPACK_PACKAGE_VERSION@/bin/activate

\note For Windows users the command is rather:

	C:\> call .\test_pyfaust-@CPACK_PACKAGE_VERSION@\Scripts\activate

Then, if all worked correctly you should see the test_pyfaust-@CPACK_PACKAGE_VERSION@ prompt of your virtual environment.

It remains now to install pyfaust in this environment.

If you downloaded the ``whl`` package the command is:
	$ pip install pyfaust-*.whl

Otherwise, a simpler way is to directly install pyfaust from the public PYPI repository:
	$ pip install pyfaust

All the dependencies will be downloaded and installed by the command above.

You can launch this python command to validate the installation:

	$ python -c 'from pyfaust import version; print(version())'

As a result, you should see the version of pyfaust installed if all went properly.

\note In the virtual environment python is necessarily python version 3 (because the creation was made through python3 above).

\note On Mac OS X you might need to install OpenMP for pyfaust to work, for further information please look at this [FAQ entry](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/FAQ.html#py_four).

\anchor anaconda
2. Anaconda
===========

Within Anaconda you can also create a virtual environment and install pyfaust there. It's quite similar to the Python way described in the previous section. In that purpose you'll need to use the command conda-create. Please type the following commands in an Anaconda prompt or a shell prompt:

Create the virtual environment:

    conda create -n pyfaust_venv python==@PY3_VER@

Add conda-forge channel (it is necessary for the pyfaust dependencies):

    conda config --add channels conda-forge

Load the virtual environment created earlier and install pyfaust:

    conda activate pyfaust_venv
    conda install -c pyfaust pyfaust

\note the -c flag specifies the channel from which to retrieve the pyfaust package.

Try if it works:

    python -c "import pyfaust as pf; print(pf.rand(5,5))"
    Faust size 5x5, density 5, nnz_sum 125, 5 factor(s):
    - FACTOR 0 (double) SPARSE, size 5x5, density 1, nnz 25
    - FACTOR 1 (double) SPARSE, size 5x5, density 1, nnz 25
    - FACTOR 2 (double) SPARSE, size 5x5, density 1, nnz 25
    - FACTOR 3 (double) SPARSE, size 5x5, density 1, nnz 25
    - FACTOR 4 (double) SPARSE, size 5x5, density 1, nnz 25


For further information please rely on the documentation [here](https://docs.conda.io/en/latest/) or about the specific conda-create command please look here: [conda-create](https://docs.conda.io/projects/conda/en/latest/commands/create.html).
