\page install_pyfaust_in_venv Installing pyfaust in Python virtual environments

Python @PY3_VER@
================

First, you need to create your virtual environment into a directory reserved to pyfaust:

	$ python3 -m venv test_pyfaust-@PY3_VER@

Second, you need to enable the virtual environment:

	$ source ./test_pyfaust-@PY3_VER@/bin/activate

\note For Windows users the command is rather:

	C:\> call .\test_pyfaust-@PY3_VER@\Scripts\activate

Then, if all worked correctly you must see the test_pyfaut-@PY3_VER@ prompt of your virtual environment.

It remains now to install pyfaust in this environment.

	$ pip install pyfaust-*.whl

All the dependencies will be downloaded and installed by the command above.

You can launch this python command to validate the installation:

	$ python -c 'from pyfaust import version; print(version())'

As a result, you should see the version of pyfaust installed if all went properly.

\note In the virtual environment python is necessarily python version 3 (because the creation was made through python3 above).

Python 2.7.15+
=============

You need at first to install virtualenv:

	$ pip install virtualenv

\note Take care to use the good version of pip (corresponding to the python 2.7 version you want to use).

Then you can create your python 2.7 virtual environment:

	$ virtualenv test_pyfaust-2.7

The next is all the same as for python 3, you enter the virtual env.:

	$ source ./test_pyfaust-2.7/bin/activate

\note For Windows users the command is rather:

	C:\> call .\test_pyfaust-3.6\Scripts\activate

Install the pip wheel package:

	$ pip install pyfaust-*.whl

And finally launch the version test command to check the whole process:

	$ python -c 'from pyfaust import version; print(version())'


Anaconda
========

Within Anaconda you can also create virtual environments and use pyfaust into it. It's quite similar to the Python way described in previous sections. For that purpose you'll use the command conda create.

Please rely on the documentation here: https://docs.conda.io/en/latest/

Or for more details about this specific command please look here: [conda-create](https://docs.conda.io/projects/conda/en/latest/commands/create.html)

Please note that if you're using conda, you still need to install pyfaust (whl) package through pip. You can install pip with conda and then install pyfaust pip package as described above with a pip install command (in a conda virtual environment or not). And final point, don't forget to use the good version of python3: @PY3_VER@.x.
