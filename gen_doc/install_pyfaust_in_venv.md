\page install_pyfaust_in_venv Installing pyfaust in Python virtual environments

Python @PY3_VER@
================

First, you need to create your virtual environment into a directory reserved to pyfaust:

	$ python3 -m venv test_pyfaust-3.6

Second, you need to enable the virtual environment:

	$ source ./test_pyfaust-3.6/bin/activate

\note For Windows users the command is rather:

	C:\> call .\test_pyfaust-3.6\Scripts\activate

Then, if all worked correctly you must see the test_pyfaut-3.6 prompt of your virtual environment.
It remains now to install pyfaust in this environment.

	$ pip install pyfaust-*.whl

All the dependencies will be downloaded and installed by the command above.

You can launch a test python command to validate the installation:

	$ python -c 'from pyfaust import version; print(version())'

As as result, you should see the version of pyfaust installed if all went properly.

Python 2.7.15+
=============

You need at first to install virtualenv:

	$ pip install virtualenv

NOTE: take care to use the good version of pip (corresponding to the python 2.7 version you want use).

Then you can create your python 2.7 virtual environment:

	$ virtualenv test_pyfaust-2.7

The next is all the same as for python 3, you enter your virtual env.:

	$ source ./test_pyfaust-2.7/bin/activate

\note For Windows users the command is rather:

	C:\> call .\test_pyfaust-3.6\Scripts\activate

Install the pip wheel package:

	$ pip install pyfaust-*.whl

And finally launch the version test command to check the whole process:

	$ python -c 'from pyfaust import version; print(version())'



