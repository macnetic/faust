# Quick Guide to Install and Use FAµST Wrappers

[1. Prerequisites](#prerequisites)

[2. Installation](#installation)

[3. Installation Testing](#installation_testing)

[4. Usage](#usage)

## <a name="prerequisites">1. Prerequisites<a/>

### 1.1 For Python

FAµST supports at least Python 2.7 and is also compiled for Python @PY3_VER@.

If you want to use FAµST with Python 3 you must use this exact version @PY3_VER@ because the FAµST Python wrapper delivered within the binary package is compiled for that version only.

FAµST is designed for the Python ecosystem which as usual relies on numpy and scipy packages. So you have to ensure these packages are installed on your system.

### 1.2 For Matlab

The FAµST wrapper for Matlab has been tested on several Matlab versions starting from version 2016 up to version 2018. However it's not totally excluded that FAµST works with older versions.

## <a name="installation"> 2. Installation</a>

@OS_SPECIFIC_INSTALL_GUIDE_INSTRUCS@

## <a name="installation_testing">3. Installation Testing<a/>

Normally, after installing, nothing is left to do. The installation process should have seamlessly set up the Faust wrappers for Python and Matlab into your environment.
Neverthless, it could be useful to check that it really worked and set the environment manually if needed like explained below.


### 3.1 Testing Python Wrapper Auto-Setup

Here is how to test the FaµST Python wrapper has been setup properly from a Bash terminal:

		$ python2 -c "import pyfaust" && echo "It works"

Obviously, if the wrapper is set up properly you'll see the message "It works" as a result from the command above.

If you see the following error message, then it didn't work:

		Traceback (most recent call last):
		  File "<string>", line 1, in <module>
		ModuleNotFoundError: No module named 'pyfaust'

Likewise, for python3 you can use the previous command for setup checking (currently on macOS only FAµST for Python 2 is supported).

If the auto-setup has failed you need to add the Python wrapper manually in your Python path, here is how to proceed:

		$ export PYTHONPATH=@CMAKE_INSTALL_PREFIX@/python
		# and test again
		$ python2 -c "import pyfaust" && echo "It works"

If it fails again you are likely on a different version of Python. Returns to the Prerequisites section and check your Python versions match FAµST requirements.

### 3.2 Testing Matlab Wrapper Auto-Setup

		$ matlab -nodisplay -nojvm -r "import matfaust.Faust;F = Faust.rand(1, 10, .5, 'dense', false);exit"

It works only if you see an output similar to:

		F =

		Faust size 10x10, density 0.56999, nnz_sum 57, 1 factor(s):
		- FACTOR 0 (complex) DENSE, size 10x10, density 0.57, nnz 57

In other case it didn't work. So here is how to setup the wrapper manually :

First, you launch Matlab terminal, then you go in the FAµST directory :

		>> cd @CMAKE_INSTALL_PREFIX@/matlab

And finally you launch the script that is responsible to add FAµST path in your Matlab path.

		>> setup_FAUST

OK! You can follow the quick start usage now.


## <a name="usage">4. Usage<a/>

### 4.1 Matlab Wrapper

Let's test FAµST with the quickstart script:

		>> quick_start

And then type further instructions to test a bit of FAµST API:

Let's see what variables quickstart script has added.

		>> whos

Now call some functions on FAµST object A:

		>> rcg(A)
		>> density(A)
		>> get_num_factors(A)

Retrieve the product factors (as full matrices):

		>> F1 = get_fact(A,1);
		>> F2 = get_fact(A,2);

Check the sizes:

		>> size(F1)
		>> size(F2)
		>> size(A)

For more information on FAµST API, and whole function listing, consult the doc:

		>> doc('Faust')

N.B.: to access the documentation, you need to be in the graphical interface of matlab.

You should rather consult the doxygen doc for the API in [@CMAKE_INSTALL_PREFIX@/doc/html/index.html](file://@CMAKE_INSTALL_PREFIX@/doc/html/index.html) from your web browser.


### 4.2 Python Wrapper

In the same idea, you can execute the quick start script for Python.

		$ python2.7 @CMAKE_INSTALL_PREFIX@/python/quickstart.py

You can also go through the python terminal to build a FAµST product and call its object methods.

		$ python2.7
		>>> import pyfaust
		>>> import pyfaust
		>>> A = pyfaust.Faust('@CMAKE_INSTALL_PREFIX@/python/A.mat') # A is the FAµST created through quickstart script
		>>> A.rcg()
		6.666666666666667
		>>> A.density()
		0.15
		>>> A.get_num_factors()
		2
		>>> F1 = A.get_factor(0)
		>>> F2 = A.get_factor(1)
		>>> A.shape
		(1000, 2000)


