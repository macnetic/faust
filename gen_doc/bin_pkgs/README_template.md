# Quick Guide to Install and Use the FAµST Wrappers

[1. Prerequisites](#prerequisites)

[2. Installation](#installation)

[3. Testing and Troubleshooting](#installation_testing)

[4. Usage](#usage)

## <a name="prerequisites">1. Prerequisites<a/>

### 1.1 For Python

FAµST is designed for the Python ecosystem which as usual relies on the numpy and scipy packages.

Please ensure these packages are installed on your system. One way to install them is to use the pip program delivered with your python version.

	pip install numpy scipy
	# it could be pip3 instead of pip

Note that you probably have to install those packages for all versions of Python you want to use (each one have normally its associated pip executable).<br/>
Note also that you can rely on your system package manager to install the Python packages (e.g. dnf/yum repositories on Fedora/Centos Linux systems).

FAµST supports at least Python 2.7 and is also compiled for Python @PY3_VER@.

If you want to use FAµST with Python 3 you must use precisely the @PY3_VER@ version because the FAµST Python wrapper delivered within the binary package is compiled for that version only.


### 1.2 For Matlab

The FAµST wrapper for Matlab has been tested on several Matlab versions starting from version 2016a,b up to version 2018a,b.<br/>
However it's not totally excluded that FAµST works with older or newer versions.

## <a name="installation"> 2. Installation</a>

@OS_SPECIFIC_INSTALL_GUIDE_INSTRUCS@

## <a name="installation_testing">3. Testing and Troubleshooting<a/>

Normally, after installing, nothing is left to do. The installation process should have seamlessly set up the Faust wrappers for Python and Matlab into your environment.
Neverthless, it could be useful to check that it really worked and set the environment manually if needed like explained below.

### 3.1 Testing the Matlab Wrapper

To test whether the FAµST Matlab wrapper auto-setup succeeded at install stage, you can open a terminal and type:

	matlab -nodisplay -nojvm -r "import matfaust.FaustFactory;F = FaustFactory.rand(1, 10, .5, 'dense', 'complex');disp(F);exit"

Note: if Matlab is not set in your PATH environment variable you need to replace `matlab' with its full path
	(e.g. on macOS /Applications/Matlab/MATLAB_R2018b.app/bin/matlab)

It works only if you see an output similar to:

	Faust size 10x10, density 0.57, nnz_sum 57, 1 factor(s):
	- FACTOR 0 (complex) DENSE, size 10x10, density 0.57, nnz 57
	>> % other values are possible for density, etc. because of the random generation

Otherwise it didn't work. So here is how to setup the wrapper manually.

First, launch a Matlab terminal, then go in the FAµST directory:

	>> cd @CMAKE_INSTALL_PREFIX@/matlab

Then launch the script that is responsible to add FAµST location in your Matlab path.

	>> setup_FAUST
	>> % then test again FAµST
	>> import matfaust.FaustFactory;F = FaustFactory.rand(1, 10, .5, 'dense', 'complex');disp(F);exit

For that change to be applied permanently, you need to automatize the `addpath()' call made by setup_FAUST.<br/>
For that purpose:

<ol>
<li>Look into your matlab installation directory (e.g. on macOS /Applications/Matlab/MATLAB_R2018b.app/bin/matlab), next:
<li>Go into the sub-directory toolbox/local
<li>Edit the file startup.m by adding the follwing line:
<pre>
	addpath(genpath('@CMAKE_INSTALL_PREFIX@/matlab'))
</pre>
</ol>
OK! You can follow the [quick start usage](#usage) now.

### 3.2 Testing the Python Wrapper

To test whether the FaµST Python wrapper has been setup properly, simply open a terminal and type:

	python2 -c "import pyfaust; print('It works.')"

	# it could be python2.7, python3, python@PY3_VER@ or just python,
	# depending your configuration and the python version you want to use

It goes without saying that if the wrapper is set up properly you'll see the message "It works." as a result of the command above.

On the contrary, the following error message

	Traceback (most recent call last):
	  File "<string>", line 1, in <module>
	ModuleNotFoundError: No module named 'pyfaust'

indicates that you need to add the Python wrapper manually in your Python path as demonstrated below.

- For Linux and macOS in a Bash terminal:

	$ export PYTHONPATH=$PYTHONPATH:@CMAKE_INSTALL_PREFIX@/python
	# and test again
	$ python2 -c "import pyfaust; print('It works.')"

- For Windows in the command prompt:

	set PYTHONPATH=%PYTHONPATH%;@CMAKE_INSTALL_PREFIX@/python
	:: and test again
	python -c "import pyfaust; print('It works.')"

If it fails again you are likely on a different version of Python or the auto-setup script failed for any reason during installation. Return to the [Prerequisites](#prerequisites) section and check your Python environment matches FAµST requirements.
Otherwise it works but you'll need to set the `export' command manually in one of your startup scripts to have it set for once and for all (e.g. in ~/.bashrc if you are a Bash user).
<br/>On Windows you can set the PYTHONPATH environment variable in the configuration panel (system, advanced settings).

Finally, note that you may get another error message indicating for example that either numpy or scipy is not available.

	ImportError: No module named numpy

For that issue look at [Prerequisites](#prerequisites).

OK! You can follow the [quick start usage](#usage) now.


## <a name="usage">4. Usage<a/>

### 4.1 Matlab Wrapper

Let's test FAµST with the quickstart script, from a matlab terminal type:

	>> import matfaust.demo.quickstart.*
	>> quick_start

And then type further instructions to test a bit of the FAµST API:

Let's see what variables quickstart script has added.

	>> whos

Now call some functions on FAµST object A:

	>> rcg(A)
	>> density(A)
	>> get_num_factors(A)

Retrieve the product factors:

	>> F1 = get_factor(A,1);
	>> F2 = get_factor(A,2);

Check the sizes:

	>> size(F1)
	>> size(F2)
	>> size(A)

Check storage organization of the factors:

	>> issparse(F1)
	>> issparse(F2)

For more information on the FAµST API, and a whole function listing, consult the doc:

	>> doc('matfaust.Faust')
	>> % or:
	>> import matfaust.*;doc('Faust')

N.B.: to access the documentation, you need to be in the graphical interface of Matlab.

Note that the doxygen documentation for the Matlab API is also available locally after installing Faust. You can consult it from your web browser: [@CMAKE_INSTALL_PREFIX@/doc/html/namespacematfaust.html](file://@CMAKE_INSTALL_PREFIX@/doc/html/namespacematfaust.html).

### 4.2 Python Wrapper

In the same spirit than the Matlab tutorial showed in the previous section, you can execute the quick start script for Python.

	$ python2.7 @CMAKE_INSTALL_PREFIX@/python/quickstart.py

You can also go through the Python terminal to build a FAµST product and call its object methods.

	$ python2.7
	>>> import pyfaust
	>>> A = pyfaust.Faust(filepath='@CMAKE_INSTALL_PREFIX@/python/A.mat') # A is the FAµST created through quickstart script
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

Note that the doxygen documentation for the Python API is also available locally after installing Faust. You can consult it from your web browser: [@CMAKE_INSTALL_PREFIX@/doc/html/namespacepyfaust.html](file://@CMAKE_INSTALL_PREFIX@/doc/html/namespacepyfaust.html).


