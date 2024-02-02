# Quick Guide to Install and Use the FAµST Wrappers

[1. Prerequisites](#prerequisites)

[2. Installation](#installation)

[3. Testing and Troubleshooting](#installation_testing)

[4. Usage](#usage)

\anchor prerequisites 
## 1. Prerequisites

### 1.1 For Python

@note  This section treats only of system packages/installers but note that some PIP and Anaconda packages are also available. Installing those packages is recommended, in particular because of the easier weak dependency management. You might install them in virtual environments:
\ref install_pyfaust_in_venv.<br/> 

FAµST is designed for the Python ecosystem which as usual relies on the numpy and scipy packages.

If you want to run the FAµST demos (module pyfaust.demo), you'll also need the matplotlib package.

Please ensure these packages are installed on your system. One way to proceed is to use the pip program delivered with your python version.

	pip install numpy scipy matplotlib # pygsp
	# it could be pip3 instead of pip

Note that you probably have to install those packages for all versions of Python you want to use (each one has normally its associated pip executable).<br/>
Note also that you can rely on your system package manager to install the Python packages (e.g. dnf/yum repositories on Fedora/Centos Linux systems).

FAµST is compiled for Python @PY3_VER@.

If you want to use FAµST with Python 3 you must use precisely the @PY3_VER@ version because the FAµST Python wrapper delivered within the binary package is compiled for that minor version only.

@note pygsp is an optional python package to install in order to generate graphs and their Laplacians for testing the FGFT/eigen decomposition algorithms. It's not mandatory but can be used in a pyfaust related [notebook](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/Faust_factorization.html).

### 1.2 For Matlab

The FAµST wrapper for Matlab has been tested on several Matlab versions starting from version 2016a,b up to version 2022b.<br/>
However it's not totally excluded that FAµST works with older and moreover newer versions.

\anchor installation 
## 2. Installation

@OS_SPECIFIC_INSTALL_GUIDE_INSTRUCS@

\anchor installation_testing 
## 3. Testing and Troubleshooting

Normally, after installing, nothing is left to do. The installation process should have seamlessly set up the Faust wrappers for Python and Matlab into your environment.
Nevertheless, it could be useful to check that it really worked and set the environment manually if needed like explained below.

### 3.1 Testing the Matlab Wrapper

To test whether the FAµST Matlab wrapper auto-setup succeeded at install stage, you can open a terminal and type:

	matlab -nodisplay -nojvm -r "matfaust.rand(10, 10, 'num_factors', 1, 'density', .5, 'fac_type', 'dense', 'field', 'complex')"

Note: if Matlab is not set in your PATH environment variable you need to replace `matlab' with its full path
	(e.g. on macOS /Applications/Matlab/MATLAB_R2018b.app/bin/matlab)

It works only if you see an output similar to:

	Faust size 10x10, density 0.50, nnz_sum 50, 1 factor(s):
	- FACTOR 0 (complex) DENSE, size 10x10, density 0.5, nnz 50
	% other values are possible for density, etc. because of the random generation

Otherwise it didn't work. So here is how to setup the wrapper manually.

First, launch a Matlab terminal, then go in the FAµST directory:

	>> cd @FAUST_INSTALL_PATH@/matlab

Then launch the script that is responsible to add FAµST location in your Matlab path.

	>> setup_FAUST
	>> % then test again FAµST
	>> matfaust.rand(10, 10, 'num_factors', 1, 'dim_sizes', 10, 'density', .5, 'fac_type', 'dense', 'field', 'complex')

### 3.1.1 Permament path using savepath

Once matfaust path has been loaded into matlab environment you can save the
configuration for future uses by simply typing:

    >> savepath

The next time you'll run Matlab the matfaust package will be directly
available.

### 3.1.2 Permanent path using Matlab's startup.m

Another way for that change of matlab path to be applied permanently, is to automatize the `addpath()' call made by ``setup_FAUST``.<br/>
For that purpose:

<ol>
<li>Look into your matlab installation directory (e.g. on macOS /Applications/Matlab/MATLAB_R2018b.app/bin/matlab), next:
<li>Go into the sub-directory toolbox/local
<li>Edit the file startup.m by adding the following lines:
<pre>
    oldpwd=pwd;
    cd @FAUST_INSTALL_PATH@/matlab; % for Windows users this is: C:\Program Files\faust\matlab
    setup_FAUST
    cd oldpwd;
</pre>
</ol>

Finally type <a href="https://fr.mathworks.com/help/matlab/ref/rehash.html#d120e1067468">`rehash toolbox'</a> in your current matlab terminal and restart Matlab in order to verify the configuration is permanent.

Note: you can also edit a user's specific startup.m script instead of system's startup.m. Look the Matlab documentation [here](https://fr.mathworks.com/help/matlab/matlab_env/startup-options.html).

OK! You can follow the [quick start usage](#usage) now.

### 3.2 Testing the Python Wrapper

To test whether the FaµST Python wrapper has been setup properly, simply open a terminal and type:

	python3 -c "import pyfaust; print('It works.')"

	# it could be python@PY3_VER@ or just python,
	# depending on your configuration and the python version you want to use

It goes without saying that if the wrapper is set up properly you'll see the message "It works." as a result of the command above.

On the contrary, the following error message

	Traceback (most recent call last):
	  File "<string>", line 1, in <module>
	ModuleNotFoundError: No module named 'pyfaust'

indicates that you need to add the Python wrapper manually in your Python path as demonstrated below.

- For Linux and macOS in a Bash terminal:

	$ export PYTHONPATH=$PYTHONPATH:@FAUST_INSTALL_PATH@/python3
	# and test again
	$ python3 -c "import pyfaust; print('It works.')"

- For Windows in the command prompt:

	set PYTHONPATH=%PYTHONPATH%;@FAUST_INSTALL_PATH@/python"
	:: and test again
	python3 -c "import pyfaust; print('It works.')"

If it fails again you are likely on a different version of Python or the auto-setup script failed for any reason during installation. Return to the [Prerequisites](#prerequisites) section and check your Python environment matches FAµST requirements.
Otherwise it works but you'll need to set the `export' command manually in one of your startup scripts to have it set for once and for all (e.g. in ~/.bashrc if you are a Bash user).
<br/>On Windows you can set the PYTHONPATH environment variable in the configuration panel (system, advanced settings).

Finally, note that you may get another error message indicating for example that either numpy or scipy is not available.

	ImportError: No module named numpy

For that issue look at [Prerequisites](#prerequisites).

OK! You can follow the [quick start usage](#usage) now.

\anchor usage 
## 4. Usage

### 4.1 Matlab Wrapper

Let's test FAµST with the quickstart script, from a matlab terminal type:

	>> import matfaust.demo.quickstart.*
	>> A = quick_start

And then call some functions on Faust object A to test a bit of the FAµST API:

	>> rcg(A)
	>> density(A)
	>> numfactors(A)

@note if you're wondering what are these functions just consult the inline doc:

	>> help A.rcg

Retrieve the product factors:

	>> F1 = factors(A,1);
	>> F2 = factors(A,2);

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

Note that the doxygen documentation for the Matlab API is also available locally after installing Faust (in the installation folder). You can consult it online: [namespacematfaust.html](@API_DOC_BASE_URL@/html/namespacematfaust.html).

### 4.2 Python Wrapper

In the same spirit as the Matlab tutorial showed in the previous section, you can execute the quick start script for Python.

	$ python3 -c "from pyfaust.demo import quickstart; quickstart.run()"
	dimension of the Faust :  (1000, 2000)
	multiplication SPEED-UP using Faust
	Faust is 1.83845941093 faster than a full matrix
	Faust nnz: 300000
	Faust density: 0.15
	Faust RCG: 6.66666666667
	Faust norm: 55156456.373
	Faust nb of factors: 2
	end quickstart.py


You can also go through the Python terminal to build a FAµST product and call its object methods.

	$ python3
	>>> import pyfaust
	>>> A = pyfaust.Faust(filepath='A.mat') # A is the FAµST created through quickstart script
	>>> A.rcg()
	6.666666666666667
	>>> A.density()
	0.15
	>>> A.numfactors()
	2
	>>> F1 = A.factors(0)
	>>> F2 = A.factors(1)
	>>> A.shape
	(1000, 2000)

Note that the doxygen documentation for the Python API is also available locally after installing Faust (in the installation folder). You can consult it online: [namespacepyfaust.html](@API_DOC_BASE_URL@/html/namespacepyfaust.html).



