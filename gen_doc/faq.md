\page FAQ Frequently Asked Question

# Frequently Asked Question



**About pyfaust:**  
[1. How can I launch the integrated unit tests of pyfaust?](#py_one)  
[2. How to launch the demos with pyfaust?](#py_two)  
[3. How to run the PALM4MSA algorithm in a step-by-step fashion?](#py_three)  
[4. Why does PALM4MSA fail on Windows when the numpy array is complex?](#py_four)

**About matfaust:**  
[1. Why did I get a file-not-found error when running demos or examples?](#mat_one)  
[2. How to fix mex not found error: "Undefined function or variable 'mexFaustReal'"?](#mat_two)  
[3. How to launch the demos with matfaust?](#mat_three)  
[4. How can I launch the integrated unit tests of matfaust?](#mat_four)  
[5. How to run the PALM4MSA algorithm in a step-by-step fashion?](#mat_five)  

# About matfaust

\anchor mat_one
## 1. Why did I get a file-not-found error when running demos or examples?


For example, running this matlab command if quickstart.mat file is not found you'll get the following error:

	$ matlab -nojvm -nodisplay -r "import matfaust.demo.quickstart; quickstart.quick_start()"

						 < M A T L A B (R) >
				       Copyright 1984-2017 The MathWorks, Inc.
				       R2017a (9.2.0.556344) 64-bit (glnxa64)
						   March 27, 2017

	For online documentation, see http://www.mathworks.com/support
	For product information, visit www.mathworks.com.

	Error using load
	Unable to read file 'faust_quick_start.mat'. No such file or directory.

	Error in matfaust.Faust (line 226)
					load(filename);

	Error in matfaust.demo.quickstart.quick_start (line 19)
				A=Faust('faust_quick_start.mat')

The same kind of error might happen also with pyfaust, the python wrapper, which depends on the same data.

Normally, at installation stage the FAµST externalized data (basically a bunch of matlab .mat files) is downloaded from a remote web server and unarchived in FAµST installation path.
Nevertheless, it might not work properly for any reason (e.g. network issue happening at installation), so here are two ways to download the data manually.

#### 1.1. The Easy Way:

Just reinstall FAµST! It will relaunch the data download if your network connection is properly enabled. If it doesn't work, repeat the operation after having deleted the data folder located in FAµST installation path or python site-packages/pyfaust folder (read the 1.2 below to determine this path).

#### 1.2. The Manual Way:

This is assumed that you installed the pyfaust wrapper from a pip package or through one of the installers/packages.

Here is an example of commands you can type to download the data (it's on Linux bash but similar if not totally the same commands apply to other systems) :


First, find where pyfaust is installed (if you installed pyfaust for python 3, run python3 not 2.7):

	$ python -c "import pyfaust; print(pyfaust.__file__)"
	/usr/lib64/python2.7/site-packages/pyfaust/__init__.pyc

Second, run these commands to download and uncompress the data:

- If you installed pyfaust from a pip package (NOTE: the path is not the same if you use python3, look above to get the prefix path to which append the data folder):

	$ DATA_DEST=/usr/lib64/python2.7/site-packages/pyfaust/data; sudo rm -Rf $DATA_DEST; mkdir $DATA_DEST; python /usr/lib64/python2.7/site-packages/pyfaust/datadl.py $DATA_DEST
	Downloading FAµST data: 100 %
	====== data downloaded: /tmp/faust_data-2.4.2.zip
	Uncompressing zip archive to /usr/lib64/python2.7/site-packages/pyfaust/data

- Or if you installed it another way (e.g with a rpm or .exe installer):

	# get the matlab wrapper path:
	$ matlab -nojvm -r "import matfaust.Faust;which Faust"
	/opt/local/faust/matlab/+matfaust/@Faust/Faust.m  % matfaust.Faust constructor
	# the result command indicates we have to download the data in the matfaust wrapper path here: /opt/local/faust/matlab/data
	$ DATA_DEST=/opt/local/faust/matlab/data
	$ sudo rm -Rf $DATA_DEST; mkdir $DATA_DEST; python /usr/lib64/python2.7/site-packages/pyfaust/datadl.py $DATA_DEST
        Downloading FAµST data: 100 %
        ====== data downloaded: /tmp/faust_data-2.4.2.zip
        Uncompressing zip archive to /opt/local/faust/matlab/data

If finally you still don't manage to retrieve the data, please write an [email](index.html) with all the details (at least the version of FAµST, the installer used and of course your system).

\anchor mat_two
## 2. How to fix mex not found error: "Undefined function or variable 'mexFaustReal'"?

If something went wrong, for example at install stage, it is possible that Matlab doesn't find FAµST (in particular the mex or even matlab files `.m`).  
In this case, an error similar to two examples might be raised:

	>> import matfaust.rand         
	Error using import
	Import argument 'matfaust.rand' cannot be found or cannot be imported.

	>> import matfaust.rand
	>> rand(5,4)
	Undefined function or variable 'mexFaustReal'.

	Error in matfaust.rand (line 320)
				core_obj = mexFaustReal('rand', num_rows, num_cols,
				fac_type, min_num_factors, max_num_factors,
				min_dim_size, max_dim_size, density, per_row);

To fix this issue, you have to update the Matlab path manually, however a script (named setup_FAUST.m) is here to help.

	% go in the matlab installation path : /opt/local/faust/matlab on Linux and MacOS X, C:\Program Files\Faust\matlab on Windows
	>> cd /path/to/faust/matlab
	% then run the setup_FAUST.m script
	>> setup_FAUST
	Welcome to the Matlab wrapper of the FAuST C++ toolbox.FAuST root directory is /path/to/faust/matlab/
	Adding path /path/to/faust/matlab/ and subdirectories to matlab path

	 To get started with the FAuST Toolbox : 

	 launch quick_start or run_all_demo.m 

\anchor mat_three
## 3. How to launch the demos with matfaust?

You might run all the demos at once or one by one. In the former case you should run this matlab code:

	>> import matfaust.demo.runall
	>> runall()

Note: the raw result are then in output files located in the directory ``./output``.

To launch the demo one by one, you can pick the command in your interest below:


For the [BSL](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/namespacematfaust_1_1demo.html) demo:

	>> import matfaust.demo.bsl.*
	>> BSL()
	>> Fig_BSL()

For the [DFT](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classmatfaust_1_1demo_1_1fft.html) demo:

	>> import matfaust.demo.fft.speed_up_fourier
	>> speed_up_fourier

For the [Hadamard](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classmatfaust_1_1demo_1_1hadamard.html) demo:

	>> import matfaust.demo.hadamard
	>> hadamard.speed_up_hadamard()

For the [runtime comparison](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classmatfaust_1_1demo_1_1runtimecmp.html) demo:

	>> import matfaust.demo.runtimecmp
	>> runtimecmp.runtime_comparison
	>> runtimecmp.Fig_runtime_comparison

And for the last and simplest demo, the [quickstart script](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classmatfaust_1_1demo_1_1quickstart.html):

	>> import matfaust.demo.quickstart
	>> quickstart.quick_start()
	>> quickstart.factorize_matrix()
	>> quickstart.construct_Faust_from_factors()

\anchor mat_four
## 4. How can I launch the integrated unit tests of matfaust?

TODO (in the meanwhile you can see the pyfaust entry)

\anchor mat_five

## 5. How to run the PALM4MSA algorithm in a step-by-step fashion?

TODO (in the meanwhile you can see the pyfaust entry)

# About pyfaust
\anchor py_one
## 1. How can I launch the integrated unit tests of pyfaust?

That's actually quite simple, if pyfaust is properly installed, you just have to type this command in a terminal:

	python -c "import pyfaust.tests; pyfaust.tests.run_tests('cpu', 'real')"
	# in some cases, it could be python3 or python3* instead of python

Note: in the above test, the Faust class is tested for the CPU C++ backend and real Faust objects (i.e. dtype == np.float). If you want to test GPU complex Faust, just replace the arguments as this: ``run_tests('gpu', 'complex')``.

\anchor py_two
## 2. How to launch the demos with pyfaust?

You might run all the demos at once or one by one. In the former case you should run this python code:

	>>> from pyfaust.demo import runall
	>>> runall()

Note: the raw result are then in output files located in the directory ``pyfaust.demo.DEFT_RESULTS_DIR``.

To generate the figures that go by default into ``pyfaust.demo.DEFT_FIG_DIR``, you have to type in a python terminal or insert in a python script the following instructions:

	>>> from pyfaust.demo import allfigs
	>>> allfigs()

To launch the demo one by one, you can pick the command in your interest below:

For the [BSL](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classpyfaust_1_1demo_1_1bsl.html) demo:

	>>> from pyfaust.demo import bsl
	>>> bsl.run() # runs the experiment
	>>> bsl.fig() # generates the figures

For the [DFT](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classpyfaust_1_1demo_1_1fft.html) demo:

	>>> from pyfaust.demo import fft
	>>> fft.speed_up_fourier()
	>>> fft.fig_speedup_fourier()

For the [Hadamard](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classpyfaust_1_1demo_1_1hadamard.html) demo:

	>>> from pyfaust.demo import hadamard
	>>> hadamard.run_fact()
	>>> hadamard.run_speedup_hadamard()
	>>> hadamard.run_norm_hadamard()
	>>> hadamard.figs()

For the [runtime comparison](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classpyfaust_1_1demo_1_1runtimecmp.html) demo:

	>>> from pyfaust.demo import runtimecmp
	>>> runtimecmp.run()
	>>> runtimecmp.fig()

And for the last and simplest demo, the [quickstart script](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/classpyfaust_1_1demo_1_1quickstart.html):

	>>> from pyfaust.demo import quickstart
	>>> quickstart.run()

\anchor py_three

## 3. How to run the PALM4MSA algorithm in step-by-step fashion?

Although the verbose mode of the PALM4MSA implementation allows displaying some info, it might be useful in order to analayze the algorithm (e.g. build the loss function or check the matrix supports evolution), to be able to run just one iteration at a time, get all the Faut layers, the scale factor (lambda), do whatever one needs with them, and then continue to the next iteration.

It implies to reinitialize the next iteration in the same state it was at the end of the past iteration. The script [step-by-step_palm4msa.py](./step-by-step_palm4msa.py) shows how to do that for a matrix factorization in two factors but it is not much different with a greater number of factors.
On the script end PALM4MSA is performed all iterations at once in order to verify the step-by-step execution was consistent.

Below is an example of output you should obtain running the script (into which you can see that the at-once and iteration-by-iteration executions match perfectly):

	python3 step-by-step_palm4msa.py | tail -3
	Relative error when running all PALM4MSA iterations at once:  0.2978799226115671
	Last relative error obtained when running PALM4MSA iteration-by-iteration:  0.2978799226115671
	Relative error comparing the final Fausts obtained either by in step-by-step PALM4MSA versus all-iterations-at-once PALM4MSA:  2.1117031467008879e-16


\anchor py_four

## 4. Why does PALM4MSA fail on Windows when the numpy array is complex?

If you run the following Python code on Windows with pyfaust, it will fail with the next error (while it works properly on Linux or Mac OS X).

```
from pyfaust.fact import palm4msa
from pyfaust.factparams import ParamsPalm4MSA, ConstraintList, StoppingCriterion
import numpy as np
from numpy.random import rand
d = 3
cons = ConstraintList('sp', d, d, d, 'sp', d, d, d)
stop_crit = StoppingCriterion(num_its=10)
param = ParamsPalm4MSA(cons, stop_crit)
M = rand(d, d)+1j*rand(d,d) # M is complex
print(M.dtype)
F1 = palm4msa(M, param, ret_lambda=True)
```

**Error**: TypeError: data type 'float128' not understood.

Unfortunately when built with Microsoft Visual Studio the ``float128`` numpy is not supported and somehow this type is used when the numpy scalar type (dtype) ``numpy.complex128`` is used for a matrix. Until numpy or VS is fixed to support this missing type pyfaust PALM4MSA won't run on complex matrices.

