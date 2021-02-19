\page FAQ Frequently Asked Question

# Frequently Asked Question

[1. Why did I get a file-not-found error when running demos or examples?](#one)  
[2. How can I launch the integrated unit tests of pyfaust?](#two)  
[3. How to launch the demos with pyfaust?](#three)  

\anchor one
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

\anchor two
## 2. How can I launch the integrated unit tests of pyfaust?

That's actually quite simple, if pyfaust is properly installed, you just have to type this command in a terminal:

	python -c "import pyfaust.tests; pyfaust.tests.run_tests('cpu', 'real')"
	# in some cases, it could be python3 or python3* instead of python

Note: in the above test, the Faust class is tested for the CPU C++ backend and real Faust objects (i.e. dtype == np.float). If you want to test GPU complex Faust, just replace the arguments as this: ``run_tests('gpu', 'complex')``.

\anchor three
## 3. How to launch the demos with pyfaust?

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

\anchor four
## 4. How can I launch the integrated unit tests of matfaust?

TODO

\anchor five
## 5. How to launch the demos with matfaust?

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

