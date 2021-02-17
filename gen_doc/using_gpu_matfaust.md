# Using The GPU FAµST API

__Sections in this tutorial:__  

1. <a href="#1" >Creating a GPU Faust object</a>  
2. <a href="#2" >Generating a GPU Faust</a>  
3. <a href="#3" >Manipulating GPU Fausts and CPU interoperability</a>  
4. <a href="#4" >Benchmarking your GPU with matfaust!</a>  
5. <a href="#5" >Running some FAµST algorithms on GPU</a>  
6. <a href="#6" >Manually loading the pyfaust GPU module</a>

In this tutorial we'll see quickly how to leverage the GPU computing power with matfaust.
Since matfaust 3.0.0 the API has been modified to make the GPU available directly from the Matlab wrapper.
Indeed, a independent GPU module (aka ``gpu_mod``) has been developed for this purpose.

The first question you might ask is: does it work on my computer? Here is the answer: the loading of this module is quite transparent, if an NVIDIA GPU is available and CUDA is properly installed on your system, you have normally nothing to do except installing matfaust to get the GPU implementations at your fingertips. We'll see at the end of this tutorial how to load manually the module and how to get further information in case of an error.

It is worthy to note two drawbacks about the matfaust GPU support:

- Mac OS X is not supported because NVIDIA has stopped to support this OS.
- On Windows and Linux, the matfaust GPU support is currently limited to CUDA 9.2 version.

In addition to these drawbacks, please notice that the GPU module support is still considered in beta status as the code is relatively young and still evolving. However the API shouldn't evolve that much in a near future.

### <a name="1">1. Creating a GPU Faust object</a>

Let's start with some basic Faust creation on the GPU. Almost all the ways of creating a Faust object in CPU memory are also available to create a GPU Faust.
First of all, creating a Faust using the constructor works seamlessly on GPU, the only need is to specify the ``dev`` keyword argument, as follows:

			import matfaust.Faust
			M = rand(10,10);
			N = rand(10, 15);
			gpuF = Faust({M, N}, 'dev', 'gpu')

			Output:
				gpuF =

				- GPU FACTOR 0 (real) DENSE size 10 x 10, addr: 0x154e9a7d8870, density 1.000000, nnz 100
				- GPU FACTOR 1 (real) DENSE size 10 x 15, addr: 0x154e99ccb500, density 1.000000, nnz 150

It's clearly indicated in the output that the Faust object is instantiated in GPU memory (the N and M numpy arrays are copied from the CPU to the GPU memory). However it's also possible to check this programmatically:

			device(gpuF)

			Output:

				ans =

				    'gpu'

While for a CPU Faust you'll get:

			device(Faust({M, N}, 'dev', 'cpu'))

			Output:

			ans =

			    'cpu'

In ``gpuF`` the factors are dense matrices but it's totally possible to instantiate sparse matrices on the GPU as you can do on CPU side.

			S = sprand(10, 15, .25);
			T = sprand(15, 10, .05);
			sparse_gpuF = Faust({S, T}, 'dev', 'gpu') 

			Output:

			sparse_gpuF = 

			- GPU FACTOR 0 (real) SPARSE size 10 x 15, addr: 0x14a5a7e882d0, density 0.220000, nnz 33
			- GPU FACTOR 1 (real) SPARSE size 15 x 10, addr: 0x14a5a7e88910, density 0.053333, nnz 8

You can also create a GPU Faust by explicitly copying a CPU Faust to the GPU memory. Actually, at anytime you can copy a CPU Faust to GPU and conversely. The ``clone()`` member function is here precisely for this purpose. Below we copy ``gpuF`` to CPU and back again to GPU in the new Faust ``gpuF2``.


			cpuF = clone(gpuF, 'cpu')

			Output:

			cpuF = 

			Faust size 10x15, density 1.66667, nnz_sum 250, 2 factor(s): 
			- FACTOR 0 (real) DENSE,  size 10x10, density 1, nnz 100
			- FACTOR 1 (real) DENSE,  size 10x15, density 1, nnz 150

			gpuF2 = clone(cpuF, 'gpu')

			Output:

			gpuF2 = 

			- GPU FACTOR 0 (real) DENSE size 10 x 10, addr: 0x154e43dda130, density 1.000000, nnz 100
			- GPU FACTOR 1 (real) DENSE size 10 x 15, addr: 0x154e43dd7b70, density 1.000000, nnz 150


### <a name="2">2. Generating a GPU Faust</a>

Many of the functions for generating a Faust object on CPU are available on GPU too. It is always the same, you precise the ``dev`` argument by assigning the ``'gpu'`` value and you'll get a GPU Faust instead of a CPU Faust.

For example, the code below will successively create a random GPU Faust, a Hadamard transform GPU Faust, identity GPU Faust and finally a DFT GPU Faust.


			import matfaust.wht
			import matfaust.dft
			% Random GPU Faust
			matfaust.rand(10, 10, 'num_factors', 11, 'dev', 'gpu')

			Output:
			ans =

			- GPU FACTOR 0 (real) SPARSE size 10 x 10, addr: 0x154e43acb150, density 0.500000, nnz 50
			- GPU FACTOR 1 (real) SPARSE size 10 x 10, addr: 0x154e43ddb830, density 0.500000, nnz 50
			- GPU FACTOR 2 (real) SPARSE size 10 x 10, addr: 0x154e43ddbce0, density 0.500000, nnz 50
			- GPU FACTOR 3 (real) SPARSE size 10 x 10, addr: 0x154e43ddce90, density 0.500000, nnz 50
			- GPU FACTOR 4 (real) SPARSE size 10 x 10, addr: 0x154e43de38c0, density 0.500000, nnz 50
			- GPU FACTOR 5 (real) SPARSE size 10 x 10, addr: 0x154e43de4710, density 0.500000, nnz 50
			- GPU FACTOR 6 (real) SPARSE size 10 x 10, addr: 0x154e43de55a0, density 0.500000, nnz 50
			- GPU FACTOR 7 (real) SPARSE size 10 x 10, addr: 0x154e43de6490, density 0.500000, nnz 50
			- GPU FACTOR 8 (real) SPARSE size 10 x 10, addr: 0x154e43de7380, density 0.500000, nnz 50
			- GPU FACTOR 9 (real) SPARSE size 10 x 10, addr: 0x154e43de8270, density 0.500000, nnz 50
			- GPU FACTOR 10 (real) SPARSE size 10 x 10, addr: 0x154e43de9160, density 0.500000, nnz 50
			% Hadamard GPU Faust
			wht(32, 'dev', 'gpu')

			Output:
			ans =

			- GPU FACTOR 0 (real) SPARSE size 32 x 32, addr: 0x154e43df0c20, density 0.062500, nnz 64
			- GPU FACTOR 1 (real) SPARSE size 32 x 32, addr: 0x154e43dee3d0, density 0.062500, nnz 64
			- GPU FACTOR 2 (real) SPARSE size 32 x 32, addr: 0x154e43df3780, density 0.062500, nnz 64
			- GPU FACTOR 3 (real) SPARSE size 32 x 32, addr: 0x154e43df4480, density 0.062500, nnz 64
			- GPU FACTOR 4 (real) SPARSE size 32 x 32, addr: 0x154e43df51c0, density 0.062500, nnz 64
			% Identity GPU Faust
			matfaust.eye(16, 'dev', 'gpu')

			Output:
			ans =

			- GPU FACTOR 0 (real) SPARSE size 16 x 16, addr: 0x154e43de7410, density 0.062500, nnz 16
			% DFT GPU Faust
			dft(32, 'dev', 'gpu')

			Output:
			ans =

			- GPU FACTOR 0 (complex) SPARSE size 32 x 32, addr: 0x154e43df5720, density 0.062500, nnz 64
			- GPU FACTOR 1 (complex) SPARSE size 32 x 32, addr: 0x154e43df21e0, density 0.062500, nnz 64
			- GPU FACTOR 2 (complex) SPARSE size 32 x 32, addr: 0x154e43df3640, density 0.062500, nnz 64
			- GPU FACTOR 3 (complex) SPARSE size 32 x 32, addr: 0x154e43de5520, density 0.062500, nnz 64
			- GPU FACTOR 4 (complex) SPARSE size 32 x 32, addr: 0x154e43de3880, density 0.062500, nnz 64
			- GPU FACTOR 5 (complex) SPARSE size 32 x 32, addr: 0x154e43df2da0, density 0.031250, nnz 32

### <a name="3">3. Manipulating GPU Fausts and CPU interoperability</a>

Once you've created GPU Faust objects, you can perform operations on them staying in GPU world (that is, with no array transfer to CPU memory). That's of course not always possible.

For example, let's consider Faust-scalar multiplication and Faust-matrix product. In the first case the scalar is copied to the GPU memory and likewise in the second case the matrix is copied from CPU to GPU in order to proceed to the computation. However in both cases the Faust factors stay into GPU memory and don't move during the computation.

			% Faust-scalar multiplication
			2*gpuF

			Output:
			ans =

			- GPU FACTOR 0 (real) DENSE size 10 x 10, addr: 0x154e34980120, density 1.000000, nnz 100
			- GPU FACTOR 1 (real) DENSE size 10 x 15, addr: 0x154e99ccb500, density 1.000000, nnz 150

As you see the first factor's address has changed in the result compared to what it was in ``gpuF``. Indeed, when you make a scalar multiplication only one factor is multiplied, the others don't change, they are shared between the Faust being multiplied and the resulting Faust. This is an optimization and to go further in this direction the factor chosen to be multiplied is the smallest in memory (not necessarily the first one).

	% Faust-matrix product (the matrix is copied to GPU
	% then the multiplication is performed on GPU)

			gpuF*rand(size(gpuF,2), 15)

			Output:
			ans =

			   20.0361   18.9839   19.2403   15.6763   17.6467   14.2024   22.0221   16.8871   14.9518   18.0550   16.9907   16.2313   19.1798   17.7156   15.0071
			   17.5964   16.4156   16.3757   13.8127   15.6474   12.0306   18.6480   15.2072   13.2615   14.7832   14.3821   13.4478   15.9126   14.7898   12.3165
			   24.2206   22.4507   22.6955   18.7368   21.6624   16.4971   26.3582   21.3575   18.2987   20.9960   20.2449   19.6180   22.7556   20.8585   17.5761
			   21.9159   21.0406   20.5890   16.8110   19.0772   14.9778   23.8053   19.0746   16.5460   18.4001   18.4172   17.2334   19.9098   18.6816   16.0544
			   24.6674   23.2190   23.0114   18.8039   21.8146   16.8671   27.4397   21.7257   18.6095   22.0368   20.2354   19.1477   22.9723   21.2252   18.0958
			   17.5488   16.5573   16.7642   13.8421   15.7318   12.2374   19.3189   15.1602   13.1840   15.6283   15.1275   14.6058   17.0519   15.7529   13.1074
			   21.8256   20.3676   20.1061   16.9731   19.4686   14.6711   23.4962   18.4808   16.2128   18.6823   18.3134   17.7391   20.4453   19.1173   15.9421
			   19.1147   17.7332   17.3880   14.7597   16.9095   12.8873   20.0457   16.0832   14.1063   15.8157   15.7174   14.8118   17.1714   15.9060   13.5531
			   23.3360   21.6335   21.7350   18.3694   20.7213   16.1601   24.4399   18.9016   17.0501   19.8258   19.6005   18.7066   21.6264   20.0042   16.9067
			   19.1308   17.9028   17.7825   14.9782   16.8418   13.1389   20.1918   16.3237   14.3236   15.9297   15.7401   14.6461   17.1327   15.7691   13.5463

On the contrary, and that matters for optimization, there is no CPU-GPU transfer at all when you create another GPU Faust named for example ``gpuF2`` on the GPU and decide to multiply the two of them like this:

			gpuF2 = matfaust.rand(size(gpuF, 2), 18, 'dev', 'gpu')

			gpuF2 = 

			- GPU FACTOR 0 (real) SPARSE size 15 x 15, addr: 0x154e4bbbcae0, density 0.333333, nnz 75
			- GPU FACTOR 1 (real) SPARSE size 15 x 18, addr: 0x154e34ac9810, density 0.277778, nnz 75
			- GPU FACTOR 2 (real) SPARSE size 18 x 18, addr: 0x154e34b088d0, density 0.277778, nnz 90
			- GPU FACTOR 3 (real) SPARSE size 18 x 18, addr: 0x154e9a8e1df0, density 0.277778, nnz 90
			- GPU FACTOR 4 (real) SPARSE size 18 x 18, addr: 0x154e9af31d30, density 0.277778, nnz 90
			gpuF3 = gpuF*gpuF2

			gpuF3 = 

			- GPU FACTOR 0 (real) DENSE size 10 x 10, addr: 0x154e9a7d8870, density 1.000000, nnz 100
			- GPU FACTOR 1 (real) DENSE size 10 x 15, addr: 0x154e99ccb500, density 1.000000, nnz 150
			- GPU FACTOR 2 (real) SPARSE size 15 x 15, addr: 0x154e4bbbcae0, density 0.333333, nnz 75
			- GPU FACTOR 3 (real) SPARSE size 15 x 18, addr: 0x154e34ac9810, density 0.277778, nnz 75
			- GPU FACTOR 4 (real) SPARSE size 18 x 18, addr: 0x154e34b088d0, density 0.277778, nnz 90
			- GPU FACTOR 5 (real) SPARSE size 18 x 18, addr: 0x154e9a8e1df0, density 0.277778, nnz 90
			- GPU FACTOR 6 (real) SPARSE size 18 x 18, addr: 0x154e9af31d30, density 0.277778, nnz 90

Besides, it's important to note that ``gpuF3`` factors are not duplicated in memory because they already exist for ``gpuF`` and ``gpuF2``, that's an extra optimization: ``gpuF3`` is just a memory view of the factors of ``gpuF`` and ``gpuF2`` (the same GPU arrays are shared between ``Faust`` objects). That works pretty well the same for CPU ``Faust`` objects.

Finally, please notice that CPU Faust objects are not directly interoperable with GPU Fausts objects. You can try, it'll end up with an error.


			cpuF = matfaust.rand(5, 5, 'num_factors', 5, 'dev', 'cpu');
			gpuF = matfaust.rand(5, 5, 'num_factors', 6, 'dev', 'gpu');
			% A first try to multiply a CPU Faust with a GPU one...
			cpuF*gpuF

			Output:
			Error using mexFaustReal
			Handle not valid.

			Error in matfaust.Faust/call_mex (line 3007)
								[varargout{1:nargout}] = mexFaustReal(func_name, F.matrix.objectHandle, varargin{:});

			Error in matfaust.Faust/mtimes_trans (line 714)
										C = matfaust.Faust(F, call_mex(F, 'mul_faust', A.matrix.objectHandle));

			Error in * (line 652)
							G = mtimes_trans(F, A, 0);

			% it doesn't work, you must either convert cpuF to a GPU Faust or gpuF to a CPU Faust before multiplying.
			% A second try using conversion as needed...
			clone(cpuF, 'gpu')*gpuF

			Output:
			ans =

			- GPU FACTOR 0 (real) SPARSE size 5 x 5, addr: 0x154e34aca8b0, density 1.000000, nnz 25
			- GPU FACTOR 1 (real) SPARSE size 5 x 5, addr: 0x154e34aca910, density 1.000000, nnz 25
			- GPU FACTOR 2 (real) SPARSE size 5 x 5, addr: 0x154e34aca7f0, density 1.000000, nnz 25
			- GPU FACTOR 3 (real) SPARSE size 5 x 5, addr: 0x154e34aca850, density 1.000000, nnz 25
			- GPU FACTOR 4 (real) SPARSE size 5 x 5, addr: 0x154e349df930, density 1.000000, nnz 25
			- GPU FACTOR 5 (real) SPARSE size 5 x 5, addr: 0x154e34aca6f0, density 1.000000, nnz 25
			- GPU FACTOR 6 (real) SPARSE size 5 x 5, addr: 0x154e34ae5b80, density 1.000000, nnz 25
			- GPU FACTOR 7 (real) SPARSE size 5 x 5, addr: 0x154e34ae6770, density 1.000000, nnz 25
			- GPU FACTOR 8 (real) SPARSE size 5 x 5, addr: 0x154e349a3d80, density 1.000000, nnz 25
			- GPU FACTOR 9 (real) SPARSE size 5 x 5, addr: 0x154e34adfa30, density 1.000000, nnz 25
			- GPU FACTOR 10 (real) SPARSE size 5 x 5, addr: 0x154e34ae2ba0, density 1.000000, nnz 25
			% Now it works

### <a name="4">4. Benchmarking your GPU with matfaust!</a>

Of course when we run some code on GPU rather than on CPU, it is clearly to enhance the performances. So let's try your GPU and find out if it is worth it or not compared to your CPU.
First, measure how much time it takes on CPU to compute a Faust norm and the dense array corresponding to the product of its factors:

			cpuF = matfaust.rand(1024, 1024, 'num_factors', 10, 'fac_type', 'dense');
			timeit(@() norm(cpuF, 2))

			Output:
			ans =

			    0.1673

			timeit(@() full(cpuF))

			Output:
			ans =

			    2.5657

Now let's make some GPU heat with norms and matrix products!


			gpuF = clone(cpuF, 'dev', 'gpu');
			timeit(@() norm(gpuF, 2))

			Output:
			ans =

			    0.0557

			 timeit(@() full(gpuF))

			Output:
			ans =

			    0.1354

Of course not all GPUs are equal, the results above were obtained with a GTX980 NVIDIA GPU. Below are the results I got using a Tesla V100:

			timeit(@() norm(gpuF, 2))

			Output:
			ans =

			    0.0076

			timeit(@() full(gpuF))

			Output:
			ans =

			    0.0103

Likewise let's compare the performance obtained for a sparse Faust (on a GTX980):

			cpuF2 = matfaust.rand(1024, 1024, 'num_factors', 10, 'fac_type', 'sparse');
			gpuF2 = clone(cpuF2, 'dev', 'gpu');
			timeit(@() norm(gpuF2, 2))

			Output:
			ans =

			    0.0093

			 timeit(@() full(gpuF2))

			Output:
			ans =

			    0.0178

And then on a Tesla V100:

			>> timeit(@() norm(gpuF2, 2))                                                 

			Output:
			ans =

			    0.0059

			>> timeit(@() full(gpuF2))                                                    

			Output:
			ans =

			    0.0102

### <a name="5">5. Running some FAµST algorithms on GPU</a>

Some of the FAµST algorithms implemented in the C++ core are now also available in pure GPU mode.  
For example, let's compare the factorization times taken by the hierarchical factorization when launched on CPU and GPU.  
When running on GPU, the matrix to factorize is copied in GPU memory and almost all operations executed during the algorithm don't imply the CPU in any manner (the only exception at this stage of development is the proximal operators that only run on CPU).  

First please copy the following function in the appropriate filename ``factorize_MEG.m`` into in the current working directory of Matlab (or by adding the destination directory to your path by calling ``addpath``).

			function [MEG16, total_time, err] = factorize_MEG(dev)
				import matfaust.fact.hierarchical
				MEG = matrix.'
				num_facts = 9;
				k = 10;
				s = 8;
				tic
				MEG16 = hierarchical(MEG, {'rectmat', num_facts, k, s}, 'backend', 2020, 'gpu', dev == 'gpu');
				total_time = toc;
				err = norm(MEG16-MEG, 'fro')/norm(MEG, 'fro');
			end

**Warning: THE COMPUTATION CAN LAST THIRTY MINUTES OR SO ON CPU**

So executing this function on a Intel Xeon CPU (E5-260), I got the following results:

			[MEG16, total_time, err] = factorize_MEG('cpu')

			Ouput:
			Faust::hierarchical: 1
			Faust::hierarchical: 2
			Faust::hierarchical: 3
			Faust::hierarchical: 4
			Faust::hierarchical: 5
			Faust::hierarchical: 6
			Faust::hierarchical: 7
			Faust::hierarchical: 8

			MEG16 = 

			Faust size 204x8193, density 0.0631649, nnz_sum 105572, 9 factor(s): 
			- FACTOR 0 (real) SPARSE, size 204x204, density 0.293589, nnz 12218
			- FACTOR 1 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 2 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 3 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 4 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 5 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 6 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 7 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 8 (real) SPARSE, size 204x8193, density 0.0490196, nnz 81930

			total_time =

			   1.4411e+03


			err =

			    0.1289



And for the comparison here are the results I got on Tesla V100 GPU:

			[MEG16, total_time, err] = factorize_MEG('gpu')

			Ouput:
			Faust::hierarchical: 1
			Faust::hierarchical: 2
			Faust::hierarchical: 3
			Faust::hierarchical: 4
			Faust::hierarchical: 5
			Faust::hierarchical: 6
			Faust::hierarchical: 7
			Faust::hierarchical: 8

			MEG16 = 

			Faust size 204x8193, density 0.0631649, nnz_sum 105572, 9 factor(s): 
			- FACTOR 0 (real) SPARSE, size 204x204, density 0.293589, nnz 12218
			- FACTOR 1 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 2 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 3 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 4 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 5 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 6 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 7 (real) SPARSE, size 204x204, density 0.0392157, nnz 1632
			- FACTOR 8 (real) SPARSE, size 204x8193, density 0.0490196, nnz 81930

			total_time =

			  320.0904


			err =

			    0.1296

As you see it's far faster than with the CPU!

### <a name="6">6. Manually loading the pyfaust GPU module</a>

If something goes wrong when trying to use the GPU pyfaust extension, here is how to manually load the module and obtain more information.

The key is the function [enable_gpu_mod](https://faustgrp.gitlabpages.inria.fr/faust/last-doc/html/namespacematfaust.html#a75568ecea590cd9f9cd14dce87bfdc84).  
This function allows to give another try to ``gpu_mod`` loading with the verbose mode enabled.

Below I copy output that show what it should look like when it doesn't work:

			matfaust.enable_gpu_mod('silent', false)

			Output:
			WARNING: you must call enable_gpu_mod() before using GPUModHandler singleton.
			loading libgm
			libgm.so: cannot open shared object file: No such file or directory




**NOTE**: this tutorial was made upon FAµST version 3.0.8.
