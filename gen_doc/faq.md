\page FAQ Frequently Asked Question

# Frequently Asked Question

[1. why did I get a file-not-found error when running demos or examples ?](#one)

\anchor one
## 1. why did I get a file-not-found error when running demos or examples ?


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

