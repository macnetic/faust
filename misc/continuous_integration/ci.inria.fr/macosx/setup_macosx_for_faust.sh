#!/bin/bash

[[ ! $(sudo id) = 'uid=0('* ]] && echo "error: this script must run as root" && exit 1

echo "This script helps to setup a MacOS X 10.9 Maverick host to build FAµST (pytfaust and matfaust wrappers packages)"
echo "However, a certain amount of settings can't be automatized."
echo "The script must be launched from a ci.inria.fr MacOSX host VM."


###### 1. Install and configure Mac Ports
cd /Users/ci/Downloads
wget https://github.com/macports/macports-base/releases/download/v2.6.4/MacPorts-2.6.4-10.9-Mavericks.pkg
installer -pkg MacPorts-2.6.4-10.9-Mavericks.pkg -target /
export PATH=/opt/local/bin:$PATH
echo 'export PATH=/opt/local/bin:$PATH' >> /Users/ci/.bash_profile
##### 2. Install and configure gitlab-runner
curl --output /usr/local/bin/gitlab-runner https://gitlab-runner-downloads.s3.amazonaws.com/v11.9.2/binaries/gitlab-runner-darwin-amd64
chmod +x /usr/local/bin/gitlab-runner
echo 'export PATH=/usr/local/bin:$PATH' >> /Users/ci/.bash_profile # not even necessary
echo -n "enter gitlab FAµST token (https://gitlab.inria.fr/faustgrp/faust/-/settings/ci_cd):"
read TOKEN
su ci -c "echo -e 'https://gitlab.inria.fr\n'$TOKEN'\n'macosx-$(date +%m-%d-%Y)'\nmacos,matlab\nshell' | /usr/local/bin/gitlab-runner register" 
su ci -c "gitlab-runner install" # to run as ci # installed here: /Users/ci/Library/LaunchAgents/gitlab-runner.plist
su ci -c "gitlab-runner start"

### 3. Install dependency libraries to build FAµST
yes | port install eigen3 matio cmake p7zip zlib
sudo ln -sf /opt/local/lib/libmatio.a /usr/local/lib/
sudo ln -sf /opt/local/include/matio.h /usr/local/include
sudo ln -sf /opt/local/include/matio_pub.h /usr/local/include
yes | port install libomp-devel libomp 
### 4. Install and configure clang compiler environment (and OpenMP)
port -f activate libomp
sudo cp /opt/local/lib/libomp/libomp.dylib /opt/local/lib/libomp/libomp_macports.dylib # needed for ci package jobs
yes | port install clang-8.0
echo 'export OpenMP_INC_DIR=/opt/local/include/libomp' >> /Users/ci/.bash_profile
echo 'export OpenMP_gomp_LIBRARY=/opt/local/lib/libomp/libgomp.dylib' >> /Users/ci/.bash_profile
mv /usr/bin/clang /usr/bin/clang_
mv /usr/bin/clang++ /usr/bin/clang++_
ln -sf /opt/local/bin/clang-mp-8.0 /opt/local/bin/clang
ln -sf /opt/local/bin/clang++-mp-8.0 /opt/local/bin/clang++
### 5. Install Python packages and dependencies
yes | port install graphviz doxygen
sudo ln -sf /opt/local/bin/dot /usr/bin/
yes | port install python39 py39-pip
yes | port install jpeg # pillow (indirect denpendency of pyfaust needs this to build)
port select --set python python39
port select --set pip pip39
ln -sf /opt/local/bin/python3.9 /opt/local/bin/python3
yes | port install py39-cython
yes | port select --set cython cython39
yes |pip install doxypypy chardet wheel pygsp 'numpy==1.22'
# get dest dir before downloading because otherwise it would be the current path
DOXYPYPY_DIR=$(dirname $(python3 -c "import doxypypy; print(doxypypy.__file__)"))
wget https://raw.githubusercontent.com/Feneric/doxypypy/master/doxypypy/doxypypy.py
# use rather the github doxypypy instead of the pypy version
sudo mv doxypypy.py ${DOXYPYPY_DIR}
yes | pip install matplotlib # separately because it could fail 
#### 6. Install C++ libtorch
wget https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.4.0.zip
unzip libtorch-macos-1.4.0.zip
mv libtorch /opt/local/
echo "Please manually install matlab by copying directly the directory from /Volumes/Untitled/ attached to faust2-macos-2019 (advice: use scp)"
echo "First, a DATA storage volume must be created and attached to the VM instance (use command: diskutil mount)."
echo "Finally add matlab to the PATH (in .bash_profile)."
echo "It is also useful to symlink the matlab binary in /usr/bin/ for the root user env. to be OK for installing/testing the pkg (postinstall script needs to find matlab)."
echo "Other need about matlab OpenMP: you need to copy the library as for example: ciosx:~ ci$ cp /Volumes/Untitled/MATLAB_R2018b.app/sys/os/maci64/libiomp5.dylib /opt/local/lib/libomp/libiomp5_matlab.dylib"
echo "ABOUT SUDO: add this line in /etc/sudoers: ci ALL=(ALL:ALL) NOPASSWD: ALL (this way the runner won't need to type the password for running commands as root)"
echo "matio should be reinstalled too because the version from MacPorts is 1.5.21, not compatible with FAµST (at least because of the type mat_sparse_t which contain fields of different types compared to previous version). The version to install is matio-1.5.19 in /usr/local or another way is to copy files /usr/local/matio_pub.h /usr/local/matio.h and /usr/local/lib/libmatio*from faust2-macosx-2019
echo ===== manually install concurrent version of python, eg.:
echo port install py310 py310-cython py310-pip
echo pip-3.10 install doxypypy chardet wheel pygsp 'numpy==1.22'
echo DOXYPYPY_DIR=$(dirname $(python3.10 -c "import doxypypy; print(doxypypy.__file__)"))
echo wget https://raw.githubusercontent.com/Feneric/doxypypy/master/doxypypy/doxypypy.py
echo sudo mv doxypypy.py ${DOXYPYPY_DIR}
#echo "WARNING: the matlab installer link is likely to fail"
#wget "https://esd.mathworks.com/R2018b/R2018b/installers/web/matlab_R2018b_maci64.dmg.zip?__gda__=1608282997_72bccb4cda14022857f8691914aea21a&dl_id=aYudAe6N&ext=.zip" || echo "the link failed, please update."
#mv matlab_R2018b_maci64.dmg.zip\?__gda__\=1608282997_72bccb4cda14022857f8691914aea21a\&dl_id\=aYudAe6N\&ext\=.zip matlab_R2018b_maci64.dmg.zip
#unzip matlab_R2018b_maci64.dmg.zip
#hdiutil attach matlab_R2018b_maci64.dmg


