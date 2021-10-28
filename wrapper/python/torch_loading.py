# this module loads libtorch libraries FAµST needs after downloading and
# installing libtorch if necessary
from os.path import dirname, basename, join
from pathlib import Path
import sys
import ctypes
import ctypes.util

if sys.platform == 'win32':
    # TODO: verify lib names
    TORCH_LIB = 'torch'
    C10_LIB = 'c10'
    URL_LIB = 'https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-1.4.0.zip'
    # TODO: preload gomp ?
elif  sys.platform == 'linux':
    TORCH_LIB = 'libtorch.so'
    C10_LIB = 'libc10.so'
    URL_LIB = 'https://download.pytorch.org/libtorch/cpu/libtorch-cxx11-abi-shared-with-deps-1.4.0%2Bcpu.zip'
    # pre-load libgomp from the system to avoid the libtorch's gomp version to
    # be load later
    ctypes.CDLL(ctypes.util.find_library('gomp'))
else: #macOS
    # no-need to preload gomp, not used for macOS
    TORCH_LIB = 'libtorch.dylib'
    C10_LIB = 'libc10.dylib'
    URL_LIB = 'https://download.pytorch.org/libtorch/cpu/libtorch-macos-1.4.0.zip'

BASE_URL = dirname(URL_LIB)
ARCH_NAME = basename(URL_LIB)
DOWNLOAD_PATH = Path.home()

from .datadl import download_uncompress
# download/uncompress libraries (if not here already)
download_uncompress(uncompress_dir=DOWNLOAD_PATH, base_url=BASE_URL,
                    arch_name=ARCH_NAME, data_name='libtorch archive',
                    already_downloaded_msg=False,
                    extra_file_to_check_dl='libtorch')

if sys.platform == 'darwin':
    from os import system
    # modify the rpath of libtorch in order for the linker to find the libiomp5 to
    # load from libtorch lib directory (so we're not using OpenMP from MacPorts
    # here)
    system('install_name_tool -add_rpath "/Users/$USER/libtorch/lib" /Users/$USER/libtorch/lib/libtorch.dylib 2>/dev/null')

# load shared libraries in memory
for libname in [C10_LIB, TORCH_LIB]: #, GOMP_LIB]:
    print("loading", libname)
    ctypes.CDLL(join(DOWNLOAD_PATH,'libtorch','lib', libname),
                mode=ctypes.RTLD_LOCAL)#, mode=ctypes.RTLD_GLOBAL)
