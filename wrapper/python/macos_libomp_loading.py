
from sys import platform as _pf

macos_libomp_loading_verbose = True


def load_lib_omp():
    """
        Loads libomp if it is embedded in the pyfaust package.
    """
    from os.path import join, exists
    import sys
    import ctypes
    import ctypes.util
    for p in sys.path:
        lib_path = join(p, 'pyfaust/lib/libomp.dylib')
        if exists(lib_path):
            ctypes.CDLL(lib_path)
            if macos_libomp_loading_verbose:
                print(lib_path+" has been loaded successfully.")
            break
    else:
        print("failed to load libomp (didn't find any shared lib in"
              " sys.path:", sys.path)


def inform_user_how_to_install_libomp():
    """
    """
    from os.path import exists
    if not exists('/opt/local/lib/libomp/libomp.dylib'):
        print("""ERROR: OpenMP is not properly installed on your system.
              You need to install it through MacPorts.
              1) Install MacPorts: https://www.macports.org/install.php
              2) Install libomp:
                  sudo port install libomp-devel libomp
                  sudo port -f activate libomp
              """)


# load libomp pyfaust embedded library if found in pyfaust location
# and if the code runs on macOS
if _pf == 'darwin':
    # load_lib_omp()
    inform_user_how_to_install_libomp()
