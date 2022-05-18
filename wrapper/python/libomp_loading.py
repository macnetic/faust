
from sys import platform as _pf
libomp_loading_verbose = False
internal_libomp_loading = True


def load_lib_omp():
    """
        Loads libomp if it is embedded in the pyfaust package.
    """
    from os.path import join, exists
    import sys
    import ctypes
    import ctypes.util
    if _pf == 'darwin':
        ext = 'dylib'
    elif _pf == 'linux':
        ext = 'so'
    elif _pf == 'win32':
        ext = 'dll'
    for p in sys.path:
        lib_path = join(p, 'pyfaust/lib/libomp.'+ext)
        if exists(lib_path):
            ctypes.CDLL(lib_path)
            if libomp_loading_verbose:
                print(lib_path+" has been loaded successfully.")
            break
    else:
        print("failed to load libomp (didn't find any shared lib in"
              " sys.path:", sys.path)


def inform_user_how_to_install_libomp():
    """
    """
    from os.path import exists
    if _pf == 'darwin':
        if not exists('/opt/local/lib/libomp/libomp.dylib'):
            print("""ERROR: OpenMP for clang is not properly installed on your system.
                  You need to install it through MacPorts.
                  1) Install MacPorts: https://www.macports.org/install.php
                  2) Install libomp:
                      sudo port install libomp-devel libomp
                      sudo port -f activate libomp
                  """)
    elif _pf == 'linux':
        if not exists('/usr/lib64/libomp.so'):
            print("""ERROR: OpenMP for clang is not properly installed on your system.
                  You need to install the proper packages, for example:
                  - On Fedora the package is libomp-11*,
                  - On Debian the package is libomp-11*.
                  """)

def try_modify_wrapper_lib_on_macos():
    from os.path import dirname, join
    from os import system
    from glob import glob
    from sys import platform
    if platform != 'darwin':
        return
    wrapper_path = glob(join(dirname(dirname(__file__)), '_Faust*so'))[0]
    libomp_path = join(dirname(__file__), 'lib/libomp.dylib')
    system('install_name_tool -change /opt/local/lib/libomp/libomp.dylib '+libomp_path+' '+wrapper_path)


# load libomp pyfaust embedded library if found in pyfaust location
# and if the code runs on macOS
if _pf in ['darwin', 'linux']:
    try_modify_wrapper_lib_on_macos()
    if internal_libomp_loading:
        load_lib_omp()
    else:
        inform_user_how_to_install_libomp()
