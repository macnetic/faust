
from sys import platform as _pf
from glob import glob
blasomp_loading_verbose = True
internal_blasomp_loading = True

blasomp_name_noext = 'libopenblaso'


def load_lib_blaso():
    """
        Loads openblaso if it is embedded in the pyfaust package.
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
        lib_path = join(p, 'pyfaust/lib/'+blasomp_name_noext+'.'+ext)
        lib_path = glob(lib_path+"*")[0]
        if exists(lib_path):
            ctypes.CDLL(lib_path)
            if blasomp_loading_verbose:
                print(lib_path+" has been loaded successfully.")
            break
    else:
        print("failed to load " + blasomp_name_noext +
              " (didn't find any shared lib in"
              " sys.path:", sys.path)


# load liblaso pyfaust embedded library if found in pyfaust location
if _pf in ['linux']:
    load_lib_blaso()
