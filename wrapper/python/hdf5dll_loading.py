

from sys import platform as _pf
if _pf == 'win32':
    import ctypes
    import sys
    from os.path import join, exists
    for p in sys.path:
        lib_path = join(p, 'pyfaust/lib/hdf5.dll')
        if exists(lib_path):
            ctypes.CDLL(lib_path)
