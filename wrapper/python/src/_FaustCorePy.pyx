##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython source file making the links between :                     ##
##        the Cython module FaustCoreCy and Python                          ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.inria.fr>                         ##
##                                                                          ##
##                              License:                                    ##
##  Copyright (2016):   Nicolas Bellot, Adrien Leman, Thomas Gautrais,      ##
##                      Luc Le Magoarou, Remi Gribonval                     ##
##                      INRIA Rennes, FRANCE                                ##
##                      http://www.inria.fr/                                ##
##                                                                          ##
##  The FAuST Toolbox is distributed under the terms of the GNU Affero      ##
##  General Public License.                                                 ##
##  This program is free software: you can redistribute it and/or modify    ##
##  it under the terms of the GNU Affero General Public License as          ##
##  published by the Free Software Foundation.                              ##
##                                                                          ##
##  This program is distributed in the hope that it will be useful, but     ##
##  WITHOUT ANY WARRANTY; without even the implied warranty of              ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    ##
##  See the GNU Affero General Public License for more details.             ##
##                                                                          ##
##  You should have received a copy of the GNU Affero General Public        ##
##  License along with this program.                                        ##
##  If not, see <http://www.gnu.org/licenses/>.                             ##
##                                                                          ##
##                             Contacts:                                    ##
##      Nicolas Bellot  : nicolas.bellot@inria.fr                           ##
##      Adrien Leman    : adrien.leman@inria.fr                             ##
##      Thomas Gautrais : thomas.gautrais@inria.fr                          ##
##      Luc Le Magoarou : luc.le-magoarou@inria.fr                          ##
##      Remi Gribonval  : remi.gribonval@inria.fr                           ##
##                                                                          ##
##############################################################################

cimport FaustCoreCy # cython C++ wrapper
#TODO: reorganize imports
from FaustCoreCy cimport PyxConstraintGeneric, PyxConstraintInt, \
PyxConstraintMat, PyxConstraintScalar, PyxStoppingCriterion, \
PyxParamsFactPalm4MSA, PyxMHTPParams

import numpy as np
cimport numpy as np
import copy

from libc.stdlib cimport malloc, free
from libc.stdio cimport printf
from libc.string cimport memcpy, strlen
from libcpp cimport bool, complex
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free
from scipy import sparse
from scipy.sparse import csr_matrix, csc_matrix
from re import match
import sys, os, pyfaust

cdef check_matrix(isReal, M, message=""):
        if not isinstance(M, (np.ndarray) ):
            raise ValueError(message+'input must be a numpy ndarray')
        if(isReal):
            if M.dtype not in ['double', 'float32']:
                raise ValueError(message+'input numpy array dtype must be double or float32.')
        else:
#            M=M.astype(complex,'F')
            if(M.dtype not in ['complex', 'complex128'] ): #could fail if complex128 etc.
                raise ValueError('input array must be complex(128) array')
        #TODO: raise exception if not real nor complex
        if not M.flags['F_CONTIGUOUS']:
            raise ValueError(message+'input array must be Fortran contiguous (Colmajor)')

        ndim_M=M.ndim;

        if (ndim_M > 2) | (ndim_M < 1):
            raise ValueError(message+'input matrix/array invalid number of dimensions')

def _is_gpu_mod_enabled():
    return FaustCoreCy._is_gpu_mod_enabled()

def enable_gpu_mod(libpaths=None, backend='cuda', silent=False, fatal=False):
    cdef char * c_libpath
    cdef void* gm_handle
    if libpaths == None:
        # use default paths
        pyfaust_path = os.path.dirname(pyfaust.__file__)
        # don't use os.path.join or os.path.sep because anyway
        # lib name suffix and prefix depend on OS
        libpaths = []
        for gpu_backend in ['-cu11.4', '-cu9.2', '']: # the last one is
        # for descendant/backward compatibility
        # (the backend are sorted per order of preference)
            if sys.platform == 'linux':
                libpaths += [pyfaust_path+"/lib/libgm"+gpu_backend+".so"]
            elif sys.platform == 'darwin':
                libpaths += [pyfaust_path+"/lib/libgm"+gpu_backend+".dylib"]
            elif sys.platform == 'win32':
                libpaths += [pyfaust_path+"\lib\gm"+gpu_backend+".dll"]
    for libpath in libpaths:
        blibpath = libpath.encode('utf-8')
        c_libpath = blibpath
        #c_libpath = libpath
        gm_handle = FaustCoreCy._enable_gpu_mod(c_libpath, silent);
        if gm_handle != NULL:
            # the backend is successfully loaded, end the loading loop
            break
        if gm_handle == NULL:
            if fatal and libpaths[-1] == libpath:
                raise Exception("Can't load gpu_mod library, tried all the candidate paths ("
                                +libpaths+"). Maybe they are not"
                                " correct or the backends (CUDA) are not installed or"
                                " configured properly so"
                                " the libraries are not found.")

import sys, os, pyfaust
# tries to load the libgm library silently,
# if not enabled at build time it will do nothing
enable_gpu_mod(silent=False)
