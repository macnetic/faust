##############################################################################
##                              Description:                                ##
##                                                                          ##
##          FaustPy is a python module  which delivered                     ##
##          a class named Faust which represents a dense matrix             ##
##          by a product of 'sparse' factors (i.e Faust)                    ##
##          Python wrapper class implemented in C++                         ##
##                                                                          ##
##  For more information on the FAuST Project, please visit the website     ##
##  of the project : <http://faust.gforge.inria.fr>                         ##
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


import copy

import numpy as np
from scipy.io import savemat, loadmat
import FaustCorePy


class Faust:
    """ This class represents a dense matrix by a product of 'sparse' factors
    (i.e Faust).
    The aim of the Faust representation is to speed-up multiplication by this
    matrix.
    """

    def  __init__(F, list_factors=None, filepath=None, core_obj=None):
        """ Creates a Faust from a list of factors or a file.

        Args:
            list_factors : list/tuple of numpy matrices or filepath of the Faust in
        matlab binary format (.mat).
        """
        if(core_obj):
            F.m_faust = core_obj
        else:
            if(filepath and isinstance(filepath, str)):
                    contents = loadmat(filepath)
                    list_factors = contents['faust_factors'][0]
            if(list_factors is not None):
                F.m_faust = FaustCorePy.FaustCore(list_factors);
            #else:
                #TODO: manage empty Faust


    def get_nb_rows(F):
        """
        Gives the Faust's number of rows.

        Returns:
                The Faust number of rows.
        """
        return F.m_faust.shape()[0]

    def get_nb_cols(F):
        """
        Gives the Faust's number of columns.

        Returns:
                The Faust number of columns.
        """
        return F.m_faust.shape()[1]

    def size(F):
        """
        Gives the size of the Faust.

        Returns:
            The Faust size tuple: get_nb_rows(), get_nb_cols().
        """
        #return F.shape[0], F.shape[1]
        return F.m_faust.shape()

    def transpose(F):
        """
        Transposes the current Faust.
        """
        F_trans = Faust(core_obj=F.m_faust.transpose())
        return F_trans

    def conj(F):
        """
        Returns the conjugate of the current Faust.
        """
        F_conj = Faust(core_obj=F.m_faust.conj())
        return F_conj

    def getH(F):
        """
        Returns the conjugate transpose of the current Faust.
        """
        F_ctrans = Faust(core_obj=F.m_faust.getH())
        return F_ctrans

    def display(F):
        """
        Displays information describing the current Faust.
        """
        F.m_faust.display()

    def __mul__(F, M):
        """
        Multiplies the Faust by the numpy matrix M.

        This method overloads the Python operator *.

        Args:
            M: 2D numpy ndarray of double scalar, must be Fortran contiguous
            (i.e. Column-major order; `order' argument passed to
            np.ndararray() must be equal 'F').

        Returns:
            The result of the multiplication as a numpy matrix.

        Raises:
            ValueError


        Examples:
            >>> import numpy as np
            >>> x = np.random.randint(120, size=(F.get_nb_cols(), 1))
            >>> y = F*x
        """
        return F.m_faust.multiply(M)

    def todense(F):
        """
        Converts the current Faust into a numpy matrix.

        Returns:
            A numpy matrix.
        """
        identity = np.eye(F.get_nb_cols(), F.get_nb_cols())
        F_dense = F*identity
        return F_dense

    def __getitem__(F, indices):
        """
        Returns a coefficient of the Faust by its indices.

        This function overloads the Python built-in.

        Args:
            indices: array of length 2, its elements must be of kind: slice, integer or
            Ellipsis (...).

        Returns:
            The float element requested.

        Examples:
            >>> F[2, 3]
            >>> F[0:dim1, ...]
            >>> F[::-1, ::-1]
        """
        # check if indices has a 2 index (row and column)
        if (len(indices) != 2):
            raise ValueError('indices must contains 2 elements, '
                             'the row index and the col index')

        # check if the index are slice or integer or Ellipsis
        for id in indices:
            if (not isinstance(id, slice) and
                not isinstance(id, int) and
                id != Ellipsis):
                raise ValueError('indices must contains slice (1:n)'
                                 ' or Ellipsis (...) or int (2)')

        keyCol = indices[1]
        keyRow = indices[0]
        identity = np.eye(F.get_nb_cols(), F.get_nb_cols())
        if(keyCol != Ellipsis):
            identity = identity[..., keyCol]
        submatrix = F*identity
        submatrix = submatrix[keyRow, :]
        return submatrix

    def nnz(F):
        """
        Gives the number of non-zero elements in the Faust.

        The function sums together the number of non-zeros elements of
        each factor and returns the result. Note that in fact the sum is
        computed at Faust creation time and kept in cache.

        Returns:
            The number of non-zero elements as an integer.

        Examples:
            >>> F.nnz()
        """
        return F.m_faust.nnz()

    def density(F):
        """
        Calculates the normalized rate of non-zero coefficients.

        Returns:
            The density value as a float.

        Examples:
            >>> F.density()
        """
        #return float(F.nnz())/(F.shape[0]*F.shape[1])
        return float(F.nnz()/(F.get_nb_cols()*F.get_nb_rows()))

    def RCG(F):
        """
        Computes the Relative Complexity Gain (inverse of the density).

        Returns:
            The RCG value as a float.
            If the density is zero it will be float("inf").
            If the density is negative it will be -1.

        Examples:
            >>> F.RCG()
        """
        d = F.density()
        if(d > 0):
            return 1/d
        elif(d == 0):
            return float("inf")
        else:
            return float(-1)

    def norm(F, ord=2):
        """
        Computes the norm of the Faust.

        Args:
            ord: (Optional) The type of norm to return. For now, it's
            limited to L1 and L2 norms only. Set it to 1 for L1, 2 for L2
            (default value).

        Returns:
            The norm as a float.

        Raises:
            ValueError.


        Examples:
            >>> F.norm()
            >>> F.norm(2)
        """
        if(ord not in [1, 2]):
            raise ValueError("Only L1 and L2 norms are supported (values 1 or"
                             "2 expected.")
        return F.m_faust.norm(int(ord))

    def get_nb_factors(F):
        """
        Gives the Faust's number of factors.

        Returns:
            The number of factors.

        Examples:
            >>> F.get_nb_factors()
        """
        return F.m_faust.get_nb_factors()

    def get_factor(F, i):
        """
        Returns the Faust's i-th factor as a numpy.ndarray.

        Args:
            i: The integer index of the factor to get back.

        Returns:
            The i-th factor as a dense matrix (of type numpy.ndarray).

        Raises:
            ValueError.


        Examples:
            >>> F.get_factor(0)
        """
        fact = F.m_faust.get_fact(i)
        return fact

    def save(F, filepath, format="Matlab"):
        """
        Saves the Faust into file.

        Args:
            filepath: The path for saving the Faust (should ends with .mat
            if Matlab format used). It can be either an absolute or relative path.
            format: (Optional) It designates the format to use for
            writing. By default, it's "Matlab" to save the Faust in a .mat
            file (and it's currently the only format available).

        Raises:
            ValueError.


        Examples:
            >>> F.save("myFaust.mat")
        """
        if(format not in ["Matlab"]):
            raise ValueError("Only Matlab or Matlab_core format is supported.")
        if(format == "Matlab"):
            F.m_faust.save_mat_file(filepath)

    def isReal(F):
        """
        Returns True if F is a real Faust and False if it's a complex Faust.
        """
        return F.m_faust.isReal()
