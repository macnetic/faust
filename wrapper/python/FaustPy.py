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
    """ This class represents a dense matrix by a product of 'sparse' factors (i.e Faust) 
    The aim of the Faust representation is to speed-up multiplication by this matrix
    """

    def  __init__(F,list_factors):
        """ create a Faust from a list of factor.

        Parameter
        ---------
        list_factors : list/tuple of numpy matrices or filepath of the Faust in
        matlab format.
        """
        if(isinstance(list_factors, str)):
            contents = loadmat(list_factors)
            list_factors = contents['faust_factors'][0] 
        F.m_faust = FaustCorePy.FaustCore(list_factors);
        F.m_transpose_flag=0;
        F.shape=F.m_faust.shape(F.m_transpose_flag)


    def getNbRow(F):
        """ return the number of row of the current Faust. """
        return F.shape[0]

    def getNbCol(F):
        """ return the number of column of the current Faust. """
        return F.shape[1]

    def size(F):
        """ Returns the Faust size tuple: getNbRow(), getNbCol().
        """
        return F.shape[0], F.shape[1]

    def transpose(F):
        """ transpose the current Faust. """
        F_trans=copy.copy(F)
        F_trans.m_transpose_flag=not (F.m_transpose_flag)
        F_trans.shape=(F.shape[1],F.shape[0])

        return F_trans;

    def display(F):
        """ display information of the current Faust. """
        print("Struct : ")
        F.m_faust.display(F.m_transpose_flag);


    def __mul__(F,x):
        """ multiplication between the current Faust and a numpy matrix.
        (overload the python operator *, return F*x)

        Parameter
        ---------
        x : 2D numpy ndarray of double scalar, must be FORTRAN contiguous 
        """
        return F.m_faust.multiply(x,F.m_transpose_flag)

    def todense(F):
        """ convert the current Faust into a numpy matrix. 
        return a numpy matrix """


        identity=np.eye(F.getNbCol(),F.getNbCol());
        F_dense=F*identity
        return F_dense


    def __getitem__(F,list_index):
        """ Slicing : Return the value of the current Faust at index list_index.
        (overload of python built-in operator F(indexRow,indexCol)

        Parameter
        ---------
        list_index : tab of length 2 with its elements must be slice, integer or Ellipsis(...) 

        Example of use :
            F[2,3], F[0:dim1,...], F[::-1,::-1] 

        """
        #check if list_index has a 2 index (row and column) 
        if (len(list_index) != 2):
            raise ValueError('list_index must contains 2 elements, the row index and the col index')

        #check if the index are slice or integer or Ellipsis
        for id in list_index:
            if (not isinstance(id, slice)) and (not isinstance(id,int) and (id !=Ellipsis)):
                raise ValueError('list_index must contains slice (1:n) or Ellipsis (...) or int (2)')

        keyCol=list_index[1]
        keyRow=list_index[0]
        identity=np.eye(F.getNbCol(),F.getNbCol())
        if(keyCol != Ellipsis): identity=identity[...,keyCol]
        submatrix=F*identity
        submatrix=submatrix[keyRow,:]
        return submatrix

    def nnz(F):
        """ Gives the number of non-zero elements in the Faust.
        """
        return F.m_faust.nnz()

    def density(F):
        """ Returns the density of the Faust which is the normalized rate of non-zero
        coefficients.
        """
        return float(F.nnz())/(F.shape[0]*F.shape[1])

    def RCG(F):
        """ Returns the Relative Complexity Gain (inverse of the density).
        """
        d = F.density()
        if(d > 0):
            return 1/d
        elif(d == 0):
            return float("inf")
        else:
            return -1

    def norm(F, typenorm=2):
        """
        Returns the 2-norm of The Faust.
        """
        if(typenorm != 2):
            raise Exception("Only 2-norm is supported for the"
                            "Faust")
        return F.m_faust.norm()

    def get_nb_factors(F):
        """
        Returns the Faust's number of factors.
        """
        return F.m_faust.get_nb_factors()

    def get_factor(F, i):
        """
        Returns the Faust's i-th factor as a numpy.ndarray.
        """
        if(F.m_transpose_flag):
            i = F.get_nb_factors()-1-i
        fact = F.m_faust.get_fact(i)
        if(F.m_transpose_flag):
            fact = np.transpose(fact)
        return fact

    def save(F, filename, format="Matlab"):
        """
            Saves the Faust into file.
        """
        if(format != "Matlab"):
            raise Exception("Only Matlab format is supported.")
        mdict = {'faust_factors':
                 np.ndarray(shape=(1, F.get_nb_factors()), dtype=object)}
        for i in range(0, F.get_nb_factors()):
            mdict['faust_factors'][0, i] = F.get_factor(i)
        savemat(filename, mdict)
