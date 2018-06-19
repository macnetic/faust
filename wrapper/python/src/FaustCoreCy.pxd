##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython Header making the links between C++ class FaustCpp         ##
##        and Cython environment, works like C/C++ Header                   ##
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

from libcpp cimport bool

cdef extern from "FaustCoreCpp.h" :
    cdef cppclass FaustCoreCpp[FPP] :
        FaustCoreCpp();
        void Display() const;
        void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol);
        void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,
                      int nbrow_x, int nbcol_x);#,bool isTranspose*/);
        unsigned int getNbRow() const;
        unsigned int getNbCol() const;
#        void setOp(const bool isTransposed,unsigned int& nbRowOp, unsigned int& nbColOp)const;
        unsigned long long nnz() const;
        double norm(int ord) const;
        double normFro() const;
        double get_nb_factors() const;
        unsigned int get_fact_nb_rows(unsigned int& i) const;
        unsigned int get_fact_nb_cols(unsigned int& i) const;
        double get_fact(unsigned int& i, FPP* fact_ptr) const;
        void save_mat_file(const char* filepath) const;
        FaustCoreCpp[FPP]* transpose()
        FaustCoreCpp[FPP]* conjugate()
        FaustCoreCpp[FPP]* adjoint()
