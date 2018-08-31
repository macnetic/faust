##############################################################################
##                              Description:                                ##
##                                                                          ##
##        Cython Header making the links between C++ class FaustCpp         ##
##        and Cython environment, works like C/C++ Header                   ##
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

from libcpp cimport bool

cdef extern from "FaustCoreCpp.h" :
    cdef cppclass FaustCoreCpp[FPP] :
        FaustCoreCpp()
        void Display() const
        const char* to_string() const
        void push_back(FPP* valueMat,unsigned int nbrow,unsigned int nbcol)
        void push_back(FPP* data, int* row_ptr, int* id_col, int nnz, int nrows,
                       int ncols)
        void multiply(FPP* value_y,int nbrow_y,int nbcol_y,FPP* value_x,
                      int nbrow_x, int nbcol_x);#,bool isTranspose*/);
        unsigned int getNbRow() const
        unsigned int getNbCol() const
#        void setOp(const bool isTransposed,unsigned int& nbRowOp, unsigned int& nbColOp)const;
        unsigned long long nnz() const
        double norm(int ord) const
        double normFro() const
        double get_nb_factors() const
        unsigned int get_fact_nb_rows(unsigned int& i) const
        unsigned int get_fact_nb_cols(unsigned int& i) const
        double get_fact(unsigned int& i, FPP* fact_ptr) const
        FaustCoreCpp[FPP]* get_slice(unsigned int, unsigned int, unsigned int,
                                    unsigned int) const
        FaustCoreCpp[FPP]* slice(unsigned int, unsigned int, unsigned int,
                                    unsigned int)
        void save_mat_file(const char* filepath) const
        FaustCoreCpp[FPP]* transpose()
        FaustCoreCpp[FPP]* conjugate()
        FaustCoreCpp[FPP]* adjoint()
        FaustCoreCpp[FPP]* mul_faust(FaustCoreCpp[FPP]*)
        @staticmethod
        FaustCoreCpp[FPP]* randFaust(unsigned int t,
                                             unsigned int min_num_factors, unsigned int max_num_factors,
                                             unsigned int min_dim_size,
                                             unsigned int max_dim_size, float density)

cdef extern from "FaustFact.h":
    cdef cppclass PyxConstraintGeneric:
        int name
        unsigned long num_rows
        unsigned long num_cols

    cdef cppclass PyxConstraintInt(PyxConstraintGeneric):
        unsigned long parameter

    cdef cppclass PyxConstraintScalar[FPP](PyxConstraintGeneric):
        FPP parameter

    cdef cppclass PyxConstraintMat[FPP](PyxConstraintGeneric):
        FPP* parameter # shape = num_rows, num_cols

    cdef cppclass PyxStoppingCriterion[FPP]:
        bool is_criterion_error
        int num_its
        FPP error_treshold
        unsigned long max_num_its

    cdef cppclass PyxParamsFact[FPP]:
        int num_facts
        bool is_update_way_R2L
        FPP init_lambda
        FPP step_size
        PyxConstraintGeneric** constraints # (num_facts-1)*2 elts
        unsigned int num_constraints
        bool is_verbose
        bool constant_step_size

    cdef cppclass PyxParamsFactPalm4MSA[FPP](PyxParamsFact[FPP]):
        FPP** init_facts # num_facts elts
        unsigned long* init_fact_sizes
        PyxStoppingCriterion stop_crit

    cdef cppclass PyxParamsHierarchicalFact[FPP](PyxParamsFact[FPP]):
        unsigned int num_rows
        unsigned int num_cols
        PyxStoppingCriterion* stop_crits #must be of size 2
        bool is_fact_side_left

    cdef FaustCoreCpp[FPP]* fact_palm4MSA[FPP](FPP*,unsigned int, unsigned int,
                                             PyxParamsFactPalm4MSA[FPP]*)

    cdef FaustCoreCpp[FPP]* fact_hierarchical[FPP](FPP*,unsigned int, unsigned int,
                                             PyxParamsHierarchicalFact[FPP]*)

