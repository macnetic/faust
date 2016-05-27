#ifndef __FAUST_CONSTRAINT_GENERIC_H__
#define __FAUST_CONSTRAINT_GENERIC_H__
#include <string>
//#include <iostream>
#include "faust_constant.h"


enum faust_constraint_name
{
   CONSTRAINT_NAME_SP, /*!< fixed number of non zero elements, INT (frobenius norm 1) */
   CONSTRAINT_NAME_SPCOL, /*!< fixed number of non zero elements per column INT (frobenius norm 1) */
   CONSTRAINT_NAME_SPLIN, /*!< fixed number of non zero elements per line INT (frobenius norm 1) */
   CONSTRAINT_NAME_NORMCOL,/*!< 2nd norm of the colons of A REAL */
   CONSTRAINT_NAME_SPLINCOL,
   CONSTRAINT_NAME_CONST, /**< Matrix equal to A ; MAT */
   CONSTRAINT_NAME_SP_POS,/**< fixed number of non zeros coefficients, INT (frobenius norm 1) */
   CONSTRAINT_NAME_BLKDIAG,
   CONSTRAINT_NAME_SUPP, /**< Matrix which support is equal to A ; MAT ; (frobenius norm 1)*/
   CONSTRAINT_NAME_NORMLIN,/**< 2nd norm of the lines of matrix A ; REAL  */
};

template<typename FPP,Device DEVICE> class MatDense;

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust {

    // modif AL AL
    template<typename FPP,Device DEVICE>
    class MatDense;

/*!
\brief Contains the generic constraint parameters for the hierarchical factorization. See following table for more precision about the type of constraint. <br>
 <img src="../../doc/html/constraint.png" alt="constraint parameters" width=800px />
*/
    template<typename FPP,Device DEVICE>
    class ConstraintGeneric
    {
        public:
        ConstraintGeneric() : constraint_name(CONSTRAINT_NAME_SP),nb_rows(32),nb_cols(32) {} // contrainte par defaut (a voir avec Luc)

        ConstraintGeneric(
            const faust_constraint_name& constraint_name_,
            const faust_unsigned_int nb_rows_,
            const faust_unsigned_int nb_cols_) :
                constraint_name(constraint_name_),
                nb_rows(nb_rows_),
                nb_cols(nb_cols_){}

        ConstraintGeneric(const ConstraintGeneric& constraint) :
            constraint_name(constraint.constraint_name),
            nb_rows(constraint.nb_rows),
            nb_cols(constraint.nb_cols){}



        const char* getType() const;
        const char* get_constraint_name()const;
        const faust_constraint_name getConstraintType() const;
        bool isConstraintParameterInt()const;
        bool isConstraintParameterReal()const;
        bool isConstraintParameterMat()const;


        const faust_unsigned_int getRows() const {return nb_rows;}
        const faust_unsigned_int getCols() const {return nb_cols;}

        virtual void set_default_parameter()=0;
        virtual void check_constraint_name()const=0;
        virtual void project(MatDense<FPP,DEVICE> & mat)const=0;

        ~ConstraintGeneric(){};

        protected:
        /// type of constraint
        const faust_constraint_name constraint_name;
        // parameter of constraint
        //const parameter_type parameter;
        const faust_unsigned_int nb_rows;
        const faust_unsigned_int nb_cols;

        private :
        static const char * class_name;


    };

}

/////// functions useful for parsing config_file (xml,matio,txt_file) //////
bool isConstraintNameInt(const char * type);

bool isConstraintNameReal(const char * type);

bool isConstraintNameMat(const char * type);

faust_constraint_name getEquivalentConstraint(const char * type);
//! \fn getTypeConstraint
//! \brief Information about the type of the constraint
//! \param type
//! \return 0 if type correspond to an integer constraint
//! \return 1 if type correspond to a FPP (Floating Point Precision) constraint
//! \return 2 if type correspond to a Dense Matrix constraint
//! \return else throw an error
int getTypeConstraint(const char * type);

#include "faust_ConstraintGeneric.hpp"

#endif
