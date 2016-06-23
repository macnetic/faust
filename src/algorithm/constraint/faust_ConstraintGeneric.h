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
        ConstraintGeneric() : m_constraintName(CONSTRAINT_NAME_SP),m_nbRows(32),m_nbCols(32) {} // contrainte par defaut (a voir avec Luc)

        ConstraintGeneric(
            const faust_constraint_name& constraintName_,
            const faust_unsigned_int nbRows_,
            const faust_unsigned_int nbCols_) :
                m_constraintName(constraintName_),
                m_nbRows(nbRows_),
                m_nbCols(nbCols_){}

        ConstraintGeneric(const ConstraintGeneric& constraint) :
            m_constraintName(constraint.constraintName),
            m_nbRows(constraint.nbRows),
            m_nbCols(constraint.nbCols){}



        const char* get_type() const;
        const char* get_constraint_name()const;
        const faust_constraint_name get_constraint_type() const;
        bool is_constraint_parameter_int()const;
        bool is_constraint_parameter_real()const;
        bool is_constraint_parameter_mat()const;


        const faust_unsigned_int get_rows() const {return m_nbRows;}
        const faust_unsigned_int get_cols() const {return m_nbCols;}

        virtual void set_default_parameter()=0;
        virtual void check_constraint_name()const=0;
        virtual void project(MatDense<FPP,DEVICE> & mat)const=0;

        ~ConstraintGeneric(){};

        protected:
        /// type of constraint
        const faust_constraint_name m_constraintName;
        // parameter of constraint
        //const parameter_type parameter;
        const faust_unsigned_int m_nbRows;
        const faust_unsigned_int m_nbCols;

        private :
        static const char * m_className;


    };

}

/////// functions useful for parsing config_file (xml,matio,txt_file) //////
bool is_constraint_name_int(const char * type);

bool is_constraint_name_real(const char * type);

bool is_constraint_name_mat(const char * type);

faust_constraint_name get_equivalent_constraint(const char * type);
//! \fn get_type_constraint
//! \brief Information about the type of the constraint
//! \param type
//! \return 0 if type correspond to an integer constraint
//! \return 1 if type correspond to a FPP (Floating Point Precision) constraint
//! \return 2 if type correspond to a Dense Matrix constraint
//! \return else throw an error
int get_type_constraint(const char * type);

#include "faust_ConstraintGeneric.hpp"

#endif
