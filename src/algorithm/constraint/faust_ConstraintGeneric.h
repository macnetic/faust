/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
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
   CONSTRAINT_NAME_TOEPLITZ,
   CONSTRAINT_NAME_CIRC,
   CONSTRAINT_NAME_HANKEL
};

template<typename FPP,FDevice DEVICE> class MatDense;

//! \namespace Faust
//! \brief Faust namespace contains the principal class of the project.
namespace Faust {

    // modif AL AL
    template<typename FPP,FDevice DEVICE>
    class MatDense;

//! \class ConstraintGeneric
//! \brief Contains the generic constraint parameters for the hierarchical factorization. See following table for more precision about the type of constraint. <br>
//! <img src="../../doc/html/constraint.png" alt="constraint parameters" width=800px />
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
		m_constraintName(constraint.m_constraintName),
		m_nbRows(constraint.m_nbRows),
		m_nbCols(constraint.m_nbCols){}



	    template<typename FPP,FDevice DEVICE, typename FPP2=double>
		const char* get_type() const;
	    const char* get_constraint_name()const;
	    const faust_constraint_name get_constraint_type() const;
	    template<typename FPP,FDevice DEVICE, typename FPP2=double>
		bool is_constraint_parameter_int()const;
	    template<typename FPP,FDevice DEVICE, typename FPP2=double>
		bool is_constraint_parameter_real()const;
	    template<typename FPP,FDevice DEVICE, typename FPP2=double>
		bool is_constraint_parameter_mat()const;


	    const faust_unsigned_int get_rows() const {return m_nbRows;}
	    const faust_unsigned_int get_cols() const {return m_nbCols;}

		virtual void Display() const;
	    virtual void set_default_parameter()=0;
	    virtual void check_constraint_name()const=0;
		template<typename FPP, FDevice DEVICE>
		/*virtual*/ void project(MatDense<FPP,DEVICE> & mat)const;//=0; //template with (pure) virtual not authorized (otherwise it must be templates from class, not function)
		template<typename FPP, FDevice DEVICE, typename FPP2>
			/*virtual*/ void project(Faust::MatDense<FPP, DEVICE>&) const /*=0*/;
	    virtual ~ConstraintGeneric(){};

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
