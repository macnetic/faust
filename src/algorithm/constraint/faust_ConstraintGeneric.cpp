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
#include "faust_ConstraintGeneric.h"
#include "faust_ConstraintInt.h"
#include "faust_ConstraintFPP.h"
#include <typeinfo>
#include "faust_exception.h"

#include "faust_ConstraintMat.h"

// modif AL AL
#include "faust_ConstraintType.h"




bool is_constraint_name_int(const char * type)
{
	bool is_const_int =((strcmp(type,"sp") == 0) || (strcmp(type,"sppos")==0));
	is_const_int = ((is_const_int) || ((strcmp(type,"spcol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splincol") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"splin") == 0)));
//	is_const_int = ((is_const_int) || ((strcmp(type,"blkdiag") == 0)));
	is_const_int = ((is_const_int) || ((strcmp(type,"skperm") == 0)));
	return is_const_int;
}


bool is_constraint_name_real(const char * type)
{
	return ((strcmp(type,"normcol") == 0) || (strcmp(type,"normlin")==0));
}

bool is_constraint_name_mat(const char * type)
{
	return (strcmp(type,"supp") == 0) || (strcmp(type,"const")==0) || ! strcmp(type, "toeplitz") || ! strcmp(type, "circ") || ! strcmp(type, "blockdiag") || !strcmp(type, "blkdiag") || ! strcmp (type, "hankel") || ! strcmp(type, "id") || ! strcmp(type, "proj_id");
}

faust_constraint_name get_equivalent_constraint(const char * type)
{

	if (!strcmp(type,"sp"))
		return CONSTRAINT_NAME_SP;
	if (!strcmp(type,"spcol"))
		return CONSTRAINT_NAME_SPCOL;
	if  (!strcmp(type,"splin"))
		return CONSTRAINT_NAME_SPLIN;
	if  (!strcmp(type,"normcol"))
		return CONSTRAINT_NAME_NORMCOL;
	if  (!strcmp(type,"splincol"))
		return CONSTRAINT_NAME_SPLINCOL;
	if  (!strcmp(type,"const"))
		return CONSTRAINT_NAME_CONST;
	if  (!strcmp(type,"sppos"))
		return CONSTRAINT_NAME_SP_POS;
	if  (!strcmp(type,"blkdiag") || !strcmp(type, "blockdiag"))
		return CONSTRAINT_NAME_BLKDIAG;
	if  (!strcmp(type,"supp"))
		return CONSTRAINT_NAME_SUPP;
	if  (!strcmp(type,"normlin"))
		return CONSTRAINT_NAME_NORMLIN;
	if(!strcmp(type, "skperm"))
		return CONSTRAINT_NAME_SKPERM;
	if(!strcmp(type, "id"))
		return CONSTRAINT_NAME_ID;
	handleError("Faust::ConstraintGeneric","get_equivalent_constraint : Unknown type of constraint");
}


int get_type_constraint(const char * type)
{
	bool is_const_int = is_constraint_name_int(type);
	bool is_const_real = is_constraint_name_real(type);
	bool is_const_mat = is_constraint_name_mat(type);

	if (is_const_int)
		return 0;
	if (is_const_real)
		return 1;
	if (is_const_mat)
		return 2;
	handleError("Faust::ConstraintGeneric","::add_constraint : invalid constraint type");
}

const faust_constraint_name Faust::ConstraintGeneric::get_constraint_type() const
{
   return m_constraintName;
}





const char* Faust::ConstraintGeneric::get_constraint_name()const
{
   switch(m_constraintName)
   {
      case CONSTRAINT_NAME_SP:
         return "CONSTRAINT_NAME_SP";
      case CONSTRAINT_NAME_SPCOL:
         return "CONSTRAINT_NAME_SPCOL";
      case CONSTRAINT_NAME_SPLIN:
         return "CONSTRAINT_NAME_SPLIN";
      case CONSTRAINT_NAME_NORMCOL:
         return "CONSTRAINT_NAME_NORMCOL";
      case CONSTRAINT_NAME_SPLINCOL:
         return "CONSTRAINT_NAME_SPLINCOL";
	  case CONSTRAINT_NAME_SKPERM:
		 return "CONSTRAINST_NAME_SKPERM";
      case CONSTRAINT_NAME_CONST:
         return "CONSTRAINT_NAME_CONST";
      case CONSTRAINT_NAME_SP_POS:
         return "CONSTRAINT_NAME_SP_POS";
      case CONSTRAINT_NAME_BLKDIAG:
         return "CONSTRAINT_NAME_BLKDIAG";
      case CONSTRAINT_NAME_SUPP:
         return "CONSTRAINT_NAME_SUPP";
      case CONSTRAINT_NAME_NORMLIN:
         return "CONSTRAINT_NAME_NORMLIN";
	  case CONSTRAINT_NAME_CIRC:
		 return "CONSTRAINT_NAME_CIRC";
	  case CONSTRAINT_NAME_HANKEL:
		 return "CONSTRAINT_NAME_HANKEL";
	  case CONSTRAINT_NAME_TOEPLITZ:
		 return "CONSTRAINT_NAME_TOEPLITZ";
	  case CONSTRAINT_NAME_ID:
		 return "CONSTRAINT_NAME_ID";
      default:
         return "unknown constraint name";
   }
}

void Faust::ConstraintGeneric::Display() const
{
	std::cout << this->get_constraint_name();
	std::cout<<" nb_row: "<< this->get_rows();
	std::cout<<" nb_col: "<< this->get_cols();
	std::cout<<" normalized :"<< get_normalizing() << std::endl;
	std::cout<<" pos :"<< get_pos() << std::endl;
}

const char * Faust::ConstraintGeneric::m_className="Faust::ConstraintGeneric::";
