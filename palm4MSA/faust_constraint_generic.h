#ifndef __FAUST_CONSTRAINT_GENERIC_H__
#define __FAUST_CONSTRAINT_GENERIC_H__
#include <string>
#include <iostream>


enum faust_constraint_name
{
   CONSTRAINT_NAME_SP,
   CONSTRAINT_NAME_SPCOL,
   CONSTRAINT_NAME_SPLIN,
   CONSTRAINT_NAME_NORMCOL,
   CONSTRAINT_NAME_SPLINCOL,
   CONSTRAINT_NAME_L0PEN,
   CONSTRAINT_NAME_L1PEN,
   CONSTRAINT_NAME_CONST,
   CONSTRAINT_NAME_WAV,
   CONSTRAINT_NAME_SP_POS,
   CONSTRAINT_NAME_BLKDIAG,
   CONSTRAINT_NAME_SPLIN_TEST,
   CONSTRAINT_NAME_SUPP,
   CONSTRAINT_NAME_NORMLIN,
   CONSTRAINT_NAME_TOEPLITZ
};


std::string getConstraintType(faust_constraint_name);	


//template<typename parameter_type>
class faust_constraint_generic
{
   public:
      faust_constraint_generic() : constraint_name(CONSTRAINT_NAME_SP),nb_rows(32),nb_cols(32) {} // contrainte par defaut (a voir avec Luc)

      faust_constraint_generic(
         const faust_constraint_name& constraint_name_, 
         //const parameter_type& parameter_, 
         const int nb_rows_, 
         const int nb_cols_) :
            constraint_name(constraint_name_),
            //parameter(parameter_),
            nb_rows(nb_rows_),
            nb_cols(nb_cols_){}

      faust_constraint_generic(const faust_constraint_generic& constraint) :
         constraint_name(constraint.constraint_name),
         //parameter(constraint.parameter),
         nb_rows(constraint.nb_rows),
         nb_cols(constraint.nb_cols){}

      const char* get_constraint_name()const;
      const faust_constraint_name getConstraintType() const {return constraint_name;}
      //const parameter_type getParameter() const {return parameter;};
      const int getRows() const {return nb_rows;}
      const int getCols() const {return nb_cols;}
	  void Display() const {std::cout<<get_constraint_name()<<" DIM : "<<nb_rows<<" "<<nb_cols<<std::endl;}	
	
      virtual void set_default_parameter()=0;
      virtual void check_constraint_name()const=0;
 
      ~faust_constraint_generic(){};


   protected:
      // type of constraint
      const faust_constraint_name constraint_name;
      // parameter of constraint
      //const parameter_type parameter;

      const int nb_rows;
      const int nb_cols;


      
};

#endif
