#ifndef __FAUST_CONSTRAINT_H__
#define __FAUST_CONSTRAINT_H__

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

template<typename parameter_type>
class faust_constraint
{
   public:
      faust_constraint() : constraint_name(); // contrainte par defaut (a voir avec Luc)

      faust_constraint(
         const faust_constraint_name& constraint_name_, 
         const parameter_type& parameter_, 
         const int nb_rows_, 
         const int nb_cols_) :
            constraint_name(constraint_name_),
            parameter(parameter_),
            nb_rows(nb_rows_),
            nb_cols(nb_cols_){}

      faust_constraint(const faust_constraint& constraint) :
         constraint_name(constraint.constraint_name),
         parameter(constraint.parameter),
         nb_rows(constraint.nb_rows),
         nb_cols(constraint.nb_cols){}

      const faust_constraint_name getConstraintType() const {return constraint_name;}
      const parameter_type getParameter() const {return parameter;};
      const int getRows() const {return nb_rows};
      const int getCols() const {return nb_cols};
 
      ~faust_constraint(){};

   private:
      // type of constraint
      const faust_constraint_name constraint_name;
      // parameter of constraint
      const parameter_type parameter;

      const int nb_rows;
      const int nb_cols;


      
}

#endif
