#include "faust_init_from_matio_params.h"

template void add_constraint<double>(std::vector<const faust_constraint_generic*> & consS,matvar_t* cons_var);
template void add_constraint<float>(std::vector<const faust_constraint_generic*> & consS,matvar_t* cons_var);