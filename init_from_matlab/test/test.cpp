#include "faust_init_from_matio.h"
#include <iostream>
#include "faust_constant.h"

using namespace std;


int main()
{
   faust_mat A;
   faust_real B,E,G;
   int C;
   bool D,F;

   init_faust_mat_from_matio_mat(A, "test.mat", "A");
   B = init_faust_mat_from_matio_double("test.mat", "B");
   C = init_faust_mat_from_matio_int("test.mat", "C");
   D = init_faust_mat_from_matio_bool("test.mat", "D");
   F = init_faust_mat_from_matio_bool("test.mat", "F");
   G = init_faust_mat_from_matio_double("test.mat", "G");



   
   


return 0;
}


