#include "faust_init_from_matio.h"
#include <iostream>
#include "faust_constant.h"
#include "faust_mat.h"

using namespace std;

void print(const faust_mat&, const char*);



int main()
{
   faust_mat A,B,C;


   init_vec_faust_mat_from_matio_cell_array(A, "test2.mat", "s");
   init_vec_faust_mat_from_matio_cell_array(B, "test2.mat", "c");
   init_vec_faust_mat_from_matio_cell_array(C, "test2.mat", "c2");
   init_vec_faust_mat_from_matio_cell_array(C, "test2.mat", "c3");


   print(A, "A");
   print(B, "B");

   (A*=C);
   print(A, " A *=C");

   C.transpose();
   print(C, "C'");
   


return 0;
}


void print(const faust_mat& A, const char* varName)
{
  cout<<varName<<"="<<endl;
  for (int i=0 ; i<A.getNbRow() ; i++)
  { 
    for (int j=0 ; j<A.getNbCol() ;j++)
     {
        cout << A.getData()[j*A.getNbRow()+i]<<"\t";
     }
     cout<<endl;
  }
  cout<<endl;
}

