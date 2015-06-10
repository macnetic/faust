#include "faust_init_from_matio.h"
#include <iostream>
#include "faust_constant.h"
#include "faust_mat.h"

using namespace std;

void print(const faust_mat&, const char*);



int main()
{
   faust_mat A,B,C;


   init_faust_mat_from_matio_mat(A, "test.mat", "A");
   init_faust_mat_from_matio_mat(B, "test.mat", "B");
   init_faust_mat_from_matio_mat(C, "test.mat", "C");

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

