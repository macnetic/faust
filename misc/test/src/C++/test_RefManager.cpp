#include <cstdlib>
#include "faust_RefManager.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include "faust_constant.h"

using namespace Faust;

void free_cb(void* ref)
{
	cout << "\tcallback frees MatGeneric" << endl;
	delete static_cast<MatGeneric<double,Cpu>*>(ref);
}

using namespace Faust;

int main()
{
	Faust::MatSparse<double,Cpu>* mat1 = Faust::MatSparse<double,Cpu>::randMat(5,5,.2);
	Faust::MatDense<double,Cpu>* mat2 = Faust::MatDense<double,Cpu>::randMat(5,5);
//	double *mat1 = new double[5];
//	double *mat2 = new double[6];
	RefManager manager;
	manager.set_free_cb(free_cb);
	cout << "acquire mat1" << endl;
	manager.acquire(mat1);
	cout << "acquire mat2" << endl;
	manager.acquire(mat2);
	cout << "acquire mat1" << endl;
	manager.acquire(mat1);
	cout << "release mat1" << endl;
	manager.release(mat1);
	cout << "release mat2" << endl;
	manager.release(mat2);
	return EXIT_SUCCESS;
}
