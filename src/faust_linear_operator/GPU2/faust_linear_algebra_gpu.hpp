namespace Faust
{

	template<typename FPP>
		void gemm(const MatDense<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
		{
			MatDense<FPP, GPU2>::gemm(A, B, C, alpha, beta, typeA, typeB);
		}

	template<typename FPP>
		void gemv(const MatDense<FPP, GPU2> &A, const Vect<FPP, GPU2> &b, Vect<FPP, GPU2> &c, const FPP& alpha, const FPP& beta, const char opA)
		{
			MatDense<FPP, GPU2>::gemv(A, b, c, alpha, beta);
		}

	template<typename FPP>
		void spgemm(const MatSparse<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB)
		{
			MatSparse<FPP, GPU2>::spgemm(A, B, C, alpha, beta, opA, opB);
		}

	template<typename FPP>
		void spgemm(const MatDense<FPP,GPU2> & A, const MatSparse<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP & alpha, const FPP & beta, char opA, char opB, int impl_meth/* = 1*/)
		{
			//TODO: benchmark the two methods (impl_meth == 1 and 2)
			if (impl_meth == 1)
			{
				// transpose / adjoint the product to rely on other signature of spgemm (MatSparse B as lhs matrix -- i.e. A)
				char nopA, nopB;
				MatDense<FPP, GPU2> nA(A);
				MatSparse<FPP, GPU2> nB(B);
				if(opA == 'N' && opB == 'N')
				{
					nopA = 'T';
					nopB = 'T';
					C.resize(nB.getNbCol(), nA.getNbRow());
					spgemm(nB, nA, C, alpha, beta, nopB, nopA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'T')
				{
					nopA = 'T';
					C.resize(nB.getNbRow(), nA.getNbRow());
					spgemm(nB, nA, C, alpha, beta, opB, nopA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'N')
				{
					nopB = 'T';
					C.resize(nB.getNbCol(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, nopB, opA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'T')
				{
					C.resize(nB.getNbRow(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, opB, opA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'H')
				{
					nopA = 'H';
					C.resize(nB.getNbRow(), nA.getNbRow());
					spgemm(nB, nA, C, alpha, beta, opB, nopA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'N')
				{
					nopB = 'H';
					C.resize(nB.getNbCol(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, nopB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'H')
				{
					C.resize(nB.getNbRow(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, opB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'T')
				{
					nopA = 'N';
					nB.conjugate();
					nopB = 'N';
					C.resize(nB.getNbRow(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, nopB, nopA);
					C.adjoint();
				}
				else if(opA == 'T' && opB == 'H')
				{
					nA.conjugate();
					nopA = 'N';
					nopB = 'N';
					C.resize(nB.getNbRow(), nA.getNbCol());
					spgemm(nB, nA, C, alpha, beta, nopB, nopA);
					C.adjoint();
				}
			}
			else {
				spgemm(MatSparse<FPP, GPU2>(A), MatDense<FPP, GPU2>(B), C, alpha, beta, opA, opB);
			}
		}


	template<typename FPP>
		void gemm_gen(const MatGeneric<FPP,GPU2> & A, const MatGeneric<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
		{
			const MatSparse<FPP, GPU2>* spA;
			const MatSparse<FPP, GPU2>* spB;
			const MatDense<FPP, GPU2>* dsA;
			const MatDense<FPP, GPU2>* dsB;
			// downcast and call the proper function
			spA = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&A);
			if(! spA)
				dsA = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&A);
			spB = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&B);
			if(! spB)
				dsB = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&B);
			if(spA && spB)
				throw std::runtime_error("gemm on two MatSparse is not supported.");
			else if(spA)
				spgemm(*spA, *dsB, C, alpha, beta, typeA, typeB);
			else if(spB)
				spgemm(*dsA, *spB, C, alpha, beta, typeA, typeB);
			else
				gemm(*dsA, *dsB, C, alpha, beta, typeA, typeB);
		}
}
