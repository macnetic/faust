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
				if(opA == 'N' && opB == 'N')
				{
					nopA = 'T';
					nopB = 'T';
					C.resize(B.getNbCol(), A.getNbRow());
					spgemm(B, A, C, alpha, beta, nopB, nopA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'T')
				{
					nopA = 'T';
					C.resize(B.getNbRow(), A.getNbRow());
					spgemm(B, A, C, alpha, beta, opB, nopA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'N')
				{
					nopB = 'T';
					C.resize(B.getNbCol(), A.getNbCol());
					spgemm(B, A, C, alpha, beta, nopB, opA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'T')
				{
					C.resize(B.getNbRow(), A.getNbCol());
					spgemm(B, A, C, alpha, beta, opB, opA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'H')
				{
					nopA = 'H';
					C.resize(B.getNbRow(), A.getNbRow());
					spgemm(B, A, C, alpha, beta, opB, nopA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'N')
				{
					nopB = 'H';
					C.resize(B.getNbCol(), A.getNbCol());
					spgemm(B, A, C, alpha, beta, nopB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'H')
				{
					C.resize(B.getNbRow(), A.getNbCol());
					spgemm(B, A, C, alpha, beta, opB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'T')
				{
					nopA = 'N';
					MatSparse<FPP, GPU2> nB(B);
					nB.conjugate();
					nopB = 'N';
					C.resize(nB.getNbRow(), A.getNbCol());
					spgemm(nB, A, C, alpha, beta, nopB, nopA);
					C.adjoint();
				}
				else if(opA == 'T' && opB == 'H')
				{
					MatDense<FPP, GPU2> nA(A);
					nA.conjugate();
					nopA = 'N';
					nopB = 'N';
					C.resize(B.getNbRow(), nA.getNbCol());
					spgemm(B, nA, C, alpha, beta, nopB, nopA);
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
			const MatSparse<FPP, GPU2>* spA = nullptr;
			const MatSparse<FPP, GPU2>* spB = nullptr;
			const MatDense<FPP, GPU2>* dsA = nullptr;
			const MatDense<FPP, GPU2>* dsB = nullptr;
			const MatBSR<FPP, GPU2>* bsrA = nullptr;
			const MatBSR<FPP, GPU2>* bsrB = nullptr;
			// downcast and call the proper function
			spA = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&A);
			if(! spA)
			{
				dsA = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&A);
				if (! dsA)
					bsrA = dynamic_cast<const Faust::MatBSR<FPP, GPU2>*>(&A);
			}
			spB = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&B);
			if(! spB)
			{
				dsB = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&B);
				if (! dsB)
					bsrB = dynamic_cast<const Faust::MatBSR<FPP, GPU2>*>(&B);
			}
			if(spA && spB)
				spgemm(*spA, MatDense<FPP, GPU2>(*spB), C, alpha, beta, typeA, typeB);
			else if(spA && dsB)
				spgemm(*spA, *dsB, C, alpha, beta, typeA, typeB);
			else if(spB && dsA)
				spgemm(*dsA, *spB, C, alpha, beta, typeA, typeB);
			else if (dsA && dsB)
				gemm(*dsA, *dsB, C, alpha, beta, typeA, typeB);
			else if (bsrA && dsB)
				bsrgemm(*bsrA, *dsB, C, alpha, beta, typeA, typeB);
			else if (bsrA && spB)
				bsrgemm(*bsrA, MatDense<FPP, GPU2>(*spB), C, alpha, beta, typeA, typeB);
			else if (bsrA && bsrB)
				bsrgemm(*bsrA, bsrB->to_dense(), C, alpha, beta, typeA, typeB); // TODO: consider also converting bsrB to MatDense, depending on the weight
			else if (bsrB && dsA)
				bsrgemm(*dsA, *bsrB, C, alpha, beta, typeA, typeB);
			else if (bsrB && spA)
				bsrgemm(MatDense<FPP, GPU2>(*spA), *bsrB, C, alpha, beta, typeA, typeB); // TODO: consider also converting bsrB to MatDense, depending on the weight // to test
			else
				throw std::runtime_error("Unsupported matrix type in faust_linear_algebra_gpu gemm_gen");

		}

	template<typename FPP>
		void bsrgemm(const MatBSR<FPP,GPU2> & A,const MatDense<FPP,GPU2> & B, MatDense<FPP,GPU2> & C,const FPP & alpha, const FPP & beta, char opA, char opB)
		{
			MatBSR<FPP, GPU2>::bsrgemm(A, B, C, alpha, beta, opA, opB);
		}

	template<typename FPP>
		void bsrgemm(const MatDense<FPP,GPU2> & A, const MatBSR<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP & alpha, const FPP & beta, char opA, char opB, int impl_meth/* = 1*/)
		{
			//TODO: benchmark the two methods (impl_meth == 1 and 2)
			if (impl_meth == 1)
			{
				// transpose / adjoint the product to rely on other signature of bsrgemm (MatSparse B as lhs matrix -- i.e. A)
				char nopA, nopB;
				if(opA == 'N' && opB == 'N')
				{
					nopA = 'T';
					nopB = 'T';
					C.resize(B.getNbCol(), A.getNbRow());
					bsrgemm(B, A, C, alpha, beta, nopB, nopA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'T')
				{
					nopA = 'T';
					C.resize(B.getNbRow(), A.getNbRow());
					bsrgemm(B, A, C, alpha, beta, opB, nopA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'N')
				{
					nopB = 'T';
					C.resize(B.getNbCol(), A.getNbCol());
					bsrgemm(B, A, C, alpha, beta, nopB, opA);
					C.transpose();
				}
				else if(opA == 'T' && opB == 'T')
				{
					C.resize(B.getNbRow(), A.getNbCol());
					bsrgemm(B, A, C, alpha, beta, opB, opA);
					C.transpose();
				}
				else if(opA == 'N' && opB == 'H')
				{
					nopA = 'H';
					C.resize(B.getNbRow(), A.getNbRow());
					bsrgemm(B, A, C, alpha, beta, opB, nopA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'N')
				{
					nopB = 'H';
					C.resize(B.getNbCol(), A.getNbCol());
					bsrgemm(B, A, C, alpha, beta, nopB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'H')
				{
					C.resize(B.getNbRow(), A.getNbCol());
					bsrgemm(B, A, C, alpha, beta, opB, opA);
					C.adjoint();
				}
				else if(opA == 'H' && opB == 'T')
				{
					nopA = 'N';
					MatBSR<FPP, GPU2> nB(B);
					nB.conjugate();
					nopB = 'N';
					C.resize(nB.getNbRow(), A.getNbCol());
					bsrgemm(nB, A, C, alpha, beta, nopB, nopA);
					C.adjoint();
				}
				else if(opA == 'T' && opB == 'H')
				{

					MatDense<FPP, GPU2> nA(A);
					nA.conjugate();
					nopA = 'N';
					nopB = 'N';
					C.resize(B.getNbRow(), nA.getNbCol());
					bsrgemm(B, nA, C, alpha, beta, nopB, nopA);
					C.adjoint();
				}
			}
			else {
				spgemm(A, B.to_sparse(), C, alpha, beta, opA, opB);
			}

		}

}

