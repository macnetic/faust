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
		void gemm_gen(const MatGeneric<FPP,GPU2> & A, const MatGeneric<FPP,GPU2> & B, MatDense<FPP,GPU2> & C, const FPP  alpha, const FPP  beta, char  typeA, char  typeB)
		{
			const MatSparse<FPP, GPU2>* spA;
			const MatSparse<FPP, GPU2>* spB;
			const MatDense<FPP, GPU2>* dsA;
			const MatDense<FPP, GPU2>* dsB;
			// downcast an call the proper function
			spA = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&A);
			if(! spA)
				dsA = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&A);
			spB = dynamic_cast<const Faust::MatSparse<FPP,GPU2>*>(&B);
			if(! spB)
				dsB = dynamic_cast<const Faust::MatDense<FPP,GPU2>*>(&B);
			if(spA && spB)
				throw std::runtime_error("gemm on two MatSparse is not supported.");
			else if(spA)
				throw std::runtime_error("spgemm is not supported yet on GPU2."); //TODO: after spgemm impl.
//				spgemm(*spA, *dsB, C, alpha, beta, typeA, typeB);
			else if(spB)
				throw std::runtime_error("spgemm is not supported yet on GPU2."); //TODO: after spgemm impl.
//				spgemm(*dsA, *spB, C, alpha, beta, typeA, typeB);
			else
				gemm(*dsA, *dsB, C, alpha, beta, typeA, typeB);
		}
}
