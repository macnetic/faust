namespace Faust
{
	template<typename FPP>
	MatGeneric<FPP,GPU2>::~MatGeneric()
	{
	}

	template<typename FPP>
	MatGeneric<FPP,GPU2>::MatGeneric() : is_zeros(false), is_identity(false)
	{
	}

//	//TODO: (factorize with faust_MatGeneric) only one template function with FDevice DEV as second param
//	template<typename FPP>
//		MatGeneric<FPP,GPU2>* optimize(Faust::MatDense<FPP,GPU2> const & M,Faust::MatSparse<FPP,GPU2> const & S)
//		{
//			//std::cout<<"DEBUT OPTIMIZE "<<std::endl;
//
//
//			if ( (M.getNbCol() != S.getNbCol()) | (M.getNbRow() != S.getNbRow()) )
//				handleError("Faust::MatGeneric::", " Faust::optimize : matrix M and S have not the same size");
//
//			Faust::Vect<FPP,GPU2> x_dense(M.getNbCol());
//
//			for (int i=0;i<M.getNbCol();i++)
//			{
//				x_dense.set_coeff(i, i*.005);
//			}
//
//			Faust::Vect<FPP,GPU2> const x(x_dense);
//			Faust::Vect<FPP,GPU2> x_sparse(x_dense);
//
//			int nb_mult=10;
//			Faust::Timer t_dense,t_sparse;
//			for (int i=0;i<nb_mult;i++)
//			{
//				x_sparse=x;
//				x_dense=x;
//				t_sparse.start();
//				S.multiply(x_sparse,'N');
//				t_sparse.stop();
//				t_dense.start();
//				M.multiply(x_dense,'N');
//				t_dense.stop();
//
//			}
//			//float density = ((float)S.getNonZeros())/((float)(S.getNbRow()*S.getNbCol()));
//			float density = S.density();
//			//std::cout<<" density "<<density<<std::endl;
//			//std::cout<<" tps sparse "<<t_sparse.get_time()<<std::endl;
//			//std::cout<<" tps dense "<<t_dense.get_time()<<std::endl;
//
//			//if (M.getNbCol() != M.getNbRow())
//			if (t_sparse.get_time() <= t_dense.get_time())
//			{
//				//std::cout<<" CHOICE SPARSE "<<t_dense.get_time()<<std::endl;
//				return new MatSparse<FPP,GPU2>(S);
//			}else
//			{
//				//std::cout<<" CHOICE DENSE "<<t_dense.get_time()<<std::endl;
//				return new MatDense<FPP,GPU2>(M);
//			}
//
//			//std::cout<<"FIN OPTIMIZE "<<t_sparse.get_time()<<std::endl;
//
//
//		}
}
