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


	template<typename FPP>
		void MatGeneric<FPP,GPU2>::Display() const
		{//TODO: factorize with CPU code
			std::cout << to_string();
		}

	template<typename FPP>
		std::string MatGeneric<FPP, GPU2>::to_string(MatType type, const bool transpose /* set to false by default */, const bool displaying_small_mat_elts) const
		{//TODO: factorize with CPU code
			return this->to_string(getNbRow(), getNbCol(), transpose, this->density(), this->getNonZeros(), this->is_identity, type);
		}

	template<typename FPP>
		std::string MatGeneric<FPP, GPU2>::to_string(const bool transpose /* set to false by default */, const bool displaying_small_mat_elts) const
		{//TODO: factorize with CPU code
			MatType type = Dense;
			if(dynamic_cast<const MatSparse<FPP, GPU2>*>(this))
				type = Sparse;
			else if (dynamic_cast<const MatDense<FPP, GPU2>*>(this))
				type = Dense;
			else if (dynamic_cast<const MatDiag<FPP>*>(this))
				type = Diag;
			else if (dynamic_cast<const MatBSR<FPP, GPU2>*>(this))
				type = BSR;
			else if (dynamic_cast<const MatPerm<FPP, GPU2>*>(this))
				type = Perm;
			else if (dynamic_cast<const MatButterfly<FPP, GPU2>*>(this))
				type = Butterfly;
			else
				throw std::runtime_error("Unhandled matrix type in MatGeneric::to_string()"); // shouldn't happen
			return this->to_string(getNbRow(), getNbCol(), transpose, this->density(), this->getNonZeros(), this->is_identity, type);
		}

	template<typename FPP>
		std::string MatGeneric<FPP, GPU2>::get_scalar_type_str()
		{ //TODO: factorize with CPU code
			std::string type_str;
			if (! std::is_same<FPP,Real<FPP>>::value)
				type_str = "complex";
			else
			{
				if(std::is_same<FPP, float>::value)
					type_str = "float";
				else
					type_str = "double";
			}
			return type_str;
		}


template<typename FPP>
std::string Faust::MatGeneric<FPP,GPU2>::to_string(int32_t nrows, int32_t ncols, bool transpose, Real<FPP> density, int32_t nnz, bool is_identity, MatType type) const
{
	std::ostringstream str;
	str << " (" << MatGeneric<FPP,GPU2>::get_scalar_type_str() << ") ";
	str << "[GPU] ";
	if(type == Dense)
		str << "DENSE,";
	else if(type == Sparse)
		str << "SPARSE,";
	else if(type == Diag)
		str << "DIAG,";
	else if(type == BSR)
		str << "BSR,";
	else if(type == Perm)
		str << "PERM,";
	else if(type == Butterfly)
		str << "Butterfly,";
	else if(type == None)
		str << "UNKNOWN MATRIX TYPE,";
	else
		throw std::runtime_error("Invalid MatType passed to MatGeneric::to_string()");
	str << " size ";
	if(transpose)
		str << ncols << "x" << nrows;
	else
		str << nrows << "x" << ncols;
	if(type == BSR) //TODO: refactor with above almost same code into a new MatBSR member function
		str << dynamic_cast<const MatBSR<FPP, GPU2>*>(this)->to_string_blocks(transpose);
	str << ", density "<< density <<", nnz "<< nnz; 
	if(type == BSR)
		str <<  " (nnz blocks: " << dynamic_cast<const MatBSR<FPP, GPU2>*>(this)->getNBlocks() << ")";
	str << ", addr: " << this;
	//TODO: add gpu_mod data addr
	str <<std::endl;
	if (is_identity)
		str <<" identity matrix flag" << std::endl;
	return str.str();
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
