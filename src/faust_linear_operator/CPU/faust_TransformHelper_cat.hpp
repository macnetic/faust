namespace Faust
{
	template<typename FPP>
	TransformHelper<FPP,Cpu>* vertcat(const std::vector<TransformHelper<FPP,Cpu>*> & THs)
	{
		TransformHelper<FPP,Cpu>* output = nullptr;
		// find the maximum size TransformHelper id in THs
		size_t max_size = 0;
		for(auto th: THs)
			if(th->size() > max_size)
				max_size = th->size();
//		std::cout << "max num of factors is: " << max_size << endl;
		std::vector<MatGeneric<FPP,Cpu>*> facts(max_size+1);
		// list to store if the i-th factor of THs[j] has been copied (converted from MatDense to Matsparse)
		std::vector<bool> inter_fac_copied(THs.size());
		inter_fac_copied.assign(THs.size(), false);
		const MatGeneric<FPP,Cpu> *tmp_fac_gen;
		// list of diag matrices to form one factor of the concatenated faust
		std::vector<const MatSparse<FPP,Cpu>*> diag_mats(THs.size());
		std::vector<const MatSparse<FPP,Cpu>*> eye_mats(THs.size());
		eye_mats.assign(THs.size(), nullptr);
		faust_unsigned_int out_fac_nnz, out_fac_nrows, out_fac_ncols;
		int * out_fac_colind, *out_fac_rowptr;
		FPP* out_fac_values;
		MatSparse<FPP, Cpu>* out_fac;
		size_t out_fact_nnz_offset, out_fact_row_offset, out_fact_col_offset;
		for(int i=0;i < max_size; i++)
		{
			// construct the i-th blockdiag factor of the output concatenation
			out_fac_nnz = 0;
			out_fac_nrows = 0;
			out_fac_ncols = 0;
			for(int j=0; j < THs.size(); j++)
			{
				auto faust = THs[j];
//				std::cout << "faust:" << faust << std::endl;
				if(i < faust->size())
				{
					tmp_fac_gen = faust->get_gen_fact(i);
					if((inter_fac_copied[j] = (tmp_fac_gen->getType() == MatType::Dense)))
					{
						diag_mats[j] = new MatSparse<FPP,Cpu>(*dynamic_cast<const MatDense<FPP,Cpu>*>(tmp_fac_gen));
						if(faust->is_transposed) const_cast<MatSparse<FPP,Cpu>*>(diag_mats[j])->transpose();
						if (faust->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(diag_mats[j])->conjugate();
					}
					else
					{
						auto dm = dynamic_cast<const MatSparse<FPP,Cpu>*>(tmp_fac_gen);
						if(faust->is_transposed || faust->is_conjugate)
						{
							dm = new MatSparse<FPP,Cpu>(dm->getNonZeros(), dm->getNbRow(), dm->getNbCol(), dm->getValuePtr(), dm->getRowPtr(), dm->getColInd(), faust->is_transposed);
							if (faust->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(dm)->conjugate();
							inter_fac_copied[j] = true;
						}
						diag_mats[j] = dm;
					}
					out_fac_nnz += diag_mats[j]->getNonZeros();
					out_fac_ncols += faust->get_fact_nb_cols(i);
					out_fac_nrows += faust->get_fact_nb_rows(i);
				}
				else
				{ // faust number of factors is less than max_size
					//fill remaining factors with identity
					tmp_fac_gen = faust->get_gen_fact(faust->size()-1);
					// number of ones
					auto faust_ncols = faust->get_fact_nb_cols(faust->size()-1);
					out_fac_nnz += faust_ncols;//tmp_fac_gen->getNbCol();
					out_fac_ncols += faust_ncols; //F_fac->getNbCol()+tmp_fac_gen->getNbCol();
					out_fac_nrows += faust_ncols;//tmp_fac_gen->getNbCol();
					// opt. create the id factor only once
					if(i == faust->size())
					{
						eye_mats[j] = MatSparse<FPP,Cpu>::eye(faust_ncols, faust_ncols);//(tmp_fac_gen->getNbCol(), tmp_fac_gen->getNbCol());
						//TODO: delete all eye_mats after the main loop
					}
					diag_mats[j] = eye_mats[j];
				}
			}
			out_fac_values = new FPP[out_fac_nnz];
			out_fac_colind = new int[out_fac_nnz];
			out_fac_rowptr = new int[out_fac_nrows+1];
			out_fac_rowptr[0] = 0;
			out_fact_nnz_offset = 0;
			out_fact_row_offset = 1;
			out_fact_col_offset = 0;
			// copy the data from all diag_mats
			for(int j=0; j < diag_mats.size();j++)
			{
//				std::cout << "j: " << j << " " << diag_mats[j]->norm() << std::endl;
				memcpy(out_fac_values+out_fact_nnz_offset, diag_mats[j]->getValuePtr(), sizeof(FPP)*diag_mats[j]->getNonZeros());
				memcpy(out_fac_colind+out_fact_nnz_offset, diag_mats[j]->getColInd(), sizeof(int)*diag_mats[j]->getNonZeros());
				// ignore first element of getRowPtr()
				memcpy(out_fac_rowptr+out_fact_row_offset, diag_mats[j]->getRowPtr()+1, sizeof(int)*diag_mats[j]->getNbRow());
				if(j > 0)
				{
					out_fact_col_offset += diag_mats[j-1]->getNbCol();
//					std::cout << "j > 0" << std::endl;
					// the indices must be shifted by the previous mat num of cols
					for(auto k = out_fact_nnz_offset; k < diag_mats[j]->getNonZeros()+out_fact_nnz_offset; k++)
						out_fac_colind[k] += out_fact_col_offset; //out_fac_colind[out_fact_nnz_offset-1]+1;
					// shift rowptr according to previous diag_mat
					for(auto k = out_fact_row_offset; k < diag_mats[j]->getNbRow()+out_fact_row_offset;k++)
						out_fac_rowptr[k] += out_fac_rowptr[out_fact_row_offset-1]; //;diag_mats[j-1]->getRowPtr()[diag_mats[j-1]->getNbRow()];
				}
				out_fact_nnz_offset += diag_mats[j]->getNonZeros();
				out_fact_row_offset += diag_mats[j]->getNbRow();
//				out_fact_col_offset += THs[j]->get_gen_fact(i)->getNbCol();
			}
//			std::cout << "creating out fac i:" << i << std::endl;
//			std::cout << "out_fac_nrows: " << out_fac_nrows << " out_fac_ncols: " << out_fac_ncols << std::endl;
			out_fac = new MatSparse<FPP,Cpu>(out_fac_nnz, out_fac_nrows, out_fac_ncols, out_fac_values, out_fac_rowptr, out_fac_colind);
			facts[i] = out_fac;
			for(int j=0; j < THs.size(); j++)
			{
				if(inter_fac_copied[j])
				{
					delete diag_mats[j];
					inter_fac_copied[j] = false;
				}
			}
			delete[] out_fac_values;
			delete[] out_fac_rowptr;
			delete[] out_fac_colind;
		}
		// add last factor: stack of id matrices
//		std::cout << "creating last out fac:" << std::endl;
//		out_fac_ncols = THs[0]->get_gen_fact(THs[0]->size()-1)->getNbCol();
		out_fac_ncols = THs[0]->getNbCol(); // all fausts have normally same num of cols
		out_fac_nrows = out_fac_ncols*THs.size();
		out_fac_nnz = THs.size()*out_fac_ncols;
		out_fac_values = new FPP[out_fac_nnz];
		out_fac_colind = new int[out_fac_nnz];
		out_fac_rowptr = new int[out_fac_nnz+1];
		for(int i = 0; i < out_fac_nnz;i++)
			out_fac_values[i] = 1;
		for(int i = 0; i < THs.size();i++)
		{
			for(int j=0;j < out_fac_ncols;j++)
			{
				out_fac_colind[i*out_fac_ncols+j] = j;
			}
		}
		out_fac_rowptr[0] = 0;
		for(int i = 1; i < out_fac_nnz+1;i++)
			out_fac_rowptr[i] = out_fac_rowptr[i-1] + 1;
//		std::cout << "out_fac_nrows, out_fac_ncols:" << out_fac_nrows << " " << out_fac_ncols << std::endl;
		out_fac = new MatSparse<FPP,Cpu>(out_fac_nnz, out_fac_nrows, out_fac_ncols, out_fac_values, out_fac_rowptr, out_fac_colind);
		delete[] out_fac_values;
		delete[] out_fac_colind;
		delete[] out_fac_rowptr;
		facts[facts.size()-1] = out_fac;
		output = new TransformHelper<FPP,Cpu>(facts, 1.0, false, false);
		for(int j=0;j<THs.size();j++)
			if(eye_mats[j] != nullptr)
				delete eye_mats[j];
		return output;
	}

		template<typename FPP>
	TransformHelper<FPP,Cpu>* horzcat(const std::vector<TransformHelper<FPP,Cpu>*> & THs)
	{
		TransformHelper<FPP,Cpu> *C, *Ct;
		std::vector<TransformHelper<FPP,Cpu>*> tTHs;
		for(auto t: THs)
			tTHs.push_back(t->transpose());
		C = vertcat(tTHs);
		Ct = C->transpose();
		delete C;
		for(auto t: tTHs)
			delete t;
		return Ct;
	}

		template<typename FPP>
			TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::vertcat(const TransformHelper<FPP,Cpu>* G)
			{
				// NOTE: this function is in this class and not in Faust::Transform
				// because it's simpler to handle the transpose and conjugate cases here
				// the functions get_gen_fact() and get_fact_nb_cols()/rows are
				// really helpful for that.
				//TODO: too big function, refactor
				//TODO: clean the code
				TransformHelper<FPP,Cpu>* T = nullptr;
				// take the greater number of factors from this and G as the concatened Faust number
				// of factors
				std::vector<MatGeneric<FPP,Cpu>*> facts(max(G->size(),this->size())+1);
				const MatSparse<FPP,Cpu> * F_fac, *G_fac;
				const MatGeneric<FPP,Cpu> *tmp_fac_gen;
				MatSparse<FPP,Cpu>* T_fac;
				faust_unsigned_int T_fac_nnz;
				faust_unsigned_int T_fac_nb_rows, T_fac_nb_cols;
				bool F_inter_fac_allocated, G_inter_fac_allocated;
				bool F_last_fac=false, G_last_fac=false;
				int * colind, *rowptr;
				FPP* values;

				if(this->getNbCol() != G->getNbCol()) handleError("TransformHelper::vertcat()","The dimensions must agree.");
				for(faust_unsigned_int i=0; i < facts.size(); i++){
					//				cout << "factor i=" << i <<  " facts.size()=" << facts.size() << endl;
					// F (*this) factor
					if(i < this->size())
					{
						tmp_fac_gen = get_gen_fact(i); //this->transform->data[i];
						if((F_inter_fac_allocated = (tmp_fac_gen->getType() == MatType::Dense)))
						{
							F_fac = new MatSparse<FPP,Cpu>(*dynamic_cast<const MatDense<FPP,Cpu>*>(tmp_fac_gen));
							if(this->is_transposed) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->transpose();
							if(this->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->conjugate();
						}
						else
						{
							F_fac = dynamic_cast<const MatSparse<FPP,Cpu>*>(tmp_fac_gen);
							if(this->is_transposed || this->is_conjugate)
							{
								F_fac = new MatSparse<FPP,Cpu>(F_fac->nnz, F_fac->getNbRow(), F_fac->getNbCol(), F_fac->getValuePtr(), F_fac->getRowPtr(), F_fac->getColInd(), this->is_transposed);
								if(this->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(F_fac)->conjugate();

								F_inter_fac_allocated = true;
							}
						}
						T_fac_nnz = F_fac->getNonZeros(); //+G_fac->getNonZeros()
						T_fac_nb_cols = this->get_fact_nb_cols(i);//F_fac->getNbCol(); //+G_fac part
						T_fac_nb_rows = this->get_fact_nb_rows(i);//F_fac->getNbRow(); //+G_fac part
					}
					else
					{
						//fill remaining factors by identity
						tmp_fac_gen = get_gen_fact(this->size()-1);//this->transform->data[this->size()-1];
						// number of ones
						T_fac_nnz = this->get_fact_nb_cols(this->size()-1);//tmp_fac_gen->getNbCol(); //+G_fac->getNonZeros()
						T_fac_nb_rows = T_fac_nb_cols = T_fac_nnz; //+G_fac part
						// opt. create the id factor only once
						if(!F_last_fac)
						{
							F_last_fac = true;
							F_fac = MatSparse<FPP,Cpu>::eye(T_fac_nnz, T_fac_nnz);
						}
					}
					// G factor
					if(i < G->size())
					{
						tmp_fac_gen = G->get_gen_fact(i);//G->transform->data[i];
						if((G_inter_fac_allocated = (tmp_fac_gen->getType() == MatType::Dense)))
						{
							G_fac = new MatSparse<FPP,Cpu>(*dynamic_cast<const MatDense<FPP,Cpu>*>(tmp_fac_gen));
							if(G->is_transposed) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->transpose();
							if (G->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->conjugate();

						}
						else
						{
							G_fac = dynamic_cast<const MatSparse<FPP,Cpu>*>(tmp_fac_gen);
							if(G->is_transposed || G->is_conjugate)
							{
								G_fac = new MatSparse<FPP,Cpu>(G_fac->nnz, G_fac->getNbRow(), G_fac->getNbCol(), G_fac->getValuePtr(), G_fac->getRowPtr(), G_fac->getColInd(), G->is_transposed);
								if (G->is_conjugate) const_cast<MatSparse<FPP,Cpu>*>(G_fac)->conjugate();
								G_inter_fac_allocated = true;
							}

						}
						T_fac_nnz += G_fac->getNonZeros();
						T_fac_nb_cols += G->get_fact_nb_cols(i);//G_fac->getNbCol();
						T_fac_nb_rows += G->get_fact_nb_rows(i);//G_fac->getNbRow();
					}
					else
					{ //G has less or equal factors than F

						//fill remaining factors by identity
						tmp_fac_gen = G->get_gen_fact(G->size()-1);//G->transform->data[G->size()-1];
						// number of ones
						T_fac_nnz += G->get_fact_nb_cols(G->size()-1);//tmp_fac_gen->getNbCol();
						if(i < facts.size()-1 || !F_last_fac) {
							T_fac_nb_cols = F_fac->getNbCol()+G->get_fact_nb_cols(G->size()-1);//F_fac->getNbCol()+tmp_fac_gen->getNbCol();
						}
						//G_fac will be square id. => nb cols == nb rows
						T_fac_nb_rows = F_fac->getNbRow()+G->get_fact_nb_cols(G->size()-1);//tmp_fac_gen->getNbCol();
						// opt. create the id factor only once
						if(!G_last_fac)
						{
							G_last_fac = true;
							G_fac = MatSparse<FPP,Cpu>::eye(G->get_fact_nb_cols(G->size()-1),G->get_fact_nb_cols(G->size()-1));//(tmp_fac_gen->getNbCol(), tmp_fac_gen->getNbCol());
						}
					}
					values = new FPP[T_fac_nnz];
					colind = new int[T_fac_nnz];
					rowptr = new int[T_fac_nb_rows+1];
					// group the data from F and G
					memcpy(values, F_fac->getValuePtr(), sizeof(FPP)*F_fac->getNonZeros());
					memcpy(values+F_fac->getNonZeros(), G_fac->getValuePtr(),
							sizeof(FPP)*G_fac->getNonZeros());
					assert(T_fac_nnz == F_fac->getNonZeros()+G_fac->getNonZeros());
					assert(T_fac_nb_rows == F_fac->getNbRow()+G_fac->getNbRow());
					/****** col indices ******/
					// F indices don't change
					memcpy(colind, F_fac->getColInd(), sizeof(int)*F_fac->getNonZeros());
					// G indices are shifted by F->getNbCol()
					memcpy(colind+F_fac->getNonZeros(), G_fac->getColInd(), sizeof(int)*G_fac->getNonZeros());
					//shift G indices
					if(i < facts.size()-1)
						for(faust_unsigned_int j=F_fac->getNonZeros();j<F_fac->getNonZeros()+G_fac->getNonZeros();j++)
							colind[j] += F_fac->getNbCol();
					/***** row indices *****/
					memcpy(rowptr, F_fac->getRowPtr(), sizeof(int)*(F_fac->getNbRow()+1));
					//ignore first ele == 0 of G_fac->getRowPtr()
					memcpy(rowptr+F_fac->getNbRow()+1, G_fac->getRowPtr()+1, sizeof(int)*G_fac->getNbRow());
					// shift number of elements to take account of already added F_fac elts
					for(faust_unsigned_int j=F_fac->getNbRow()+1;j<F_fac->getNbRow()+G_fac->getNbRow()+1;j++)
						rowptr[j] += F_fac->getRowPtr()[F_fac->getNbRow()];
					// concatened Faust factor
					//				cout << "T_fac_nb_rows:" << T_fac_nb_rows << endl;
					//				cout << "T_fac_nb_cols:" << T_fac_nb_cols << endl;
					//				cout << "nnz:" << T_fac_nnz << endl;
					//				cout << "rowptr=";
					//				for(faust_unsigned_int j=0;j<T_fac_nb_rows+1;j++)
					//					cout << rowptr[j] << ",";
					//				cout << endl;
					//				cout << "colind=";
					//				for(faust_unsigned_int j=0;j<T_fac_nnz;j++)
					//					cout << colind[j] << ",";
					//				cout << endl;
					//				cout << "values=";
					//				for(faust_unsigned_int j=0;j<T_fac_nnz;j++)
					//					cout << values[j] << ",";
					//				cout << endl;
					T_fac = new MatSparse<FPP,Cpu>(T_fac_nnz, T_fac_nb_rows, T_fac_nb_cols, values, rowptr, colind);
					//				cout << "stage 4: ok" << endl;
					//				cout << "T_fac:"<< endl;
					//				T_fac->Display();
					facts[i] = T_fac;
					if(F_inter_fac_allocated)
					{
						delete F_fac;
						F_inter_fac_allocated = ! F_inter_fac_allocated;
					}
					if(G_inter_fac_allocated)
					{
						delete G_fac;
						G_inter_fac_allocated = ! G_inter_fac_allocated;
					}
					delete values;
					delete colind;
					delete rowptr;
				}
				T = new TransformHelper<FPP,Cpu>(facts, 1.0, false, false);
				//			cout << "final stage ok" << endl;
				//			T->display();
				//delete last factors (identity for each Faust)
				if(!F_inter_fac_allocated) delete F_fac;
				if(!F_inter_fac_allocated) delete G_fac;
				return T;
			}

		template<typename FPP>
			TransformHelper<FPP,Cpu>* TransformHelper<FPP,Cpu>::horzcat(const TransformHelper<FPP,Cpu>* G)
			{
				TransformHelper<FPP,Cpu> *Ft, *Gt, *C, *Ct;
				Ft = this->transpose();
				Gt = const_cast<TransformHelper<FPP,Cpu> *>(G)->transpose(); //no harm
				C = Ft->vertcat(Gt);
				Ct = C->transpose();
				delete Ft;
				delete Gt;
				delete C;
				return Ct;
			}
}
