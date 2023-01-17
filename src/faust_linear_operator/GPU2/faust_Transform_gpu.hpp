namespace Faust
{
	//TODO: refactor to generic CPU/GPU code (using if needed a non-member function on Transform<FPP, DEV>)
	template<typename FPP>
		MatDense<FPP,GPU2> Transform<FPP,GPU2>::multiply(const MatDense<FPP,GPU2> &A, const char opThis) /*const*/ //TODO: should be const
		{
			if (size() == 0)
				handleWarning("Transform<FPP,GPU2> : multiply : empty Transform<FPP,GPU2>");

			MatDense<FPP,GPU2> mat(A);


			if (opThis == 'N')
				for (int i=this->size()-1; i >= 0; i--)
					data[i]->multiply(mat, opThis);
			else
				for (int i=0; i < this->size(); i++)
					data[i]->multiply(mat, opThis);

			return mat;
		}
}
