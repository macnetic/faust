/****************************************************************************/
/*                              Description:                                */
/*  For more information on the FAuST Project, please visit the website     */
/*  of the project : <http://faust.inria.fr>                         */
/*                                                                          */
/*                              License:                                    */
/*  Copyright (2019):    Hakim Hadj-Djilani, Nicolas Bellot, Adrien Leman, Thomas Gautrais,      */
/*                      Luc Le Magoarou, Remi Gribonval                     */
/*                      INRIA Rennes, FRANCE                                */
/*                      http://www.inria.fr/                                */
/*                                                                          */
/*  The FAuST Toolbox is distributed under the terms of the GNU Affero      */
/*  General Public License.                                                 */
/*  This program is free software: you can redistribute it and/or modify    */
/*  it under the terms of the GNU Affero General Public License as          */
/*  published by the Free Software Foundation.                              */
/*                                                                          */
/*  This program is distributed in the hope that it will be useful, but     */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of              */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                    */
/*  See the GNU Affero General Public License for more details.             */
/*                                                                          */
/*  You should have received a copy of the GNU Affero General Public        */
/*  License along with this program.                                        */
/*  If not, see <http://www.gnu.org/licenses/>.                             */
/*                                                                          */
/*                             Contacts:                                    */
/*      Hakim H. hakim.hadj-djilani@inria.fr                                */
/*      Nicolas Bellot  : nicolas.bellot@inria.fr                           */
/*      Adrien Leman    : adrien.leman@inria.fr                             */
/*      Thomas Gautrais : thomas.gautrais@inria.fr                          */
/*      Luc Le Magoarou : luc.le-magoarou@inria.fr                          */
/*      Remi Gribonval  : remi.gribonval@inria.fr                           */
/*                                                                          */
/*                              References:                                 */
/*  [1] Le Magoarou L. and Gribonval R., "Flexible multi-layer sparse       */
/*  approximations of matrices and applications", Journal of Selected       */
/*  Topics in Signal Processing, 2016.                                      */
/*  <https://hal.archives-ouvertes.fr/hal-01167948v1>                       */
/****************************************************************************/
#ifndef __FAUST_PROX_CU_HPP__
#define __FAUST_PROX_CU_HPP__


// const char * interface_prox_name="prox : ";

namespace Faust
{
	template<typename FPP>
		void prox_sp(MatDense<FPP,GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
		{
			if(pure_gpu)
			{
				M.prox_sp(k, normalized, pos);
			}
			else
			{
				MatDense<FPP,Cpu> cpuM = M.tocpu();
				prox_sp(cpuM, k, normalized, pos);
				M = cpuM;
			}
		}

	template<typename FPP>
		void prox_sp_pos(MatDense<FPP, GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false */)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_sp_pos(cpuM, k, normalized, pos);
			M = cpuM;
		}

	template<typename FPP>
		void prox_spcol(MatDense<FPP,GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
		{
			if(pure_gpu)
			{
				M.prox_spcol(k, normalized, pos);
			}
			else
			{
				MatDense<FPP,Cpu> cpuM = M.tocpu();
				prox_spcol(cpuM, k, normalized, pos);
				M = cpuM;
			}
		}

	template<typename FPP>
		void prox_splin(MatDense<FPP,GPU2> & M,faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
		{
			if(pure_gpu)
			{
				M.prox_splin(k, normalized, pos);
			}
			else
			{
				MatDense<FPP,Cpu> cpuM = M.tocpu();
				prox_splin(cpuM, k, normalized, pos);
				M = cpuM;
			}
		}

	template<typename FPP>
		void prox_splincol(MatDense<FPP,GPU2> &M,faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_splincol(cpuM, k, normalized, pos);
			M = cpuM;
		}

	template<typename FPP>
		void prox_supp(MatDense<FPP,GPU2> & M, const MatDense<FPP,GPU2> & supp, const bool normalized/*=true*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			MatDense<FPP,Cpu> cpu_supp = supp.tocpu();
			prox_supp(cpuM, cpu_supp, normalized, pos);
			M = cpuM;
		}

	template<typename FPP>
		void prox_const(MatDense<FPP,GPU2> & M, const MatDense<FPP,GPU2> & supp, const bool normalized/*=false*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			MatDense<FPP,Cpu> cpu_supp = supp.tocpu();
			prox_const(cpuM, cpu_supp, normalized, pos);
			M = cpuM;
		}

	template<typename FPP>
		void prox_id(MatDense<FPP,GPU2> & M, const bool normalized/*=false*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_id(cpuM, normalized, pos);
			M = cpuM;
		}

	template<typename FPP, typename FPP2>
		void prox_normcol(MatDense<FPP,GPU2> & M,FPP2 s, const bool normalized/*=false*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_normcol(cpuM, s, normalized, pos);
			M = cpuM;
		}
	template<typename FPP, typename FPP2>
		void prox_normlin(MatDense<FPP,GPU2> & M,FPP2 s, const bool normalized/*=false*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_normlin(cpuM, s, normalized, pos);
			M = cpuM;
		}
	//	template<typename FPP>
	//		void prox_blkdiag(MatDense<FPP,GPU2> & M,faust_unsigned_int k)
	//	{
	//		MatDense<FPP,Cpu> cpuM = M.tocpu();
	//		prox_gpu(cpuM, k);
	//		M = cpuM;
	//	}
	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,GPU2> & M, std::vector<faust_unsigned_int>& m_vec, std::vector<faust_unsigned_int>& n_vec, const bool normalized /*= false*/, const bool pos /*= false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_blockdiag(cpuM, m_vec, n_vec, normalized, pos);
			M = cpuM;
		}

	template<typename FPP>
		void prox_blockdiag(MatDense<FPP,GPU2> & M, MatDense<FPP,GPU2> mn_vec, const bool normalized/*=true*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			MatDense<FPP,Cpu> cpu_mn_vec = mn_vec.tocpu();
			prox_blockdiag(cpuM, cpu_mn_vec, normalized, pos);
			M = cpuM;
		}

	template<typename FPP> void pre_prox_pos(MatDense<FPP,GPU2> & M)
	{
		MatDense<FPP,Cpu> cpuM = M.tocpu();
		pre_prox_pos(cpuM);
		M = cpuM;
	}

	template<typename FPP> void prox_hankel(MatDense<FPP, GPU2> & M, const bool normalized /*= true*/, const bool pos /*= false*/)
	{
		MatDense<FPP,Cpu> cpuM = M.tocpu();
		prox_hankel(cpuM, normalized, pos);
		M = cpuM;
	}

	template<typename FPP> void prox_toeplitz(MatDense<FPP, GPU2> & M, const bool normalized /*= true*/, const bool pos /*= false*/)
	{
		MatDense<FPP,Cpu> cpuM = M.tocpu();
		prox_toeplitz(cpuM, normalized, pos);
		M = cpuM;
	}

	template<typename FPP> void prox_circ(MatDense<FPP, GPU2> & M, const bool normalized /*= true*/, const bool pos /*= false*/)
	{
		MatDense<FPP,Cpu> cpuM = M.tocpu();
		prox_circ(cpuM, normalized, pos);
		M = cpuM;
	}

	template<typename FPP> void prox_anticirc(MatDense<FPP, GPU2> & M, const bool normalized /*= true*/, const bool pos /*= false*/)
	{
		MatDense<FPP,Cpu> cpuM = M.tocpu();
		prox_anticirc(cpuM, normalized, pos);
		M = cpuM;
	}

	template<typename FPP>
		void prox_skperm(MatDense<FPP, GPU2> & M,const unsigned int k,  const bool normalized/*=true*/, const bool pos/*=false*/)
		{
			MatDense<FPP,Cpu> cpuM = M.tocpu();
			prox_skperm(cpuM, k, normalized, pos);
			M = cpuM;
		}

  template<typename FPP>
  void prox_triu_sp(MatDense<FPP,GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
  {
    prox_tri_sp(M, k, true, normalized, pos, pure_gpu);
  }

  template<typename FPP>
  void prox_tril_sp(MatDense<FPP,GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
  {
    prox_tri_sp(M, k, false, normalized, pos, pure_gpu);
  }

  template<typename FPP>
  void prox_tri_sp(MatDense<FPP,GPU2> & M, faust_unsigned_int k, bool upper, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
  {
    if(pure_gpu)
//      M.prox_tri_sp(k, upper, normalized, pos);
		throw std::runtime_error("prox_tri_sp is not implemented on GPU");
    else
      {
	MatDense<FPP,Cpu> cpuM = M.tocpu();
	prox_tri_sp(cpuM, k, upper, normalized, pos);
	M = cpuM;
      }
  }

  template<typename FPP>
  void prox_symm_sp(MatDense<FPP,GPU2> & M, faust_unsigned_int k, const bool normalized/*=true*/, const bool pos/*=false*/, const bool pure_gpu/*=true*/)
  {
    if(pure_gpu)
//      M.prox_symm_sp(k, normalized, pos);
		throw std::runtime_error("prox_symm_sp is not implemented on GPU");
    else
      {
	MatDense<FPP,Cpu> cpuM = M.tocpu();
	prox_symm_sp(cpuM, k, normalized, pos);
	M = cpuM;
      }
  }

}
#endif


