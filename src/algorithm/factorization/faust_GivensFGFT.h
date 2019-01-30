
#ifndef __GIVENS_FGFT__
#define __GIVENS_FGFT__

#include "faust_constant.h"
#include "faust_MatSparse.h"
#include "faust_MatDense.h"
#include <vector>

namespace Faust {

	template<typename FPP, Device DEVICE, typename FPP2 = float>
		class GivensFGFT {

			vector<Faust::MatSparse<FPP,DEVICE>> facts;
			Faust::MatSparse<FPP,DEVICE> D;
			Faust::MatSparse<FPP,DEVICE> C;
			vector<FPP2> errs;
			vector<pair<faust_unsigned_int,faust_unsigned_int>> coord_choices;
			Faust::MatDense<FPP, DEVICE> Lap;
			Faust::MatDense<FPP, DEVICE> L;

			faust_unsigned_int p, q;

			public:
			GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, faust_unsigned_int J);
			/**
			 * \brief Algo. main step.
			 */
			void next_step();
			void compute_facts();
			/** \brief Algo. step 2.1.
			*/
			void choose_pivot();
			/** \brief Algo. step 2.1.1
			 *
			 */
			void sort_L_in_C();
			/** \brief Algo. step 2.2.
			*/
			void calc_theta();
			/**
			 * \brief Algo. step 2.3.
			 */
			void update_fact();

			/**
			 * \brief Algo. step 2.4
			 */
			void update_L();

			/**
			 * \brief Algo. step 2.5.
			 */
			void update_D();

			/**
			 * \brief Algo. step 2.6.
			 */
			void update_err();

			/**
			 *
			 */
			const vector<FPP2>& get_errs() const;

			/**
			 *
			 */
			FPP2 get_err(faust_unsigned_int j) const;

			/**
			 *
			 */
			const MatDense<FPP,DEVICE>& get_D() const;

			/**
			 *
			 */
			const MatDense<FPP,DEVICE>& get_L() const ;

			/**
			 *
			 */
			const vector<pair<faust_unsigned_int,faust_unsigned_int>>& get_coord_choices() const;

			/**
			 *
			 */
			void get_coord_choice(faust_unsigned_int j, faust_unsigned_int& p, faust_unsigned_int& q) const;

			const MatDense<FPP,DEVICE>& get_Lap() const;

			const vector<MatSparse<FPP,DEVICE>>& get_facts() const;
		};

}
#include "faust_GivensFGFT.hpp"
#endif
