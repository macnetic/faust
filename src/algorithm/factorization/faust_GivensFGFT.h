
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
			Faust::MatDense<FPP,DEVICE> C;
			Faust::Vect<FPP,DEVICE> C_min_row;
			int* q_candidates; //default IndexType for underlying eigen matrix is int
			vector<FPP2> errs;
			vector<pair<int,int>> coord_choices;
			Faust::MatDense<FPP, DEVICE> Lap;
			Faust::MatDense<FPP, DEVICE> L;
			FPP2 theta;

			// model to blank fact before update (in update_fact())
			vector<int> fact_mod_row_ids;
			vector<int> fact_mod_col_ids;
			vector<FPP> fact_mod_values;

			/**
			 * Matrix pivot and column indices.
			 */
			int p, q;

			/**
			 * Current iteration index (pointing to the current factor to update).
			 */
			unsigned int ite;

			/**
			 * Function pointer to any step of the algorithm.
			 */
			typedef void (GivensFGFT<FPP,DEVICE,FPP2>::*substep_fun)();


			public:
			GivensFGFT(Faust::MatDense<FPP,DEVICE>& Lap, int J);
			~GivensFGFT() {delete[] q_candidates;};

			/**
			 * \brief Algo. main step.
			 */
			void next_step();

			/**
			 * \brief Algo. main loop (facts.size() iterations).
			 */
			void compute_facts();

			/** \brief Algo. step 2.1.
			*/
			void choose_pivot();

			/** \brief Algo. step 2.1.1
			 *
			 */
			void max_L_into_C();

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
			FPP2 get_err(int j) const;

			/**
			 *
			 */
			const MatDense<FPP,DEVICE> get_D() const;

			/**
			 *
			 */
			const MatDense<FPP,DEVICE>& get_L() const ;

			/**
			 *
			 */
			const vector<pair<int,int>>& get_coord_choices() const;

			/**
			 *
			 */
			void get_coord_choice(int j, int& p, int& q) const;

			const MatDense<FPP,DEVICE>& get_Lap() const;

			const vector<MatSparse<FPP,DEVICE>>& get_facts() const;
		};

}
#include "faust_GivensFGFT.hpp"
#endif
