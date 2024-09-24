/**
 * tlibs2 -- magnetic dynamics
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * References:
 *   - (Toth 2015) S. Toth and B. Lake, J. Phys.: Condens. Matter 27 166002 (2015):
 *                 https://doi.org/10.1088/0953-8984/27/16/166002
 *                 https://arxiv.org/abs/1402.6069
 *   - (Heinsdorf 2021) N. Heinsdorf, manual example calculation for a simple
 *                      ferromagnetic case, personal communications, 2021/2022.
 *
 * @desc This file implements the formalism given by (Toth 2015).
 * @desc For further references, see the 'LITERATURE' file.
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_MAGDYN_DISP_H__
#define __TLIBS2_MAGDYN_DISP_H__

#include <algorithm>
#include <numeric>
#include <vector>
#include <string>
#include <thread>
#include <cstdint>

#include <boost/asio.hpp>

#include "../maths.h"
#include "../units.h"
#include "../phys.h"
#include "../algos.h"

#include "magdyn.h"


// --------------------------------------------------------------------
// calculation functions
// --------------------------------------------------------------------

/**
 * unite degenerate energies and their corresponding eigenstates
 */
MAGDYN_TEMPL
MAGDYN_TYPE::EnergiesAndWeights MAGDYN_INST::UniteEnergies(
	const MAGDYN_TYPE::EnergiesAndWeights& energies_and_correlations) const
{
	EnergiesAndWeights new_energies_and_correlations{};
	new_energies_and_correlations.reserve(energies_and_correlations.size());

	for(const EnergyAndWeight& curState : energies_and_correlations)
	{
		const t_real curE = curState.E;

		auto iter = std::find_if(
			new_energies_and_correlations.begin(),
			new_energies_and_correlations.end(),
			[curE, this](const EnergyAndWeight& E_and_S) -> bool
		{
			t_real E = E_and_S.E;
			return tl2::equals<t_real>(E, curE, m_eps);
		});

		if(iter == new_energies_and_correlations.end())
		{
			// energy not yet seen
			new_energies_and_correlations.push_back(curState);
		}
		else
		{
			// energy already seen: add correlation matrices and weights
			iter->S           += curState.S;
			iter->S_perp      += curState.S_perp;
			iter->S_sum       += curState.S_sum;
			iter->S_perp_sum  += curState.S_perp_sum;
			iter->weight      += curState.weight;
			iter->weight_full += curState.weight_full;
		}
	}

	return new_energies_and_correlations;
}



/**
 * get the energies and the spin-correlation at the given momentum
 * (also calculates incommensurate contributions and applies weight factors)
 * @note implements the formalism given by (Toth 2015)
 */
MAGDYN_TEMPL
MAGDYN_TYPE::EnergiesAndWeights
MAGDYN_INST::CalcEnergies(const t_vec_real& Q_rlu,
	bool only_energies) const
{
	EnergiesAndWeights EandWs;

	if(m_calc_H)
	{
		const t_mat H = CalcHamiltonian(Q_rlu);
		EandWs = CalcEnergiesFromHamiltonian(H, Q_rlu, only_energies);
	}

	if(IsIncommensurate())
	{
		// equations (39) and (40) from (Toth 2015)
		const t_mat proj_norm = tl2::convert<t_mat>(
			tl2::projector<t_mat_real, t_vec_real>(m_rotaxis, true));

		t_mat rot_incomm = tl2::unit<t_mat>(3);
		rot_incomm -= s_imag * m_phase_sign * tl2::skewsymmetric<t_mat, t_vec>(m_rotaxis);
		rot_incomm -= proj_norm;
		rot_incomm *= 0.5;

		const t_mat rot_incomm_conj = tl2::conj(rot_incomm);

		EnergiesAndWeights EandWs_p, EandWs_m;

		if(m_calc_Hp)
		{
			const t_mat H_p = CalcHamiltonian(Q_rlu + m_ordering);
			EandWs_p = CalcEnergiesFromHamiltonian(
				H_p, Q_rlu + m_ordering, only_energies);
		}

		if(m_calc_Hm)
		{
			const t_mat H_m = CalcHamiltonian(Q_rlu - m_ordering);
			EandWs_m = CalcEnergiesFromHamiltonian(
				H_m, Q_rlu - m_ordering, only_energies);
		}

		if(!only_energies)
		{
			// formula 40 from (Toth 2015)
			for(EnergyAndWeight& EandW : EandWs)
				EandW.S = EandW.S * proj_norm;
			for(EnergyAndWeight& EandW : EandWs_p)
				EandW.S = EandW.S * rot_incomm;
			for(EnergyAndWeight& EandW : EandWs_m)
				EandW.S = EandW.S * rot_incomm_conj;
		}

		// unite energies and weights
		EandWs.reserve(EandWs.size() + EandWs_p.size() + EandWs_m.size());
		for(EnergyAndWeight& EandW : EandWs_p)
			EandWs.emplace_back(std::move(EandW));
		for(EnergyAndWeight& EandW : EandWs_m)
			EandWs.emplace_back(std::move(EandW));
	}

	if(!only_energies)
		CalcIntensities(Q_rlu, EandWs);

	if(m_unite_degenerate_energies)
		EandWs = UniteEnergies(EandWs);

	if(!only_energies)
		CheckImagWeights(Q_rlu, EandWs);

	return EandWs;
}



MAGDYN_TEMPL
MAGDYN_TYPE::EnergiesAndWeights
MAGDYN_INST::CalcEnergies(t_real h, t_real k, t_real l,
	bool only_energies) const
{
	// momentum transfer
	const t_vec_real Qvec = tl2::create<t_vec_real>({ h, k, l });
	return CalcEnergies(Qvec, only_energies);
}



/**
 * generates the dispersion plot along the given q path
 */
MAGDYN_TEMPL
MAGDYN_TYPE::SofQEs
MAGDYN_INST::CalcDispersion(t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end,
	t_size num_qs, t_size num_threads) const
{
	// determine number of threads
	if(num_threads == 0)
		num_threads = std::max<t_size>(1, std::thread::hardware_concurrency() / 2);

	// thread pool and tasks
	using t_pool = boost::asio::thread_pool;
	using t_task = std::packaged_task<SofQE()>;
	using t_taskptr = std::shared_ptr<t_task>;

	t_pool pool{num_threads};
	std::vector<t_taskptr> tasks;
	tasks.reserve(num_qs);

	// calculate dispersion
	for(t_size i = 0; i < num_qs; ++i)
	{
		auto task = [this, i, num_qs,
			h_start, k_start, l_start,
			h_end, k_end, l_end]() -> SofQE
		{
			// get Q
			const t_real h = std::lerp(h_start, h_end, t_real(i) / t_real(num_qs - 1));
			const t_real k = std::lerp(k_start, k_end, t_real(i) / t_real(num_qs - 1));
			const t_real l = std::lerp(l_start, l_end, t_real(i) / t_real(num_qs - 1));

			// get E and S(Q, E) for this Q
			EnergiesAndWeights E_and_S = CalcEnergies(h, k, l, false);

			SofQE result
			{
				.h = h, .k = k, .l = l,
				.E_and_S = E_and_S
			};

			return result;
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		boost::asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	// collect results
	SofQEs results;
	results.reserve(tasks.size());

	for(auto& task : tasks)
	{
		const SofQE& result = task->get_future().get();
		results.push_back(result);
	}

	return results;
}

// --------------------------------------------------------------------

#endif
