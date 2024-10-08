/**
 * tlibs2 -- magnetic dynamics -- loading and saving
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

#ifndef __TLIBS2_MAGDYN_FILE_H__
#define __TLIBS2_MAGDYN_FILE_H__

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "magdyn.h"



// --------------------------------------------------------------------
// saving of the dispersion data
// --------------------------------------------------------------------
/**
 * generates the dispersion plot along the given Q path
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveDispersion(const std::string& filename,
	t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end,
	t_size num_Qs, t_size num_threads,
	const bool *stop_request) const
{
	std::ofstream ofstr{filename};
	if(!ofstr)
		return false;

	return SaveDispersion(ofstr,
		h_start, k_start, l_start,
		h_end, k_end, l_end, num_Qs,
		num_threads, stop_request);
}



/**
 * generates the dispersion plot along multiple Q paths
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveMultiDispersion(const std::string& filename,
	const std::vector<std::array<t_real, 3>>& Qs,
	t_size num_Qs, t_size num_threads, const bool *stop_request,
	const std::vector<std::string>* Q_names) const
{
	std::ofstream ofstr{filename};
	if(!ofstr)
		return false;

	return SaveMultiDispersion(ofstr,
		Qs, num_Qs, num_threads,
		stop_request, Q_names);
}



/**
 * generates the dispersion plot along the given Q path
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveDispersion(std::ostream& ostr,
	t_real h_start, t_real k_start, t_real l_start,
	t_real h_end, t_real k_end, t_real l_end,
	t_size num_Qs, t_size num_threads, const bool *stop_request) const
{
	ostr.precision(m_prec);
	int field_len = m_prec * 2.5;

	ostr	<< std::setw(field_len) << std::left << "# h" << " "
		<< std::setw(field_len) << std::left << "k" << " "
		<< std::setw(field_len) << std::left << "l" << " "
		<< std::setw(field_len) << std::left << "E" << " "
		<< std::setw(field_len) << std::left << "S(Q,E)" << " "
		<< std::setw(field_len) << std::left << "S_xx" << " "
		<< std::setw(field_len) << std::left << "S_yy" << " "
		<< std::setw(field_len) << std::left << "S_zz"
		<< std::endl;

	SofQEs results = CalcDispersion(h_start, k_start, l_start,
		h_end, k_end, l_end, num_Qs, num_threads, stop_request);

	// print results
	for(const auto& result : results)
	{
		if(stop_request && *stop_request)
			return false;

		for(const EnergyAndWeight& E_and_S : result.E_and_S)
		{
			ostr	<< std::setw(field_len) << std::left << result.h << " "
				<< std::setw(field_len) << std::left << result.k << " "
				<< std::setw(field_len) << std::left << result.l << " "
				<< std::setw(field_len) << E_and_S.E << " "
				<< std::setw(field_len) << E_and_S.weight << " "
				<< std::setw(field_len) << E_and_S.S_perp(0, 0).real() << " "
				<< std::setw(field_len) << E_and_S.S_perp(1, 1).real() << " "
				<< std::setw(field_len) << E_and_S.S_perp(2, 2).real()
				<< std::endl;
		}
	}

	return true;
}



/**
 * generates the dispersion plot along multiple Q paths
 */
MAGDYN_TEMPL
bool MAGDYN_INST::SaveMultiDispersion(std::ostream& ostr,
	const std::vector<std::array<t_real, 3>>& Qs,
	t_size num_Qs, t_size num_threads, const bool *stop_request,
	const std::vector<std::string>* Q_names) const
{
	bool ok = true;
	ostr.precision(m_prec);

	const t_size N = Qs.size();
	for(t_size i = 0; i < N - 1; ++i)
	{
		const std::array<t_real, 3>& Q1 = Qs[i];
		const std::array<t_real, 3>& Q2 = Qs[i + 1];

		ostr << "# ";
		if(Q_names && (*Q_names)[i] != "" && (*Q_names)[i + 1] != "")
		{
			ostr << (*Q_names)[i] << " -> " << (*Q_names)[i + 1]
				<< ": ";
		}
		ostr << "(" << Q1[0] << ", " << Q1[1] << ", " << Q1[2]
			<< ") -> (" << Q2[0] << ", " << Q2[1] << ", " << Q2[2]
			<< ")\n";

		if(!SaveDispersion(ostr, Q1[0], Q1[1], Q1[2],
			Q2[0], Q2[1], Q2[2], num_Qs,
			num_threads, stop_request))
		{
			ok = false;
			break;
		}

		ostr << std::endl;
	}

	return ok;
}
// --------------------------------------------------------------------

#endif
