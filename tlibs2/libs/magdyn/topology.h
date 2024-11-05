/**
 * tlibs2 -- magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
 * @license GPLv3, see 'LICENSE' file
 *
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

#ifndef __TLIBS2_MAGDYN_TOPO_H__
#define __TLIBS2_MAGDYN_TOPO_H__

#include "../maths.h"

#include "magdyn.h"



// --------------------------------------------------------------------
// topological calculations
// --------------------------------------------------------------------

/**
 * get the berry connection for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_vec> MAGDYN_INST::GetBerryConnections(const t_vec_real& /*Q_start*/, t_real /*delta*/) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		return S.evec_mat;
	};


	std::vector<t_vec> conn{};

	// TODO

	return conn;
}



/**
 * get the berry curvature for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_cplx> MAGDYN_INST::GetBerryCurvatures(const t_vec_real& /*Q_start*/, t_real /*delta*/) const
{
	//SetUniteDegenerateEnergies(false);

	std::vector<t_cplx> curv{};

	// TODO

	return curv;
}

// --------------------------------------------------------------------

#endif
