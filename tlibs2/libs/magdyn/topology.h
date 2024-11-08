/**
 * tlibs2 -- magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note Forked on 5-November-2024 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 * @desc For references, see the 'LITERATURE' file.
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

namespace tl2_mag {

/**
 * calculates the berry connection
 * @see equ. 7 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type>
std::vector<t_vec> berry_connections(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon())
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
{
	using t_size = decltype(Q.size());
	constexpr const t_cplx imag{0, 1};

	const t_mat evecs = get_evecs(Q);
	const t_size BANDS = evecs.size1();
	const t_size DIM = Q.size();

	// to ensure correct commutators
	t_mat comm = tl2::unit<t_mat>(BANDS);

	std::vector<t_vec> connections{};
	connections.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
	{
		if(band >= BANDS / 2)
			comm(band, band) = -1;

		connections.emplace_back(tl2::create<t_vec>(DIM));
	}

	for(t_size dim = 0; dim < DIM; ++dim)
	{
		t_vec_real Q1 = Q;
		Q1[dim] += delta;

		// differentiate eigenvector matrix
		t_mat evecs_diff = (get_evecs(Q1) - evecs) / delta;

		t_mat evecs_H = tl2::herm(evecs);
		t_mat C = comm * evecs_H * comm * evecs_diff;

		for(t_size band = 0; band < BANDS; ++band)
			connections[band][dim] = C(band, band) * imag;
	}

	return connections;
}



/**
 * calculates the berry curvature
 * @see equ. 8 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type>
std::vector<t_cplx> berry_curvatures(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon())
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec>
#endif
{
	using t_size = decltype(Q.size());

	const t_mat evecs = get_evecs(Q);
	const t_size BANDS = evecs.size1();
	const t_size DIM = Q.size();

	// only valid in three dimensions
	assert(DIM == 3);

	t_vec_real h = Q, k = Q;
	h[0] += delta;
	k[1] += delta;

	std::vector<t_vec> connections =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, Q, delta);
	std::vector<t_vec> connections_h =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, h, delta);
	std::vector<t_vec> connections_k =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(get_evecs, k, delta);

	std::vector<t_cplx> curvatures{};
	curvatures.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
	{
		// differentiate connection's y component with respect to h
		t_cplx curv1 = (connections_h[band][1] - connections[band][1]) / delta;
		// differentiate connection's x component with respect to k
		t_cplx curv2 = (connections_k[band][0] - connections[band][0]) / delta;

		curvatures.emplace_back(curv1 - curv2);
	}

	return curvatures;
}

}   // namespace tl2_mag



/**
 * get the berry connection for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_vec> MAGDYN_INST::GetBerryConnections(const t_vec_real& Q, t_real delta) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		return S.evec_mat;
	};

	return berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
		get_states, Q, delta);
}



/**
 * get the berry curvature for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_cplx> MAGDYN_INST::GetBerryCurvatures(const t_vec_real& Q, t_real delta) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		return S.evec_mat;
	};

	return berry_curvatures<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
		get_states, Q, delta);
}

// --------------------------------------------------------------------

#endif
