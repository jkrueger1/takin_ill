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
 * calculates the berry connections
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
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
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
 * calculates the berry curvatures
 * @see equ. 8 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type,
	typename t_size = std::size_t>
std::vector<t_cplx> berry_curvatures(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	const t_vec_real& Q, t_real delta = std::numeric_limits<t_real>::epsilon(),
	t_size dim1 = 0, t_size dim2 = 1)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
{
	const t_mat evecs = get_evecs(Q);
	const t_size BANDS = evecs.size1();

	// only valid in three dimensions
	assert(Q.size() == 3);

	t_vec_real h = Q, k = Q;
	h[dim1] += delta;
	k[dim2] += delta;

	std::vector<t_vec> connections =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
			get_evecs, Q, delta);
	std::vector<t_vec> connections_h =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
			get_evecs, h, delta);
	std::vector<t_vec> connections_k =
		berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
			get_evecs, k, delta);

	std::vector<t_cplx> curvatures{};
	curvatures.reserve(BANDS);

	for(t_size band = 0; band < BANDS; ++band)
	{
		// differentiate connection's y component with respect to h
		t_cplx curv1 = (connections_h[band][dim2] - connections[band][dim2]) / delta;
		// differentiate connection's x component with respect to k
		t_cplx curv2 = (connections_k[band][dim1] - connections[band][dim1]) / delta;

		curvatures.emplace_back(curv1 - curv2);
	}

	return curvatures;
}



/**
 * calculates the chern numbers
 * @see equ. 9 in https://doi.org/10.1146/annurev-conmatphys-031620-104715
 * @see https://en.wikipedia.org/wiki/Berry_connection_and_curvature
 */
template<class t_mat, class t_vec, class t_vec_real,
	typename t_cplx = typename t_vec::value_type,
	typename t_real = typename t_cplx::value_type,
	typename t_size = std::size_t>
std::vector<t_cplx> chern_numbers(
	const std::function<t_mat(const t_vec_real& Q)>& get_evecs,
	t_real bz = 0.5,  // brillouin zone boundary
	t_real delta_diff = std::numeric_limits<t_real>::epsilon(),
	t_real delta_int = std::cbrt(std::numeric_limits<t_real>::epsilon()),
	t_size dim1 = 0, t_size dim2 = 1, bool calc_via_boundary = true)
#ifndef SWIG  // TODO: remove this as soon as swig understands concepts
requires tl2::is_mat<t_mat> && tl2::is_vec<t_vec> && tl2::is_vec<t_vec_real>
#endif
{
	std::vector<t_cplx> chern_nums;

	auto int_boundary = [get_evecs, delta_diff, delta_int, bz, dim1, dim2, &chern_nums](
		t_size dim, t_vec_real& Q, t_real sign)
	{
		for(Q[dim] = -bz; Q[dim] < bz; Q[dim] += delta_int)
		{
			std::vector<t_vec> conns =
				berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
					get_evecs, Q, delta_diff);

			// initialise by resetting chern numbers to zeros
			if(!chern_nums.size())
				chern_nums.resize(conns.size(), t_cplx{});

			// numerically integrate along boundary segment
			for(t_size band = 0; band < conns.size(); ++band)
			{
				chern_nums[band] += conns[band][dim1] * delta_int * sign;
				chern_nums[band] -= conns[band][dim2] * delta_int * sign;
			}
		}
	};

	// calculate via boundary over berry connections
	if(calc_via_boundary)
	{
		// bottom part of boundary
		t_vec_real Q = tl2::zero<t_vec_real>(3);
		Q[dim2] -= bz;
		int_boundary(dim1, Q, 1.);

		// top part of boundary
		Q = tl2::zero<t_vec_real>(3);
		Q[dim2] += bz;
		int_boundary(dim1, Q, -1.);

		// left part of boundary
		Q = tl2::zero<t_vec_real>(3);
		Q[dim1] -= bz;
		int_boundary(dim2, Q, 1.);

		// right part of boundary
		Q = tl2::zero<t_vec_real>(3);
		Q[dim1] += bz;
		int_boundary(dim2, Q, -1.);
	}

	// calculate via area over berry curvatures
	else
	{
		t_vec_real Q = tl2::zero<t_vec_real>(3);

		for(Q[dim1] = -bz; Q[dim1] < bz; Q[dim1] += delta_int)
		for(Q[dim2] = -bz; Q[dim2] < bz; Q[dim2] += delta_int)
		{
			std::vector<t_cplx> curvs = berry_curvatures<
				t_mat, t_vec, t_vec_real, t_cplx, t_real, t_size>(
					get_evecs, Q, delta_diff, dim1, dim2);

			// initialise by resetting chern numbers to zeros
			if(!chern_nums.size())
				chern_nums.resize(curvs.size(), t_cplx{});

			// numerically integrate the brillouin zone area for each band
			for(t_size band = 0; band < curvs.size(); ++band)
				chern_nums[band] += curvs[band] * delta_int * delta_int;
		}

		for(t_cplx& num : chern_nums)
			num /= t_real(2) * tl2::pi<t_real>;
	}

	// should be integers
	return chern_nums;
}

}   // namespace tl2_mag



/**
 * get the berry connection for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_vec> MAGDYN_INST::CalcBerryConnections(
	const t_vec_real& Q, t_real delta,
	const std::vector<t_size>* perm) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this, perm](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		t_mat M = S.evec_mat;

		if(perm)
			M = tl2::reorder_cols<t_mat, t_vec>(M, *perm);
		return M;
	};

	return berry_connections<t_mat, t_vec, t_vec_real, t_cplx, t_real>(
		get_states, Q, delta);
}



/**
 * get the berry curvature for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_cplx> MAGDYN_INST::CalcBerryCurvatures(
	const t_vec_real& Q, t_real delta,
	const std::vector<t_size>* perm,
	t_size dim1, t_size dim2) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this, perm](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		t_mat M = S.evec_mat;

		if(perm)
			M = tl2::reorder_cols<t_mat, t_vec>(M, *perm);
		return M;
	};

	return berry_curvatures<
		t_mat, t_vec, t_vec_real,t_cplx, t_real, t_size>(
			get_states, Q, delta, dim1, dim2);
}



/**
 * get the berry curvature for each magnon band
 */
MAGDYN_TEMPL
std::vector<t_cplx> MAGDYN_INST::CalcChernNumbers(
	t_real bz, t_real delta_diff, t_real delta_int,
	t_size dim1, t_size dim2) const
{
	//SetUniteDegenerateEnergies(false);

	// get eigenstates at specific Q
	auto get_states = [this](const t_vec_real& Q) -> t_mat
	{
		SofQE S = CalcEnergies(Q, false);
		t_mat M = S.evec_mat;

		//if(perm)
		//	M = tl2::reorder_cols<t_mat, t_vec>(M, *perm);
		return M;
	};

	bool calc_via_boundary = true;
	return chern_numbers<
		t_mat, t_vec, t_vec_real, t_cplx, t_real, t_size>(
			get_states, bz, delta_diff, delta_int,
			dim1, dim2, calc_via_boundary);
}

// --------------------------------------------------------------------

#endif

