/**
 * tlibs2 -- (container-agnostic) math library
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2024
 * @license GPLv3, see 'LICENSE' file
 *
 * @note The present version was forked on 8-Nov-2018 from my privately developed "magtools" project (https://github.com/t-weber/magtools).
 * @note Additional functions forked on 7-Nov-2018 from my privately and TUM-PhD-developed "tlibs" project (https://github.com/t-weber/tlibs).
 * @note Further functions and updates forked on 1-Feb-2021 and 19-Apr-2021 from my privately developed "geo" and "misc" projects (https://github.com/t-weber/geo and https://github.com/t-weber/misc).
 * @note Additional functions forked on 6-Feb-2022 from my privately developed "mathlibs" project (https://github.com/t-weber/mathlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
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

#ifndef __TLIBS2_CXX20_MATH_ALGOS_POLY_H__
#define __TLIBS2_CXX20_MATH_ALGOS_POLY_H__

#include <cmath>
#include <tuple>
#include <vector>
#include <initializer_list>
#include <limits>

#include "../str.h"
#include "../algos.h"
#include "../traits.h"

#include "decls.h"


namespace tl2 {

/**
 * intersection of plane <x|n> = d and line |org> + lam*|dir>
 * @returns [position of intersection, 0: no intersection, 1: intersection, 2: line on plane, line parameter lambda]
 *
 * insert |x> = |org> + lam*|dir> in plane equation:
 * <org|n> + lam*<dir|n> = d
 * lam = (d - <org|n>) / <dir|n>
 *
 * @see http://mathworld.wolfram.com/Line-PlaneIntersection.html
 * @see (Stoecker 1999), Chapter "Analytische Geometrie"
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, int, t_real>
intersect_line_plane(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_vec& planeNorm, t_real plane_d,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	// are line and plane parallel?
	const t_real dir_n = tl2::inner<t_vec>(lineDir, planeNorm);
	if(equals<t_real>(dir_n, 0, eps))
	{
		const t_real org_n = tl2::inner<t_vec>(lineOrg, planeNorm);
		// line on plane?
		if(equals<t_real>(org_n, plane_d, eps))
			return std::make_tuple(t_vec(), 2, t_real(0));
		// no intersection
		return std::make_tuple(t_vec(), 0, t_real(0));
	}

	const t_real org_n = tl2::inner<t_vec>(lineOrg, planeNorm);
	const t_real lam = (plane_d - org_n) / dir_n;

	const t_vec vecInters = lineOrg + lam*lineDir;
	return std::make_tuple(vecInters, 1, lam);
}


/**
 * intersection of plane <x|n> = d and a polygon
 * @return vertices of plane - polygon edge intersections
 */
template<class t_vec, class t_real = typename t_vec::value_type,
	template<class...> class t_cont>
t_cont<t_vec> intersect_plane_poly(
	const t_vec& planeNorm, t_real plane_d,
	const t_cont<t_vec>& polyVerts,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	t_cont<t_vec> edgeInters;

	// intersect with each polygon edge
	for(std::size_t i = 0; i < polyVerts.size(); ++i)
	{
		std::size_t j = (i+1) % polyVerts.size();

		t_vec lineOrg = polyVerts[i];
		t_vec lineDir = polyVerts[j] - lineOrg;

		auto [pos, ty, lam] = intersect_line_plane(
			lineOrg, lineDir, planeNorm, plane_d, eps);

		if(ty == 1 && lam >= 0. && lam < 1.)
		{
			edgeInters.emplace_back(std::move(pos));
		}
		if(ty == 2)
		{
			edgeInters.push_back(polyVerts[i]);
			edgeInters.push_back(polyVerts[j]);
		}
	}

	return edgeInters;
}


/**
 * intersection of a sphere and a line |org> + lam*|dir>
 * @returns vector of intersections
 * insert |x> = |org> + lam*|dir> in sphere equation <x-mid | x-mid> = r^2
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection for the solution.
 */
template<class t_vec, template<class...> class t_cont = std::vector>
requires is_vec<t_vec>
t_cont<t_vec> intersect_line_sphere(
	const t_vec& lineOrg, const t_vec& _lineDir,
	const t_vec& sphereOrg, typename t_vec::value_type sphereRad,
	bool linedir_normalised = false, bool only_segment = false,
	typename t_vec::value_type eps = std::numeric_limits<typename t_vec::value_type>::epsilon())
{
	using T = typename t_vec::value_type;

	t_vec lineDir = _lineDir;
	T lenDir = linedir_normalised ? T(1) : tl2::norm<t_vec>(lineDir);

	if(!linedir_normalised)
		lineDir /= lenDir;

	auto vecDiff = sphereOrg - lineOrg;
	auto proj = tl2::project_scalar<t_vec>(vecDiff, lineDir, true);
	auto rt = proj*proj + sphereRad*sphereRad - tl2::inner<t_vec>(vecDiff, vecDiff);

	// no intersection
	if(rt < T(0))
		return t_cont<t_vec>{};

	// one intersection
	if(equals(rt, T(0), eps))
	{
		T lam = proj/lenDir;
		if(!only_segment || (only_segment && lam >= T(0) && lam < T(1)))
			return t_cont<t_vec>{{ lineOrg + proj*lineDir }};
		return t_cont<t_vec>{};
	}

	// two intersections
	auto val = std::sqrt(rt);
	t_cont<t_vec> inters;
	inters.reserve(2);

	T lam1 = (proj + val)/lenDir;
	T lam2 = (proj - val)/lenDir;

	if(!only_segment || (only_segment && lam1 >= T(0) && lam1 < T(1)))
		inters.emplace_back(lineOrg + (proj + val)*lineDir);
	if(!only_segment || (only_segment && lam2 >= T(0) && lam2 < T(1)))
		inters.emplace_back(lineOrg + (proj - val)*lineDir);

	// sort intersections by x
	std::sort(inters.begin(), inters.end(),
		[](const t_vec& vec1, const t_vec& vec2) -> bool
		{ return vec1[0] < vec2[0]; });

	return inters;
}


/**
 * intersection of a polygon and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& poly)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	// middle point
	const t_vec mid = tl2::mean<t_vec, t_cont>(poly);

	// calculate polygon plane
	const t_vec vec0 = poly[0] - mid;
	const t_vec vec1 = poly[1] - mid;
	t_vec planeNorm = tl2::cross<t_vec>({vec0, vec1});
	planeNorm /= tl2::norm<t_vec>(planeNorm);
	const T planeD = tl2::inner<t_vec>(poly[0], planeNorm);

	// intersection with plane
	auto [vec, intersects, lam] =
		tl2::intersect_line_plane<t_vec>(lineOrg, lineDir, planeNorm, planeD);
	if(intersects != 1)
		return std::make_tuple(t_vec(), false, T(0));

	// is intersection point contained in polygon?
	const t_vec* vecFirst = &(*poly.rbegin());
	for(auto iter=poly.begin(); iter!=poly.end(); std::advance(iter, 1))
	{
		const t_vec* vecSecond = &(*iter);
		const t_vec edge = *vecSecond - *vecFirst;

		// plane through edge
		t_vec edgeNorm = tl2::cross<t_vec>({edge, planeNorm});
		edgeNorm /= tl2::norm<t_vec>(edgeNorm);
		const T edgePlaneD = tl2::inner<t_vec>(*vecFirst, edgeNorm);

		// side of intersection
		const T ptEdgeD = tl2::inner<t_vec>(vec, edgeNorm);

		// outside polygon?
		if(ptEdgeD > edgePlaneD)
			return std::make_tuple(t_vec(), false, T(0));

		vecFirst = vecSecond;
	}

	// intersects with polygon
	return std::make_tuple(vec, true, lam);
}


/**
 * intersection of a polygon (transformed with a matrix) and a line
 * @returns [position of intersection, intersects?, line parameter lambda]
 */
template<class t_vec, class t_mat, template<class ...> class t_cont = std::vector>
std::tuple<t_vec, bool, typename t_vec::value_type>
intersect_line_poly(
	const t_vec& lineOrg, const t_vec& lineDir,
	const t_cont<t_vec>& _poly, const t_mat& mat)
requires is_vec<t_vec> && is_mat<t_mat>
{
	auto poly = _poly;

	// transform each vertex of the polygon
	// TODO: check for homogeneous coordinates!
	for(t_vec& vec : poly)
		vec = mat * vec;

	return tl2::intersect_line_poly<t_vec, t_cont>(lineOrg, lineDir, poly);
}


/**
 * intersection or closest points of lines |org1> + lam1*|dir1> and |org2> + lam2*|dir2>
 * @returns [nearest position 1, nearest position 2, dist, valid?, line parameter 1, line parameter 2]
 *
 * |org1> + lam1*|dir1>  =  |org2> + lam2*|dir2>
 * |org1> - |org2>  =  lam2*|dir2> - lam1*|dir1>
 * |org1> - |org2>  =  (dir2 | -dir1) * |lam2 lam1>
 * (dir2 | -dir1)^T * (|org1> - |org2>)  =  (dir2 | -dir1)^T * (dir2 | -dir1) * |lam2 lam1>
 * |lam2 lam1> = ((dir2 | -dir1)^T * (dir2 | -dir1))^(-1) * (dir2 | -dir1)^T * (|org1> - |org2>)
 *
 * @see https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
 * @see (Stoecker 1999), Chapter "Analytische Geometrie"
 */
template<class t_vec, class t_real = typename t_vec::value_type>
std::tuple<t_vec, t_vec, bool, t_real, t_real, t_real>
intersect_line_line(
	const t_vec& line1Org, const t_vec& line1Dir,
	const t_vec& line2Org, const t_vec& line2Dir,
	t_real eps = std::numeric_limits<t_real>::epsilon())
requires is_vec<t_vec>
{
	const t_vec orgdiff = line1Org - line2Org;

	// direction matrix (symmetric)
	const t_real d11 = inner<t_vec>(line2Dir, line2Dir);
	const t_real d12 = -inner<t_vec>(line2Dir, line1Dir);
	const t_real d22 = inner<t_vec>(line1Dir, line1Dir);

	const t_real d_det = d11*d22 - d12*d12;

	// check if matrix is invertible
	if(equals<t_real>(d_det, 0, eps))
		return std::make_tuple(t_vec(), t_vec(), false, 0, 0, 0);

	// inverse (symmetric)
	const t_real d11_i = d22 / d_det;
	const t_real d12_i = -d12 / d_det;
	const t_real d22_i = d11 / d_det;

	const t_vec v1 = d11_i*line2Dir - d12_i*line1Dir;
	const t_vec v2 = d12_i*line2Dir - d22_i*line1Dir;

	const t_real lam2 = inner<t_vec>(v1, orgdiff);
	const t_real lam1 = inner<t_vec>(v2, orgdiff);

	const t_vec pos1 = line1Org + lam1*line1Dir;
	const t_vec pos2 = line2Org + lam2*line2Dir;
	const t_real dist = norm<t_vec>(pos2 - pos1);

	return std::make_tuple(pos1, pos2, true, dist, lam1, lam2);
}


/**
 * intersection of planes <x|n1> = d1 and <x|n2> = d2
 * @returns line [org, dir, 0: no intersection, 1: intersection, 2: planes coincide]
 *
 * @see http://mathworld.wolfram.com/Plane-PlaneIntersection.html
 */
template<class t_vec>
std::tuple<t_vec, t_vec, int>
	intersect_plane_plane(
	const t_vec& plane1Norm, typename t_vec::value_type plane1_d,
	const t_vec& plane2Norm, typename t_vec::value_type plane2_d)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec lineDir = cross<t_vec>({plane1Norm, plane2Norm});
	const T lenCross = norm<t_vec>(lineDir);

	// planes parallel or coinciding
	if(equals<T>(lenCross, 0))
	{
		const bool bCoincide = equals<T>(plane1_d, plane2_d);
		return std::make_tuple(t_vec(), t_vec(), bCoincide ? 2 : 0);
	}

	lineDir /= lenCross;

	t_vec lineOrg = - cross<t_vec>({plane1Norm, lineDir}) * plane2_d
		+ cross<t_vec>({plane2Norm, lineDir}) * plane1_d;
	lineOrg /= lenCross;

	return std::make_tuple(lineOrg, lineDir, 1);
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 */
template<class t_vec>
t_vec poly_uv_ortho(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_vec<t_vec>
{
	using T = typename t_vec::value_type;

	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;

	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;


	// ----------------------------------------------------
	// orthonormalisation
	const T len12 = norm<t_vec>(vec12);
	const T len13 = norm<t_vec>(vec13);
	const T lenuv12 = norm<t_vec>(uv12);
	const T lenuv13 = norm<t_vec>(uv13);
	auto vecBasis = orthonorm_sys<t_vec, std::initializer_list, std::vector>({vec12, vec13});
	auto uvBasis = orthonorm_sys<t_vec, std::initializer_list, std::vector>({uv12, uv13});
	vec12 = vecBasis[0]*len12; vec13 = vecBasis[1]*len13;
	uv12 = uvBasis[0]*lenuv12; uv13 = uvBasis[1]*lenuv13;
	// ----------------------------------------------------


	const t_vec pt = _pt - vert1;

	// project a point onto a vector and return the fraction along that vector
	auto project_lam = [](const t_vec& vec, const t_vec& vecProj) -> T
	{
		const T len = norm<t_vec>(vecProj);
		const t_vec _vecProj = vecProj / len;
		T lam = inner<t_vec>(_vecProj, vec);
		return lam / len;
	};

	T lam12 = project_lam(pt, vec12);
	T lam13 = project_lam(pt, vec13);

	// uv coordinates at specified point
	const t_vec uv_pt = uv1 + lam12*uv12 + lam13*uv13;
	return uv_pt;
}


/**
 * uv coordinates of a point inside a polygon defined by three vertices
 * (more general version than poly_uv_ortho)
 */
template<class t_mat, class t_vec>
t_vec poly_uv(const t_vec& vert1, const t_vec& vert2, const t_vec& vert3,
	const t_vec& uv1, const t_vec& uv2, const t_vec& uv3,
	const t_vec& _pt)
requires is_mat<t_mat> && is_vec<t_vec>
{
	t_vec vec12 = vert2 - vert1;
	t_vec vec13 = vert3 - vert1;
	t_vec vecnorm = cross<t_vec>({vec12, vec13});

	// basis
	const t_mat basis = create<t_mat, t_vec>({vec12, vec13, vecnorm}, false);

	// reciprocal basis, RECI = REAL^(-T)
	const auto [basisInv, bOk] = inv<t_mat>(basis);
	if(!bOk)
		return zero<t_vec>(uv1.size());

	t_vec pt = _pt - vert1;		// real pt
	pt = basisInv * pt;		// reciprocal pt

	// uv coordinates at specified point
	t_vec uv12 = uv2 - uv1;
	t_vec uv13 = uv3 - uv1;

	// pt has components in common reciprocal basis
	// assumes that both vector and uv coordinates have the same reciprocal basis
	return uv1 + pt[0]*uv12 + pt[1]*uv13;
}
// ----------------------------------------------------------------------------



/**
 * get the normal vector to a polygon
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec get_poly_normal(t_cont<t_vec>& vecPoly)
requires is_vec<t_vec>
{
	using namespace tl2_ops;
	using t_real = typename t_vec::value_type;

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	// face normal
	t_vec vecNormBest = zero<t_vec>(vecCentre.size());
	t_real tBestCross = t_real(0);

	// find non-collinear vectors
	for(std::size_t iVecPoly=1; iVecPoly<vecPoly.size(); ++iVecPoly)
	{
		t_vec vecNorm = tl2::cross<t_vec>({vecPoly[0]-vecCentre, vecPoly[1]-vecCentre});
		t_real tCross = tl2::norm(vecNorm);
		if(tCross > tBestCross)
		{
			tBestCross = tCross;
			vecNormBest = vecNorm;
		}
	}

	// nothing found
	if(vecNormBest.size() < vecCentre.size())
		return t_vec{};

	return vecNormBest;
}


/**
 * sort vertices in a convex polygon
 */
template<class t_vec, template<class...> class t_cont = std::vector>
void sort_poly_verts_norm(t_cont<t_vec>& vecPoly, const t_vec& _vecNorm)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return;

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	const t_vec vecNorm = _vecNorm / tl2::norm(_vecNorm);

	t_vec vec0 = vecPoly[0] - vecCentre;

	std::stable_sort(vecPoly.begin(), vecPoly.end(),
	[&vecCentre, &vec0, &vecNorm](const t_vec& vertex1, const t_vec& vertex2) -> bool
	{
		t_vec vec1 = vertex1 - vecCentre;
		t_vec vec2 = vertex2 - vecCentre;

		return angle(vec0, vec1, &vecNorm) < angle(vec0, vec2, &vecNorm);
	});
}


/**
 * sort vertices in a convex polygon using an absolute centre for determining the normal
 * @returns normal
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec sort_poly_verts(t_cont<t_vec>& vecPoly, const t_vec& vecAbsCentre,
	bool make_norm_perp_to_poly = false)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return t_vec{};

	// line from centre to vertex
	const t_vec vecCentre = mean(vecPoly);
	// face normal
	t_vec vecNorm = vecCentre - vecAbsCentre;

	sort_poly_verts_norm<t_vec, t_cont>(vecPoly, vecNorm);

	if(make_norm_perp_to_poly)
	{
		t_vec normal = get_poly_normal<t_vec, t_cont>(vecPoly);

		// do they point in the same direction?
		if(inner(normal, vecNorm) < 0.)
			normal = -normal;
		vecNorm = normal;
	}

	vecNorm /= tl2::norm(vecNorm);
	return vecNorm;
}


/**
 * sort vertices in a convex polygon determining normal
 * @returns normal
 */
template<class t_vec, template<class...> class t_cont = std::vector>
t_vec sort_poly_verts(t_cont<t_vec>& vecPoly)
requires is_vec<t_vec>
{
	using namespace tl2_ops;

	if(vecPoly.size() <= 1)
		return;

	t_vec vecNorm = get_poly_normal<t_vec, t_cont>(vecPoly);
	sort_poly_verts_norm<t_vec, t_cont>(vecPoly, vecNorm, true);

	vecNorm /= tl2::norm(vecNorm);
	return vecNorm;
}


}

// ----------------------------------------------------------------------------
#endif
