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

#ifndef __TLIBS2_MATHS_SCALAR_H__
#define __TLIBS2_MATHS_SCALAR_H__

#include <cmath>

#include "../str.h"
#include "../algos.h"
#include "../traits.h"

#include "decls.h"



namespace tl2 {
// ----------------------------------------------------------------------------
// scalar algos
// ----------------------------------------------------------------------------

/**
 * are two scalars equal within an epsilon range?
 */
template<class T>
bool equals(T t1, T t2, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	return std::abs(t1 - t2) <= eps;
}


/**
 * is the given value an integer?
 */
template<class T = double>
bool is_integer(T val, T eps = std::numeric_limits<T>::epsilon())
requires is_scalar<T>
{
	T val_diff = val - std::round(val);
	return tl2::equals<T>(val_diff, T(0), eps);
}


/**
 * get next multiple of the given granularity
 */
template<typename t_num = unsigned int>
t_num next_multiple(t_num num, t_num granularity)
requires is_scalar<t_num>
{
	t_num div = num / granularity;
	bool rest_is_0 = 1;

	if constexpr(std::is_floating_point_v<t_num>)
	{
		div = std::floor(div);
		t_num rest = std::fmod(num, granularity);
		rest_is_0 = tl2::equals(rest, t_num{0});
	}
	else
	{
		t_num rest = num % granularity;
		rest_is_0 = (rest==0);
	}

	return rest_is_0 ? num : (div+1) * granularity;
}


/**
 * mod operation, keeping result positive
 */
template<class t_real>
t_real mod_pos(t_real val, t_real tomod = t_real{2}*pi<t_real>)
requires is_scalar<t_real>
{
	val = std::fmod(val, tomod);
	if(val < t_real(0))
		val += tomod;

	return val;
}


/**
 * are two angles equal within an epsilon range?
 */
template<class T>
bool angle_equals(T t1, T t2,
	T eps = std::numeric_limits<T>::epsilon(),
	T tomod = T{2}*pi<T>)
requires is_scalar<T>
{
	t1 = mod_pos<T>(t1, tomod);
	t2 = mod_pos<T>(t2, tomod);

	return std::abs(t1 - t2) <= eps;
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2,
	typename T::value_type eps = std::numeric_limits<typename T::value_type>::epsilon())
{
	return (std::abs(t1.real() - t2.real()) <= eps) &&
		(std::abs(t1.imag() - t2.imag()) <= eps);
}


/**
 * are two complex numbers equal within an epsilon range?
 */
template<class T> requires is_complex<T>
bool equals(const T& t1, const T& t2, const T& eps)
{
	return tl2::equals<T>(t1, t2, eps.real());
}

// ----------------------------------------------------------------------------

}

// ----------------------------------------------------------------------------
#endif
