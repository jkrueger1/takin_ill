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

#ifndef __TLIBS2_MATHS_HELPERS_H__
#define __TLIBS2_MATHS_HELPERS_H__

#include <cmath>
#include <vector>
#include <limits>
#include <random>

#include "../str.h"
#include "../algos.h"
#include "../traits.h"

#include "decls.h"


namespace tl2 {
// ----------------------------------------------------------------------------
// helpers and constants
// ----------------------------------------------------------------------------

// constants
#ifdef __TLIBS2_USE_NUMBERS__
	template<typename T=double> constexpr T golden{std::numbers::phi_v<T>};
	template<typename T=double> constexpr T pi{std::numbers::pi_v<T>};
#else
	// see: https://en.wikipedia.org/wiki/Golden_ratio
	template<typename T=double> constexpr T golden{1.618033988749895};
	template<typename T=double> constexpr T pi{M_PI};
#endif

template<typename INT=int> bool is_even(INT i) { return (i%2 == 0); }
template<typename INT=int> bool is_odd(INT i) { return !is_even<INT>(i); }

template<class T=double> constexpr T r2d(T rad) { return rad/pi<T>*T(180); }	// rad -> deg
template<class T=double> constexpr T d2r(T deg) { return deg/T(180)*pi<T>; }	// deg -> rad
template<class T=double> constexpr T r2m(T rad) { return rad/pi<T>*T(180*60); }	// rad -> min
template<class T=double> constexpr T m2r(T min) { return min/T(180*60)*pi<T>; }	// min -> rad

/**
 * Gaussian around 0: f(x) = exp(-1/2 * (x/sig)^2)
 * at hwhm: f(x_hwhm) = 1/2
 *          exp(-1/2 * (x_hwhm/sig)^2) = 1/2
 *          -1/2 * (x_hwhm/sig)^2 = ln(1/2)
 *          (x_hwhm/sig)^2 = -2*ln(1/2)
 *          x_hwhm^2 = sig^2 * 2*ln(2)
 */
template<class T=double> static constexpr T SIGMA2FWHM = T(2)*std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static constexpr T SIGMA2HWHM = std::sqrt(T(2)*std::log(T(2)));
template<class T=double> static constexpr T FWHM2SIGMA = T(1)/SIGMA2FWHM<T>;
template<class T=double> static constexpr T HWHM2SIGMA = T(1)/SIGMA2HWHM<T>;


template<typename T> T sign(T t)
{
	if(t<0.) return -T(1);
	return T(1);
}


template<typename T> T cot(T t)
{
	return std::tan(T(0.5)*pi<T> - t);
}


template<typename T> T coth(T t)
{
	return T(1) / std::tanh(t);
}


template<typename T=double>
T log(T tbase, T tval)
{
	return T(std::log(tval)/std::log(tbase));
}


template<typename T=double>
T nextpow(T tbase, T tval)
{
	return T(std::pow(tbase, std::ceil(log(tbase, tval))));
}


/**
 * unsigned angle between two vectors
 * <q1|q2> / (|q1| |q2|) = cos(alpha)
 */
template<class t_vec>
typename t_vec::value_type angle_unsigned(const t_vec& q1, const t_vec& q2)
requires is_basic_vec<t_vec>
{
	using t_real = typename t_vec::value_type;

	if(q1.size() != q2.size())
		return t_real(0);

	t_real dot = t_real(0);
	t_real len1 = t_real(0);
	t_real len2 = t_real(0);

	for(std::size_t i=0; i<q1.size(); ++i)
	{
		dot += q1[i]*q2[i];

		len1 += q1[i]*q1[i];
		len2 += q2[i]*q2[i];
	}

	len1 = std::sqrt(len1);
	len2 = std::sqrt(len2);

	dot /= len1;
	dot /= len2;

	return std::acos(dot);
}


/**
 * x = 0..1, y = 0..1
 * @see https://en.wikipedia.org/wiki/Bilinear_interpolation
 */
template<typename T=double>
T bilinear_interp(T x0y0, T x1y0, T x0y1, T x1y1, T x, T y)
{
	T top = std::lerp(x0y1, x1y1, x);
	T bottom = std::lerp(x0y0, x1y0, x);

	return std::lerp(bottom, top, y);
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> linspace(const T& tmin, const T& tmax, std::size_t iNum)
{
	t_vec<T> vec;
	vec.reserve(iNum);

	if(iNum == 1)
	{
		// if just one point is requested, use the lower limit
		vec.push_back(tmin);
		return vec;
	}

	for(std::size_t i=0; i<iNum; ++i)
		vec.push_back(std::lerp(tmin, tmax, REAL(i)/REAL(iNum-1)));
	return vec;
}


template<typename T=double, typename REAL=double,
template<class...> class t_vec = std::vector>
t_vec<T> logspace(const T& tmin, const T& tmax, std::size_t iNum, T tBase=T(10))
{
	t_vec<T> vec = linspace<T, REAL>(tmin, tmax, iNum);

	for(T& t : vec)
		t = std::pow(tBase, t);
	return vec;
}


template<class T>
bool is_in_range(T val, T centre, T pm)
{
	pm = std::abs(pm);

	if(val < centre-pm) return false;
	if(val > centre+pm) return false;
	return true;
}


/**
 * point contained in linear range?
 */
template<class T = double>
bool is_in_linear_range(T dStart, T dStop, T dPoint)
{
	if(dStop < dStart)
		std::swap(dStart, dStop);

	return (dPoint >= dStart) && (dPoint <= dStop);
}


/**
 * angle contained in angular range?
 */
template<class T = double>
bool is_in_angular_range(T dStart, T dRange, T dAngle)
{
	if(dStart < T(0)) dStart += T(2)*pi<T>;
	if(dAngle < T(0)) dAngle += T(2)*pi<T>;

	dStart = std::fmod(dStart, T(2)*pi<T>);
	dAngle = std::fmod(dAngle, T(2)*pi<T>);

	T dStop = dStart + dRange;


	// if the end point is contained in the circular range
	if(dStop < T(2)*pi<T>)
	{
		return is_in_linear_range<T>(dStart, dStop, dAngle);
	}
	// else end point wraps around
	else
	{
		return is_in_linear_range<T>(dStart, T(2)*pi<T>, dAngle) ||
			is_in_linear_range<T>(T(0), dRange-(T(2)*pi<T>-dStart), dAngle);
	}
}


/**
 * get a random number in the given range
 */
template<class t_num>
t_num get_rand(t_num min=1, t_num max=-1)
{
	static std::mt19937 rng{std::random_device{}()};

	if(max <= min)
	{
		min = std::numeric_limits<t_num>::lowest() / 10.;
		max = std::numeric_limits<t_num>::max() / 10.;
	}

	if constexpr(std::is_integral_v<t_num>)
		return std::uniform_int_distribution<t_num>(min, max)(rng);
	else
		return std::uniform_real_distribution<t_num>(min, max)(rng);
}
// ----------------------------------------------------------------------------


}

// ----------------------------------------------------------------------------
#endif
