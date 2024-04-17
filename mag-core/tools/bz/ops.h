/**
 * brillouin zone tool, symop helpers
 * @author Tobias Weber <tweber@ill.fr>
 * @date May-2022
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
 * Copyright (C) 2017-2021  Tobias WEBER (privately developed).
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

#ifndef __BZTOOL_OPS_H__
#define __BZTOOL_OPS_H__

#include <iostream>
#include "tlibs2/libs/phys.h"
#include "tlibs2/libs/str.h"


/**
 * converts a symmetry operation matrix to a string
 */
template<class t_mat>
std::string op_to_str(const t_mat& op)
{
	using namespace tl2_ops;

	std::ostringstream ostr;
	ostr.precision(g_prec);

	for(std::size_t row = 0; row < op.size1(); ++row)
	{
		for(std::size_t col = 0; col < op.size2(); ++col)
		{
			t_real elem = op(row, col);
			tl2::set_eps_0(elem);

			// write periodic values as fractions to avoid
			// calculation errors involving epsilon
			if(tl2::equals<t_real>(elem, 1./3., g_eps))
				ostr << "1/3";
			else if(tl2::equals<t_real>(elem, 2./3., g_eps))
				ostr << "2/3";
			else if(tl2::equals<t_real>(elem, 1./6., g_eps))
				ostr << "1/6";
			else if(tl2::equals<t_real>(elem, 5./6., g_eps))
				ostr << "5/6";
			else if(tl2::equals<t_real>(elem, -1./3., g_eps))
				ostr << "-1/3";
			else if(tl2::equals<t_real>(elem, -2./3., g_eps))
				ostr << "-2/3";
			else if(tl2::equals<t_real>(elem, -1./6., g_eps))
				ostr << "-1/6";
			else if(tl2::equals<t_real>(elem, -5./6., g_eps))
				ostr << "-5/6";
			else
				ostr << elem;

			if(col != op.size2()-1)
				ostr << " ";
		}

		if(row != op.size1()-1)
			ostr << " \n";
	}

	return ostr.str();
}


/**
 * convert a string to a symmetry operation matrix
 */
template<class t_mat>
t_mat str_to_op(const std::string& str)
{
	using namespace tl2_ops;
	t_mat op = tl2::unit<t_mat>(4);

	std::istringstream istr(str);
	for(std::size_t row = 0; row < op.size1(); ++row)
	{
		for(std::size_t col = 0; col < op.size2(); ++col)
		{
			std::string op_str;
			istr >> op_str;

			auto [ok, val] = tl2::eval_expr<std::string, t_real>(op_str);
			op(row, col) = val;

			if(!ok)
			{
				std::cerr << "Error evaluating symop component \""
					<< op_str << "\"" << std::endl;
			}
		}
	}

	return op;
}


/**
 * get the properties of a symmetry operation
 */
template<class t_mat>
std::string get_op_properties(const t_mat& op)
{
	using namespace tl2_ops;
	std::string prop;

	if(tl2::is_unit<t_mat>(op, g_eps))
	{
		if(prop.size())
			prop += ", ";
		prop += "identity";
	}

	if(tl2::hom_is_centring<t_mat>(op, g_eps))
	{
		if(prop.size())
			prop += ", ";
		prop += "centring";
	}

	return prop;
}


#endif
