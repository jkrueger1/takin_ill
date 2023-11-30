/**
 * magnon dynamics -- type definitions and setting variables
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * "misc" project
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

#ifndef __MAGDYN_DEFS__
#define __MAGDYN_DEFS__

#include <QtCore/QString>

#include <string>
#include <array>
#include <variant>
#include <optional>

#include "tlibs2/libs/maths.h"
#include "tlibs2/libs/qt/gl.h"



// ----------------------------------------------------------------------------
// type definitions
// ----------------------------------------------------------------------------
using t_size = std::size_t;
using t_real = /*float*/ double;
using t_cplx = std::complex<t_real>;

using t_vec_real = tl2::vec<t_real, std::vector>;
using t_mat_real = tl2::mat<t_real, std::vector>;

using t_vec = tl2::vec<t_cplx, std::vector>;
using t_mat = tl2::mat<t_cplx, std::vector>;

using t_real_gl = tl2::t_real_gl;
using t_vec2_gl = tl2::t_vec2_gl;
using t_vec3_gl = tl2::t_vec3_gl;
using t_vec_gl = tl2::t_vec_gl;
using t_mat_gl = tl2::t_mat_gl;
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// maximum number of recent files
extern unsigned int g_maxnum_recents;

// number precisions
extern int g_prec, g_prec_gui;

// epsilons
extern t_real g_eps;

// optional features
extern int g_allow_ortho_spin, g_allow_general_J;

// gui theme and font
extern QString g_theme, g_font;

// use native menubar and dialogs?
extern int g_use_native_menubar, g_use_native_dialogs;

// sorting of sites
extern int g_allow_sorting_sites;

// transfer the setting from the takin core program
void get_settings_from_takin_core();
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// variables register
// ----------------------------------------------------------------------------
#include "settings.h"

constexpr std::array<SettingsVariable, 7> g_settingsvariables
{{
	// epsilons and precisions
	{
		.description = "Calculation epsilon.",
		.key = "eps",
		.value = &g_eps,
	},
	{
		.description = "Number precision.",
		.key = "prec",
		.value = &g_prec,
	},
	{
		.description = "GUI number precision.",
		.key = "prec_gui",
		.value = &g_prec_gui,
	},

	// file options
	{
		.description = "Maximum number of recent files.",
		.key = "maxnum_recents",
		.value = &g_maxnum_recents,
	},

	// optional features
	{
		.description = "Allow setting of orthogonal spins.",
		.key = "ortho_spins",
		.value = &g_allow_ortho_spin,
		.editor = SettingsVariableEditor::YESNO,
	},
	{
		.description = "Allow setting of general exchange matrix J.",
		.key = "allow_gen_J",
		.value = &g_allow_general_J,
		.editor = SettingsVariableEditor::YESNO,
	},

	// sorting of sites
	{
		.description = "Allow sorting of sites (careful: index depends on order!).",
		.key = "allow_sorting_sites",
		.value = &g_allow_sorting_sites,
		.editor = SettingsVariableEditor::YESNO,
	},
}};
// ----------------------------------------------------------------------------


#endif
