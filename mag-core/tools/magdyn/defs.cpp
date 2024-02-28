/**
 * magnon dynamics -- type definitions and setting variables
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2023
 * @license GPLv3, see 'LICENSE' file
 * @desc The present version was forked on 28-Dec-2018 from my privately developed "misc" project (https://github.com/t-weber/misc).
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#include <QtCore/QSettings>

#include "defs.h"



// ----------------------------------------------------------------------------
// global settings variables
// ----------------------------------------------------------------------------
// maximum number of recent files
unsigned int g_maxnum_recents = 16;

// epsilons and precisions
int g_prec = 6;
int g_prec_gui = 3;
t_real g_eps = 1e-6;

// optional features
int g_allow_ortho_spin = 0;
int g_allow_general_J = 1;

// gui theme
QString g_theme = "Fusion";

// gui font
QString g_font = "";

// use native menu bar?
int g_use_native_menubar = 0;

// use native dialogs?
int g_use_native_dialogs = 0;
// ----------------------------------------------------------------------------



/**
 * transfer the setting from the takin core program
 */
void get_settings_from_takin_core()
{
	QSettings sett_core("takin", "core");

	if(sett_core.contains("main/font_gen"))
	{
		g_font = sett_core.value("main/font_gen").toString();
	}

	if(sett_core.contains("main/prec"))
	{
		g_prec = sett_core.value("main/prec").toInt();
		g_eps = std::pow(t_real(10), -t_real(g_prec));
	}

	if(sett_core.contains("main/prec_gfx"))
	{
		g_prec_gui = sett_core.value("main/prec_gfx").toInt();
	}

	if(sett_core.contains("main/gui_style_value"))
	{
		g_theme = sett_core.value("main/gui_style_value").toString();
		//std::cout << g_theme.toStdString() << std::endl;
	}
}
