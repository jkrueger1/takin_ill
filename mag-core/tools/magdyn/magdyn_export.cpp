/**
 * magnetic dynamics -- exporting the magnetic structure to other magnon tools
 * @author Tobias Weber <tweber@ill.fr>
 * @date 26-june-2024
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

// these need to be included before all other things on mingw
#include <boost/scope_exit.hpp>
#include <boost/property_tree/xml_parser.hpp>
namespace pt = boost::property_tree;

#include "magdyn.h"

#include <QtCore/QString>

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>

#include "tlibs2/libs/str.h"
#include "tlibs2/libs/units.h"


// precision
extern int g_prec;



/**
 * export magnetic structure to the sunny tool (https://github.com/SunnySuite/Sunny.jl)
 */
void MagDynDlg::ExportToSunny()
{
	QString dirLast = m_sett->value("dir_export_sun", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save As Jl File", dirLast, "jl files (*.jl)");
	if(filename == "")
		return;

	if(ExportToSunny(filename))
		m_sett->setValue("dir_export_sun", QFileInfo(filename).path());
}



/**
 * export magnetic structure to the sunny tool (https://github.com/SunnySuite/Sunny.jl)
 */
bool MagDynDlg::ExportToSunny(const QString& filename)
{
	std::ofstream ofstr(filename.toStdString());
	if(!ofstr)
	{
		QMessageBox::critical(this, "Magnetic Dynamics",
			"Cannot open file for writing.");
		return false;
	}

	ofstr.precision(g_prec);

	ofstr << "using Sunny\nusing Printf\n\n";


	ofstr << "# variables\n";

	// internal variables
	ofstr << "g_e = " << tl2::g_e<t_real> << "\n";

	// user variables
	for(const auto &var : m_dyn.GetVariables())
	{
		ofstr << var.name << " = " << var.value.real();
		if(!tl2::equals_0<t_real>(var.value.imag(), g_eps))
			ofstr << " + var.value.imag()" << "im";
		ofstr << "\n";
	}


	ofstr << "\n# magnetic sites and xtal lattice\n";
	const auto& xtal = m_dyn.GetCrystalLattice();
	ofstr << "magsites = Crystal(\n"
		<< "\tlattice_vectors("
		<< xtal[0] << ", "
		<< xtal[1] << ", "
		<< xtal[2] << ", "
		<< xtal[3] / tl2::pi<t_real> * t_real(180) << ", "
		<< xtal[4] / tl2::pi<t_real> * t_real(180) << ", "
		<< xtal[5] / tl2::pi<t_real> * t_real(180) << "),\n\t[\n";

	ofstr << "\t\t# site list\n";
	for(const auto &site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t\t[ "
			<< site.pos[0] << ", "
			<< site.pos[1] << ", "
			<< site.pos[2] << " ],"
			<< " # " << site.name << "\n";
	}

	// save as the P1 space group, as we have already performed the symmetry operations
	ofstr << "\t], 1)\n\n";


	ofstr << "# spin magnitudes and magnetic system\n";
	ofstr << "magsys = System(magsites, (1, 1, 1),\n\t[\n";

	std::size_t site_idx = 1;
	for(const auto &site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t\tSpinInfo("
			<< /*"atom = " <<*/ site_idx << ", "
			<< "S = " << site.spin_mag << ", "
			<< "g = -[ g_e 0 0; 0 g_e 0; 0 0 g_e ]),"
			<< " # " << site.name << "\n";
		++site_idx;
	}
	ofstr << "\t], :dipole)\n\n";


	ofstr << "# spin directions\n";
	site_idx = 1;
	for(const auto &site : m_dyn.GetMagneticSites())
	{
		ofstr << "set_dipole!(magsys, [ "
			<< site.spin_dir[0] << ", "
			<< site.spin_dir[1] << ", "
			<< site.spin_dir[2] << " ], "
			<< "( 1, 1, 1, " << site_idx << " ))"
			<< " # " << site.name << "\n";
		++site_idx;
	}


	ofstr << "\n@printf(\"%s%s\\n\", magsites, magsys)\n";


	ofstr << "\n# magnetic couplings\n";
	// TODO

	return true;
}
