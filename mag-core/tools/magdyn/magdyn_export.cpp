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


static std::string get_str_var(const std::string& var, bool add_brackets = false)
{
	if(var == "")
	{
		return "0";
	}
	else
	{
		if(add_brackets)
			return "(" + var + ")";
		return var;
	}
}



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

	const char* user = std::getenv("USER");
	if(!user)
		user = "";

	ofstr	<< "#\n"
		<< "# Created by Takin/Magdyn\n"
		<< "# URL: https://github.com/ILLGrenoble/takin\n"
		<< "# DOI: https://doi.org/10.5281/zenodo.4117437\n"
		<< "# User: " << user << "\n"
		<< "#\n\n";

	ofstr << "using Sunny\nusing Printf\n\n";


	// --------------------------------------------------------------------
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
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
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
			<< get_str_var(site.pos[0]) << ", "
			<< get_str_var(site.pos[1]) << ", "
			<< get_str_var(site.pos[2]) << " ],"
			<< " # " << site.name << "\n";
	}

	// save as the P1 space group, as we have already performed the symmetry operations
	ofstr << "\t], 1)\n\n";


	ofstr << "# spin magnitudes and magnetic system\n";
	ofstr << "magsys = System(magsites, ( 1, 1, 1 ),\n\t[\n";

	std::size_t site_idx = 1;
	for(const auto &site : m_dyn.GetMagneticSites())
	{
		ofstr << "\t\tSpinInfo("
			<< /*"atom = " <<*/ site_idx << ", "
			<< "S = " << get_str_var(site.spin_mag) << ", "
			<< "g = -[ g_e 0 0; 0 g_e 0; 0 0 g_e ]),"
			<< " # " << site.name << "\n";
		++site_idx;
	}
	ofstr << "\t], :dipole)\n\n";


	ofstr << "# spin directions\n";
	const auto& field = m_dyn.GetExternalField();
	if(field.align_spins)
	{
		// set all spins to field direction
		ofstr << "polarize_spins!(magsys, [ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ])\n";
	}
	else
	{
		// set individual spins
		site_idx = 1;
		for(const auto &site : m_dyn.GetMagneticSites())
		{
			ofstr << "set_dipole!(magsys, [ "
				<< get_str_var(site.spin_dir[0]) << ", "
				<< get_str_var(site.spin_dir[1]) << ", "
				<< get_str_var(site.spin_dir[2]) << " ], "
				<< "( 1, 1, 1, " << site_idx << " ))"
				<< " # " << site.name << "\n";
			++site_idx;
		}
	}

	ofstr << "\n@printf(\"%s\", magsites)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";
	for(const auto& term : m_dyn.GetExchangeTerms())
	{
		std::size_t idx1 = m_dyn.GetMagneticSiteIndex(term.site1) + 1;
		std::size_t idx2 = m_dyn.GetMagneticSiteIndex(term.site2) + 1;

		// TODO: also add general exchange matrix, term.Jgen
		ofstr << "set_exchange!(magsys," << " # " << term.name
			<< "\n\t[\n"
			<< "\t\t" << get_str_var(term.J, true)
			<< "   " << get_str_var(term.dmi[2], true)
			<< "  -" << get_str_var(term.dmi[1], true) << ";"
			<< "\n\t\t-" << get_str_var(term.dmi[2], true)
			<< "   " << get_str_var(term.J, true)
			<< "   " << get_str_var(term.dmi[0], true) << ";"
			<< "\n\t\t" << get_str_var(term.dmi[1], true)
			<< "  -" << get_str_var(term.dmi[0], true)
			<< "   " << get_str_var(term.J, true)
			<< "\n\t], Bond(" << idx1 << ", " << idx2 << ", [ "
			<< get_str_var(term.dist[0]) << ", "
			<< get_str_var(term.dist[1]) << ", "
			<< get_str_var(term.dist[2])
			<< " ]))\n";
	}


	ofstr << "\n# external field\n";
	if(!tl2::equals_0<t_real>(field.mag, g_eps))
	{
		ofstr << "set_external_field!(magsys, [ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ] * " << field.mag
			<< ")\n";
	}
	// --------------------------------------------------------------------

	ofstr << "\n@printf(\"%s\\n\", magsys)\n";

	// --------------------------------------------------------------------
	ofstr << "\n# spin-wave calculation\n";
	ofstr << "calc = SpinWaveTheory(magsys)\n";

	// TODO

	/*ofstr << R"BLOCK(
# output the dispersion and spin-spin correlation
@printf("# %8s %10s %10s %10s %10s\n",
	"h (rlu)", "k (rlu)", "l (rlu)", "E (meV)", "S(Q, E)")
for q_idx in 1:length(qpath)
	Q = qpath[q_idx]
	for e_idx in 1:length(energies[q_idx, :])
		@printf("%10.4f %10.4f %10.4f %10.4f %10.4f\n",
			Q[1], Q[2], Q[3],
			energies[q_idx, e_idx],
			correlations[q_idx, e_idx])
	end
end)BLOCK";*/

	return true;
}
