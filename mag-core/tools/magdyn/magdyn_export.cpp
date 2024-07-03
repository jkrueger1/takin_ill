/**
 * magnetic dynamics -- exporting the magnetic structure to other magnon tools
 * @author Tobias Weber <tweber@ill.fr>
 * @date 26-june-2024
 * @license GPLv3, see 'LICENSE' file
 *
 * ----------------------------------------------------------------------------
 * mag-core (part of the Takin software suite)
 * Copyright (C) 2018-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
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
 * export the magnetic structure to the sunny tool
 *   (https://github.com/SunnySuite/Sunny.jl)
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
 * export the magnetic structure to the sunny tool
 *   (https://github.com/SunnySuite/Sunny.jl)
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
	ofstr << "@printf(stderr, \"Setting up magnetic sites...\\n\")\n";

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
	// (you can also manually set the crystal's space group and delete all
	//  symmetry-equivalent positions and couplings in the generated file)
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

	ofstr << "\n@printf(stderr, \"%s\", magsites)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# magnetic couplings\n";
	ofstr << "@printf(stderr, \"Setting up magnetic couplings...\\n\")\n";

	for(const auto& term : m_dyn.GetExchangeTerms())
	{
		std::size_t idx1 = m_dyn.GetMagneticSiteIndex(term.site1) + 1;
		std::size_t idx2 = m_dyn.GetMagneticSiteIndex(term.site2) + 1;

		ofstr	<< "set_exchange!(magsys," << " # " << term.name
			<< "\n\t[\n"
			<< "\t\t" << get_str_var(term.J, true)             // 0,0
			<< "   " << get_str_var(term.dmi[2], true)         // 0,1
			<< "  -" << get_str_var(term.dmi[1], true) << ";"  // 0,2
			<< "\n\t\t-" << get_str_var(term.dmi[2], true)     // 1,0
			<< "   " << get_str_var(term.J, true)              // 1,1
			<< "   " << get_str_var(term.dmi[0], true) << ";"  // 1,2
			<< "\n\t\t" << get_str_var(term.dmi[1], true)      // 2,0
			<< "  -" << get_str_var(term.dmi[0], true)         // 2,1
			<< "   " << get_str_var(term.J, true)              // 2,2
			<< "\n\t]";

		if(!tl2::equals_0(term.Jgen_calc))
		{
			ofstr	<< " +\n\t[\n"
				<< "\t\t" << get_str_var(term.Jgen[0][0], true)
				<< "  " << get_str_var(term.Jgen[0][1], true)
				<< "  " << get_str_var(term.Jgen[0][2], true) << ";"
				<< "\n\t\t" << get_str_var(term.Jgen[1][0], true)
				<< "  " << get_str_var(term.Jgen[1][1], true)
				<< "  " << get_str_var(term.Jgen[1][2], true) << ";"
				<< "\n\t\t" << get_str_var(term.Jgen[2][0], true)
				<< "  " << get_str_var(term.Jgen[2][1], true)
				<< "  " << get_str_var(term.Jgen[2][2], true)
				<< "\n\t]";
		}

		ofstr	<< ", Bond(" << idx1 << ", " << idx2 << ", [ "
			<< get_str_var(term.dist[0]) << ", "
			<< get_str_var(term.dist[1]) << ", "
			<< get_str_var(term.dist[2])
			<< " ]))\n";
	}


	if(!tl2::equals_0<t_real>(field.mag, g_eps))
	{
		ofstr << "\n# external field\n";
		ofstr << "set_external_field!(magsys, [ "
			<< field.dir[0] << ", "
			<< field.dir[1] << ", "
			<< field.dir[2] << " ] * " << field.mag
			<< ")\n";
	}
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	if(m_dyn.IsIncommensurate())
	{
		ofstr << "\n# supercell for incommensurate structure\n";

		const t_vec_real& prop = m_dyn.GetOrderingWavevector();
		const t_vec_real& axis = m_dyn.GetRotationAxis();
		t_vec_real s0 = tl2::cross(prop, axis);
		s0 /= tl2::norm(s0);

		int sc_x = tl2::equals_0(prop[0], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[0]));
		int sc_y = tl2::equals_0(prop[1], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[1]));
		int sc_z = tl2::equals_0(prop[2], g_eps)
			? 1 : int(std::ceil(t_real(1) / prop[2]));

		ofstr << "magsys = reshape_supercell(magsys, [ "
			<< sc_x << " 0 0; 0 " << sc_y << " 0; 0 0 " << sc_z
			<< " ])\n";

		ofstr << "set_spiral_order!(magsys; "
			<< "k = [ " << prop[0] << ", " << prop[1] << ", " << prop[2] << " ], "
			<< "axis = [ " << axis[0] << ", " << axis[1] << ", " << axis[2] << " ], "
			<< "S0 = [ " << s0[0] << ", " << s0[1] << ", " << s0[2] << " ])\n";
	}

	ofstr << "\n@printf(stderr, \"%s\\n\", magsys)\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << "\n# spin-wave calculation\n";
	ofstr << "@printf(stderr, \"Calculating S(Q, E)...\\n\")\n";

	ofstr << "calc = SpinWaveTheory(magsys; apply_g = true)\n";

	t_real h1 = (t_real)m_Q_start[0]->value();
	t_real k1 = (t_real)m_Q_start[1]->value();
	t_real l1 = (t_real)m_Q_start[2]->value();
	t_real h2 = (t_real)m_Q_end[0]->value();
	t_real k2 = (t_real)m_Q_end[1]->value();
	t_real l2 = (t_real)m_Q_end[2]->value();

	ofstr << "momenta = collect(range([ "
		<< h1 << ", " << k1 << ", " << l1 << " ], [ "
		<< h2 << ", " << k2 << ", " << l2 << " ], "
		<< m_num_points->value() << "))\n";
	std::string proj = m_use_projector->isChecked() ? ":perp" : ":trace";
	ofstr << "energies, correlations = intensities_bands(calc, momenta,\n"
		<< "\tintensity_formula(calc, " << proj
		<< "; kernel = delta_function_kernel))\n";
	// --------------------------------------------------------------------


	// --------------------------------------------------------------------
	ofstr << R"BLOCK(
# output the dispersion and spin-spin correlation
@printf(stderr, "Outputting data, save to \"disp.dat\" and plot with: \n\t%s\n",
	"gnuplot -p -e \"plot \\\"disp.dat\\\" u 1:4:(\\\$5/4) w p pt 7 ps var\"")

@printf(stdout, "# %8s %10s %10s %10s %10s\n",
	"h (rlu)", "k (rlu)", "l (rlu)", "E (meV)", "S(Q, E)")
for q_idx in 1:length(momenta)
	for e_idx in 1:length(energies[q_idx, :])
		@printf(stdout, "%10.4f %10.4f %10.4f %10.4f %10.4f\n",
			momenta[q_idx][1], momenta[q_idx][2], momenta[q_idx][3],
			energies[q_idx, e_idx],
			correlations[q_idx, e_idx])
	end
end
)BLOCK";
	// --------------------------------------------------------------------

	return true;
}
