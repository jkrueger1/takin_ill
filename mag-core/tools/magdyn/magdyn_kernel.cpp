/**
 * magnetic dynamics -- synchronisation with magdyn kernel
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
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

#include "magdyn.h"
#include "tlibs2/libs/units.h"

#include <iostream>
#include <boost/scope_exit.hpp>

namespace pt = boost::property_tree;



/**
 * get the magnetic sites from the kernel
 * and add them to the table
 */
void MagDynDlg::SyncSitesFromKernel(boost::optional<const pt::ptree&> extra_infos)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		this_->m_ignoreSitesCalc = false;
	} BOOST_SCOPE_EXIT_END

	// prevent syncing before the new sites are transferred
	m_ignoreCalc = true;
	m_ignoreSitesCalc = true;

	// clear old sites
	DelTabItem(m_sitestab, -1);

	for(t_size site_index = 0; site_index < m_dyn.GetMagneticSitesCount(); ++site_index)
	{
		const auto &site = m_dyn.GetMagneticSite(site_index);

		// default colour
		std::string rgb = "auto";

		// get additional data from exchange term entry
		if(extra_infos && site_index < extra_infos->size())
		{
			auto siteiter = (*extra_infos).begin();
			std::advance(siteiter, site_index);

			// read colour
			rgb = siteiter->second.get<std::string>("colour", "auto");
		}

		std::string spin_ortho_x = site.spin_ortho[0];
		std::string spin_ortho_y = site.spin_ortho[1];
		std::string spin_ortho_z = site.spin_ortho[2];

		if(spin_ortho_x == "")
			spin_ortho_x = "auto";
		if(spin_ortho_y == "")
			spin_ortho_y = "auto";
		if(spin_ortho_z == "")
			spin_ortho_z = "auto";

		AddSiteTabItem(-1,
			site.name,
			site.pos[0], site.pos[1], site.pos[2],
			site.spin_dir[0], site.spin_dir[1], site.spin_dir[2], site.spin_mag,
			spin_ortho_x, spin_ortho_y, spin_ortho_z,
			rgb);
	}

	m_ignoreSitesCalc = false;
	SyncSiteComboBoxes();
}



/**
 * get the exchange terms from the kernel
 * and add them to the table
 */
void MagDynDlg::SyncTermsFromKernel(boost::optional<const pt::ptree&> extra_infos)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
	} BOOST_SCOPE_EXIT_END

	// prevent syncing before the new terms are transferred
	m_ignoreCalc = true;

	// clear old terms
	DelTabItem(m_termstab, -1);

	for(t_size term_index = 0; term_index < m_dyn.GetExchangeTermsCount(); ++term_index)
	{
		const auto& term = m_dyn.GetExchangeTerm(term_index);

		// default colour
		std::string rgb = "#0x00bf00";

		// get additional data from exchange term entry
		if(extra_infos && term_index < extra_infos->size())
		{
			auto termiter = (*extra_infos).begin();
			std::advance(termiter, term_index);

			// read colour
			rgb = termiter->second.get<std::string>("colour", "#0x00bf00");
		}

		AddTermTabItem(-1,
			term.name, term.site1, term.site2,
			term.dist[0], term.dist[1], term.dist[2],
			term.J,
			term.dmi[0], term.dmi[1], term.dmi[2],
			term.Jgen[0][0], term.Jgen[0][1], term.Jgen[0][2],
			term.Jgen[1][0], term.Jgen[1][1], term.Jgen[1][2],
			term.Jgen[2][0], term.Jgen[2][1], term.Jgen[2][2],
			rgb);
	}
}



/**
 * get the sites, exchange terms, and variables from the table
 * and transfer them to the dynamics calculator
 */
void MagDynDlg::SyncToKernel()
{
	if(m_ignoreCalc)
		return;
	m_dyn.Clear();

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_sitestab->blockSignals(false);
		this_->m_termstab->blockSignals(false);
		this_->m_varstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_sitestab->blockSignals(true);
	m_termstab->blockSignals(true);
	m_varstab->blockSignals(true);

	// get variables
	for(int row=0; row<m_varstab->rowCount(); ++row)
	{
		auto *name = m_varstab->item(row, COL_VARS_NAME);
		auto *val_re = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_varstab->item(row, COL_VARS_VALUE_REAL));
		auto *val_im = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_varstab->item(row, COL_VARS_VALUE_IMAG));

		if(!name || !val_re || !val_im)
		{
			std::cerr << "Invalid entry in variables table row "
				<< row << "." << std::endl;
			continue;
		}

		t_magdyn::Variable var;
		var.name = name->text().toStdString();
		var.value = val_re->GetValue() + val_im->GetValue() * t_cplx(0, 1);

		m_dyn.AddVariable(std::move(var));
	}

	// get ordering vector and rotation axis
	{
		t_vec_real ordering = tl2::create<t_vec_real>(
		{
			(t_real)m_ordering[0]->value(),
			(t_real)m_ordering[1]->value(),
			(t_real)m_ordering[2]->value(),
		});

		t_vec_real rotaxis = tl2::create<t_vec_real>(
		{
			(t_real)m_normaxis[0]->value(),
			(t_real)m_normaxis[1]->value(),
			(t_real)m_normaxis[2]->value(),
		});

		m_dyn.SetOrderingWavevector(ordering);
		m_dyn.SetRotationAxis(rotaxis);
	}

	// get external field
	if(m_use_field->isChecked())
	{
		t_magdyn::ExternalField field;
		field.dir = tl2::create<t_vec_real>(
		{
			(t_real)m_field_dir[0]->value(),
			(t_real)m_field_dir[1]->value(),
			(t_real)m_field_dir[2]->value(),
		});

		field.mag = m_field_mag->value();
		field.align_spins = m_align_spins->isChecked();

		m_dyn.SetExternalField(field);
	}

	m_dyn.CalcExternalField();

	// get temperature
	if(m_use_temperature->isChecked())
	{
		t_real temp = m_temperature->value();
		m_dyn.SetTemperature(temp);
	}

	// get magnetic sites
	for(int row=0; row<m_sitestab->rowCount(); ++row)
	{
		auto *name = m_sitestab->item(row, COL_SITE_NAME);
		auto *pos_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_X));
		auto *pos_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_Y));
		auto *pos_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_Z));
		auto *spin_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_SPIN_X));
		auto *spin_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_SPIN_Y));
		auto *spin_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_SPIN_Z));
		auto *spin_mag = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_SPIN_MAG));

		tl2::NumericTableWidgetItem<t_real> *spin_ortho_x = nullptr;
		tl2::NumericTableWidgetItem<t_real> *spin_ortho_y = nullptr;
		tl2::NumericTableWidgetItem<t_real> *spin_ortho_z = nullptr;
		if(m_allow_ortho_spin)
		{
			spin_ortho_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_sitestab->item(row, COL_SITE_SPIN_ORTHO_X));
			spin_ortho_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_sitestab->item(row, COL_SITE_SPIN_ORTHO_Y));
			spin_ortho_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_sitestab->item(row, COL_SITE_SPIN_ORTHO_Z));
		}

		if(!name || !pos_x || !pos_y || !pos_z ||
			!spin_x || !spin_y || !spin_z || !spin_mag)
		{
			std::cerr << "Invalid entry in sites table row "
				<< row << "." << std::endl;
			continue;
		}

		t_magdyn::MagneticSite site;
		site.name = name->text().toStdString();
		site.g_e = tl2::g_e<t_real> * tl2::unit<t_mat>(3);

		site.pos[0] = pos_x->text().toStdString();
		site.pos[1] = pos_y->text().toStdString();
		site.pos[2] = pos_z->text().toStdString();

		site.spin_mag = spin_mag->text().toStdString();
		site.spin_dir[0] = spin_x->text().toStdString();
		site.spin_dir[1] = spin_y->text().toStdString();
		site.spin_dir[2] = spin_z->text().toStdString();

		site.spin_ortho[0] = "";
		site.spin_ortho[1] = "";
		site.spin_ortho[2] = "";

		if(m_allow_ortho_spin)
		{
			if(spin_ortho_x && spin_ortho_x->text() != "" && spin_ortho_x->text() != "auto")
				site.spin_ortho[0] = spin_ortho_x->text().toStdString();
			if(spin_ortho_y && spin_ortho_y->text() != "" && spin_ortho_y->text() != "auto")
				site.spin_ortho[1] = spin_ortho_y->text().toStdString();
			if(spin_ortho_z && spin_ortho_z->text() != "" && spin_ortho_z->text() != "auto")
				site.spin_ortho[2] = spin_ortho_z->text().toStdString();
		}

		m_dyn.AddMagneticSite(std::move(site));
	}

	m_dyn.CalcMagneticSites();

	// get exchange terms
	for(int row=0; row<m_termstab->rowCount(); ++row)
	{
		auto *name = m_termstab->item(row, COL_XCH_NAME);
		auto *dist_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DIST_X));
		auto *dist_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DIST_Y));
		auto *dist_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DIST_Z));
		auto *interaction = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_INTERACTION));
		auto *dmi_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DMI_X));
		auto *dmi_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DMI_Y));
		auto *dmi_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_termstab->item(row, COL_XCH_DMI_Z));
		QComboBox *site_1 = reinterpret_cast<QComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM1_IDX));
		QComboBox *site_2 = reinterpret_cast<QComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM2_IDX));

		tl2::NumericTableWidgetItem<t_real>* gen_xx = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_xy = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_xz = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_yx = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_yy = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_yz = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_zx = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_zy = nullptr;
		tl2::NumericTableWidgetItem<t_real>* gen_zz = nullptr;
		if(m_allow_general_J)
		{
			gen_xx = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_XX));
			gen_xy = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_XY));
			gen_xz = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_XZ));
			gen_yx = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_YX));
			gen_yy = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_YY));
			gen_yz = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_YZ));
			gen_zx = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_ZX));
			gen_zy = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_ZY));
			gen_zz = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
				m_termstab->item(row, COL_XCH_GEN_ZZ));
		}

		if(!name || !site_1 || !site_2 ||
			!dist_x || !dist_y || !dist_z ||
			!interaction || !dmi_x || !dmi_y || !dmi_z)
		{
			std::cerr << "Invalid entry in couplings table row "
				<< row << "." << std::endl;
			continue;
		}

		t_magdyn::ExchangeTerm term;
		term.name = name->text().toStdString();

		term.site1 = site_1->currentText().toStdString();
		term.site2 = site_2->currentText().toStdString();

		term.dist[0] = dist_x->text().toStdString();
		term.dist[1] = dist_y->text().toStdString();
		term.dist[2] = dist_z->text().toStdString();

		term.J = interaction->text().toStdString();

		if(m_use_dmi->isChecked())
		{
			term.dmi[0] = dmi_x->text().toStdString();
			term.dmi[1] = dmi_y->text().toStdString();
			term.dmi[2] = dmi_z->text().toStdString();
		}

		if(m_allow_general_J && m_use_genJ->isChecked())
		{
			term.Jgen[0][0] = gen_xx->text().toStdString();
			term.Jgen[0][1] = gen_xy->text().toStdString();
			term.Jgen[0][2] = gen_xz->text().toStdString();
			term.Jgen[1][0] = gen_yx->text().toStdString();
			term.Jgen[1][1] = gen_yy->text().toStdString();
			term.Jgen[1][2] = gen_yz->text().toStdString();
			term.Jgen[2][0] = gen_zx->text().toStdString();
			term.Jgen[2][1] = gen_zy->text().toStdString();
			term.Jgen[2][2] = gen_zz->text().toStdString();
		}

		m_dyn.AddExchangeTerm(std::move(term));
	}

	m_dyn.CalcExchangeTerms();

	// ground state energy
	std::ostringstream ostrGS;
	ostrGS.precision(g_prec_gui);
	ostrGS << "E0 = " << m_dyn.CalcGroundStateEnergy() << " meV";
	m_statusFixed->setText(ostrGS.str().c_str());
}
