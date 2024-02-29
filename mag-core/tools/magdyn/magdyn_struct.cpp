/**
 * magnetic dynamics -- calculations for sites and coupling terms
 * @author Tobias Weber <tweber@ill.fr>
 * @date Jan-2022
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
#include <QtWidgets/QMessageBox>

#include <iostream>
#include <boost/scope_exit.hpp>

namespace pt = boost::property_tree;



/**
 * flip the coordinates of the magnetic site positions
 * (e.g. to get the negative phase factor for weights)
 */
void MagDynDlg::MirrorAtoms()
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// iterate the magnetic sites
	for(int row=0; row<m_sitestab->rowCount(); ++row)
	{
		auto *pos_x = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_X));
		auto *pos_y = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_Y));
		auto *pos_z = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_sitestab->item(row, COL_SITE_POS_Z));

		if(!pos_x || !pos_y || !pos_z)
		{
			std::cerr << "Invalid entry in sites table row "
				<< row << "." << std::endl;
			continue;
		}

		pos_x->SetValue(-pos_x->GetValue());
		pos_y->SetValue(-pos_y->GetValue());
		pos_z->SetValue(-pos_z->GetValue());
	}
}



/**
 * rotate the direction of the magnetic field
 */
void MagDynDlg::RotateField(bool ccw)
{
	t_vec_real axis = tl2::create<t_vec_real>(
	{
		(t_real)m_rot_axis[0]->value(),
		(t_real)m_rot_axis[1]->value(),
		(t_real)m_rot_axis[2]->value(),
	});

	t_vec_real B = tl2::create<t_vec_real>(
	{
		(t_real)m_field_dir[0]->value(),
		(t_real)m_field_dir[1]->value(),
		(t_real)m_field_dir[2]->value(),
	});

	t_real angle = m_rot_angle->value() / 180.*tl2::pi<t_real>;
	if(!ccw)
		angle = -angle;

	t_mat_real R = tl2::rotation<t_mat_real, t_vec_real>(
		axis, angle, false);
	B = R*B;
	tl2::set_eps_0(B, g_eps);

	for(int i=0; i<3; ++i)
	{
		m_field_dir[i]->blockSignals(true);
		m_field_dir[i]->setValue(B[i]);
		m_field_dir[i]->blockSignals(false);
	}

	if(m_autocalc->isChecked())
		CalcAll();
};



/**
 * set selected field as current
 */
void MagDynDlg::SetCurrentField()
{
	if(m_fields_cursor_row < 0 || m_fields_cursor_row >= m_fieldstab->rowCount())
		return;

	const auto* Bh = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_H));
	const auto* Bk = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_K));
	const auto* Bl = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_L));
	const auto* Bmag = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
		m_fieldstab->item(m_fields_cursor_row, COL_FIELD_MAG));

	if(!Bh || !Bk || !Bl || !Bmag)
		return;

	m_ignoreCalc = true;
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	m_field_dir[0]->setValue(Bh->GetValue());
	m_field_dir[1]->setValue(Bk->GetValue());
	m_field_dir[2]->setValue(Bl->GetValue());
	m_field_mag->setValue(Bmag->GetValue());
}



/**
 * generate magnetic sites form the space group symmetries
 */
void MagDynDlg::GenerateSitesFromSG()
{
	try
	{
		// symops of current space group
		auto sgidx = m_comboSG->itemData(m_comboSG->currentIndex()).toInt();
		if(sgidx < 0 || t_size(sgidx) >= m_SGops.size())
		{
			QMessageBox::critical(this, "Magnetic Dynamics",
				"Invalid space group selected.");
			return;
		}

		SyncToKernel();
		m_dyn.SymmetriseMagneticSites(m_SGops[sgidx]);
		SyncSitesFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", ex.what());
	}
}



/**
 * generate exchange terms from space group symmetries
 */
void MagDynDlg::GenerateCouplingsFromSG()
{
	try
	{
		// symops of current space group
		auto sgidx = m_comboSG->itemData(m_comboSG->currentIndex()).toInt();
		if(sgidx < 0 || t_size(sgidx) >= m_SGops.size())
		{
			QMessageBox::critical(this, "Magnetic Dynamics",
				"Invalid space group selected.");
			return;
		}

		SyncToKernel();
		m_dyn.SymmetriseExchangeTerms(m_SGops[sgidx]);
		SyncTermsFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", ex.what());
	}
}



/**
 * generate possible couplings up to a certain distance
 */
void MagDynDlg::GeneratePossibleCouplings()
{
	try
	{
		t_real dist_max = m_maxdist->value();
		t_size sc_max = m_maxSC->value();
		t_size couplings_max = m_maxcouplings->value();

		t_real a = m_xtallattice[0]->value();
		t_real b = m_xtallattice[1]->value();
		t_real c = m_xtallattice[2]->value();
		t_real alpha = m_xtalangles[0]->value() / 180. * tl2::pi<t_real>;
		t_real beta = m_xtalangles[1]->value() / 180. * tl2::pi<t_real>;
		t_real gamma = m_xtalangles[2]->value() / 180. * tl2::pi<t_real>;

		SyncToKernel();
		m_dyn.GeneratePossibleExchangeTerms(
			a, b, c, alpha, beta, gamma,
			dist_max, sc_max, couplings_max);
		SyncTermsFromKernel();

		if(m_autocalc->isChecked())
			CalcAll();
	}
	catch(const std::exception& ex)
	{
		QMessageBox::critical(this, "Magnetic Dynamics", ex.what());
	}
}



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

	// get temperature
	if(m_use_temperature->isChecked())
	{
		t_real temp = m_temperature->value();
		m_dyn.SetTemperature(temp);
	}

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
		site.g = -2. * tl2::unit<t_mat>(3);

		site.pos = tl2::create<t_vec_real>(
		{
			pos_x->GetValue(),
			pos_y->GetValue(),
			pos_z->GetValue(),
		});

		site.spin_mag = spin_mag->GetValue();
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
		term.dist = tl2::create<t_vec_real>(
		{
			dist_x->GetValue(),
			dist_y->GetValue(),
			dist_z->GetValue(),
		});

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



/**
 * open the table import dialog
 */
void MagDynDlg::ShowTableImporter()
{
	if(!m_table_import_dlg)
	{
		m_table_import_dlg = new TableImportDlg(this, m_sett);

		connect(m_table_import_dlg, &TableImportDlg::SetAtomsSignal,
			this, &MagDynDlg::ImportAtoms);
		connect(m_table_import_dlg, &TableImportDlg::SetCouplingsSignal,
				this, &MagDynDlg::ImportCouplings);
	}

	m_table_import_dlg->show();
	m_table_import_dlg->raise();
	m_table_import_dlg->activateWindow();
}



/**
 * import magnetic site positions from table dialog
 */
void MagDynDlg::ImportAtoms(const std::vector<TableImportAtom>& atompos_vec)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	DelTabItem(m_sitestab, -1);  // remove original sites

	for(const TableImportAtom& atompos : atompos_vec)
	{
		t_real pos_x = 0, pos_y = 0, pos_z = 0;
		std::string spin_x = "0", spin_y = "0", spin_z = "1";
		t_real spin_mag = 1;

		std::string name = "";
		if(atompos.name) name = *atompos.name;
		if(atompos.x) pos_x = *atompos.x;
		if(atompos.y) pos_y = *atompos.y;
		if(atompos.z) pos_z = *atompos.z;
		if(atompos.Sx) spin_x = tl2::var_to_str(*atompos.Sx);
		if(atompos.Sy) spin_y = tl2::var_to_str(*atompos.Sy);
		if(atompos.Sz) spin_z = tl2::var_to_str(*atompos.Sz);
		if(atompos.Smag) spin_mag = *atompos.Smag;

		if(name == "")
		{
			// default site name
			std::ostringstream ostrName;
			ostrName << "site_" << (++m_curSiteCtr);
			name = ostrName.str();
		}

		AddSiteTabItem(-1, name,
			pos_x, pos_y, pos_z,
			spin_x, spin_y, spin_z, spin_mag);
	}
}



/**
 * import magnetic couplings from table dialog
 */
void MagDynDlg::ImportCouplings(const std::vector<TableImportCoupling>& couplings)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	DelTabItem(m_termstab, -1);  // remove original couplings

	for(const TableImportCoupling& coupling : couplings)
	{
		t_size atom_1 = 0, atom_2 = 0;
		t_real dist_x = 0, dist_y = 0, dist_z = 0;
		std::string J = "0";
		std::string dmi_x = "0", dmi_y = "0", dmi_z = "0";

		std::string name = "";
		if(coupling.name) name = *coupling.name;
		if(coupling.atomidx1) atom_1 = *coupling.atomidx1;
		if(coupling.atomidx2) atom_2 = *coupling.atomidx2;
		if(coupling.dx) dist_x = *coupling.dx;
		if(coupling.dy) dist_y = *coupling.dy;
		if(coupling.dz) dist_z = *coupling.dz;
		if(coupling.J) J = tl2::var_to_str(*coupling.J);
		if(coupling.dmix) dmi_x = tl2::var_to_str(*coupling.dmix);
		if(coupling.dmiy) dmi_y = tl2::var_to_str(*coupling.dmiy);
		if(coupling.dmiz) dmi_z = tl2::var_to_str(*coupling.dmiz);

		if(name == "")
		{
			// default coupling name
			std::ostringstream ostrName;
			ostrName << "coupling_" << (++m_curCouplingCtr);
			name = ostrName.str();
		}

		std::string atom_1_name = tl2::var_to_str(atom_1);
		std::string atom_2_name = tl2::var_to_str(atom_2);

		// get the site names from the table
		if((int)atom_1 < m_sitestab->rowCount())
		{
			if(auto *name = m_sitestab->item(atom_1, COL_SITE_NAME); name)
				atom_1_name = name->text().toStdString();
		}
		if((int)atom_2 < m_sitestab->rowCount())
		{
			if(auto *name = m_sitestab->item(atom_2, COL_SITE_NAME); name)
				atom_2_name = name->text().toStdString();
		}

		AddTermTabItem(-1, name, atom_1_name, atom_2_name,
			dist_x, dist_y, dist_z, J, dmi_x, dmi_y, dmi_z);
	}
}
