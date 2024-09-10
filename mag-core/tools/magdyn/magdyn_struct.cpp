/**
 * magnetic dynamics -- calculations for sites and coupling terms
 * @author Tobias Weber <tweber@ill.fr>
 * @date 2022 - 2024
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

#include "magdyn.h"
#include <QtWidgets/QMessageBox>

#include <iostream>
#include <boost/scope_exit.hpp>



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
		const auto& symops = GetSymOpsForCurrentSG();

		SyncToKernel();
		m_dyn.SymmetriseMagneticSites(symops);
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
		const auto& symops = GetSymOpsForCurrentSG();

		SyncToKernel();
		m_dyn.SymmetriseExchangeTerms(symops);
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
		const auto& symops = GetSymOpsForCurrentSG();

		t_real dist_max = m_maxdist->value();
		t_size sc_max = m_maxSC->value();
		t_size couplings_max = m_maxcouplings->value();

		SyncToKernel();
		m_dyn.GeneratePossibleExchangeTerms(dist_max, sc_max, couplings_max);
		m_dyn.CalcSymmetryIndices(symops);
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
 * get the symmetry operators for the currently selected space group
 */
const std::vector<t_mat_real>& MagDynDlg::GetSymOpsForCurrentSG(bool show_err) const
{
	// current space group index
	int sgidx = m_comboSG->itemData(m_comboSG->currentIndex()).toInt();

	if(sgidx < 0 || t_size(sgidx) >= m_SGops.size())
	{
		if(show_err)
		{
			QMessageBox::critical(const_cast<MagDynDlg*>(this),
				"Magnetic Dynamics",
				"Invalid space group selected.");
		}

		// return empty symop list
		static const std::vector<t_mat_real> nullvec{};
		return nullvec;
	}

	return m_SGops[sgidx];
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
void MagDynDlg::ImportAtoms(const std::vector<TableImportAtom>& atompos_vec,
	bool clear_existing)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// remove existing sites
	if(clear_existing)
		DelTabItem(m_sitestab, -1);

	for(const TableImportAtom& atompos : atompos_vec)
	{
		std::string pos_x = "0", pos_y = "0", pos_z = "0";
		std::string spin_x = "0", spin_y = "0", spin_z = "1";
		std::string spin_mag = "1";

		std::string name = "";
		t_size sym_idx = 0;  // TODO

		if(atompos.name) name = *atompos.name;
		if(atompos.x) pos_x = tl2::var_to_str(*atompos.x);
		if(atompos.y) pos_y = tl2::var_to_str(*atompos.y);
		if(atompos.z) pos_z = tl2::var_to_str(*atompos.z);
		if(atompos.Sx) spin_x = tl2::var_to_str(*atompos.Sx);
		if(atompos.Sy) spin_y = tl2::var_to_str(*atompos.Sy);
		if(atompos.Sz) spin_z = tl2::var_to_str(*atompos.Sz);
		if(atompos.Smag) spin_mag = tl2::var_to_str(*atompos.Smag);

		AddSiteTabItem(-1, name, sym_idx,
			pos_x, pos_y, pos_z,
			spin_x, spin_y, spin_z, spin_mag);
	}
}



/**
 * import magnetic couplings from table dialog
 */
void MagDynDlg::ImportCouplings(const std::vector<TableImportCoupling>& couplings,
	bool clear_existing)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_ignoreCalc = false;
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END
	m_ignoreCalc = true;

	// remove existing couplings
	if(clear_existing)
		DelTabItem(m_termstab, -1);

	for(const TableImportCoupling& coupling : couplings)
	{
		t_size atom_1 = 0, atom_2 = 0;
		t_size sym_idx = 0;  // TODO

		std::string dist_x = "0", dist_y = "0", dist_z = "0";
		std::string J = "0";
		std::string dmi_x = "0", dmi_y = "0", dmi_z = "0";

		std::string name = "";
		if(coupling.name) name = *coupling.name;
		if(coupling.atomidx1) atom_1 = *coupling.atomidx1;
		if(coupling.atomidx2) atom_2 = *coupling.atomidx2;
		if(coupling.dx) dist_x = tl2::var_to_str(*coupling.dx);
		if(coupling.dy) dist_y = tl2::var_to_str(*coupling.dy);
		if(coupling.dz) dist_z = tl2::var_to_str(*coupling.dz);
		if(coupling.J) J = tl2::var_to_str(*coupling.J);
		if(coupling.dmix) dmi_x = tl2::var_to_str(*coupling.dmix);
		if(coupling.dmiy) dmi_y = tl2::var_to_str(*coupling.dmiy);
		if(coupling.dmiz) dmi_z = tl2::var_to_str(*coupling.dmiz);

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

		AddTermTabItem(-1, name, sym_idx, atom_1_name, atom_2_name,
			dist_x, dist_y, dist_z, J, dmi_x, dmi_y, dmi_z);
	}
}



/**
 * assign symmetry groups to sites and couplings
 */
void MagDynDlg::CalcSymmetryIndices()
{
	const auto& symops = GetSymOpsForCurrentSG();
	m_dyn.CalcSymmetryIndices(symops);
	SyncSymmetryIndicesFromKernel();
}
