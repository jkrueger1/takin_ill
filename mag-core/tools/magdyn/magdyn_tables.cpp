/**
 * magnetic dynamics -- table management
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

#include <iostream>
#include <unordered_set>
#include <boost/scope_exit.hpp>



/**
 * set a unique name for the given table item
 */
static void set_unique_tab_item_name(
	QTableWidget *tab, QTableWidgetItem *item, int name_col, const std::string& prefix)
{
	if(!item)
		return;

	// was the item renamed?
	if(tab->column(item) != name_col)
		return;

	int row = tab->row(item);
	if(row < 0 || row >= tab->rowCount())
		return;

	// check if the new name is unique and non-empty
	std::string new_name_base = item->text().toStdString();
	std::string new_name = new_name_base;
	if(tl2::trimmed(new_name) == "")
		new_name_base = new_name = prefix;

	// hash all other names
	std::unordered_set<std::string> used_names;
	for(int idx = 0; idx < tab->rowCount(); ++idx)
	{
		auto *name = tab->item(idx, name_col);
		if(!name || idx == row)
			continue;

		used_names.emplace(name->text().toStdString());
	}

	// possibly add a suffix to make the new name unique
	std::size_t ctr = 1;
	while(used_names.find(new_name) != used_names.end())
		new_name = new_name_base + "_" + tl2::var_to_str(ctr++);

	// rename the site if needed
	if(new_name != item->text().toStdString())
		item->setText(new_name.c_str());
}



/**
 * add an atom site
 */
void MagDynDlg::AddSiteTabItem(
	int row, const std::string& name,
	const std::string& x, const std::string& y, const std::string& z,
	const std::string& sx, const std::string& sy, const std::string& sz,
	const std::string& S,
	const std::string& sox, const std::string& soy, const std::string& soz,
	const std::string& rgb)
{
	bool bclone = false;
	m_sitestab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_sitestab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_sitestab->rowCount();
	else if(row == -2 && m_sites_cursor_row >= 0)	// use row from member variable
		row = m_sites_cursor_row;
	else if(row == -3 && m_sites_cursor_row >= 0)	// use row from member variable +1
		row = m_sites_cursor_row + 1;
	else if(row == -4 && m_sites_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_sites_cursor_row + 1;
		bclone = true;
	}

	m_sitestab->setSortingEnabled(false);
	m_sitestab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < NUM_SITE_COLS; ++thecol)
		{
			m_sitestab->setItem(row, thecol,
				m_sitestab->item(m_sites_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_sitestab->setItem(row, COL_SITE_NAME,
			new QTableWidgetItem(name.c_str()));
		m_sitestab->setItem(row, COL_SITE_POS_X,
			new tl2::NumericTableWidgetItem<t_real>(x));
		m_sitestab->setItem(row, COL_SITE_POS_Y,
			new tl2::NumericTableWidgetItem<t_real>(y));
		m_sitestab->setItem(row, COL_SITE_POS_Z,
			new tl2::NumericTableWidgetItem<t_real>(z));
		m_sitestab->setItem(row, COL_SITE_SPIN_X,
			new tl2::NumericTableWidgetItem<t_real>(sx));
		m_sitestab->setItem(row, COL_SITE_SPIN_Y,
			new tl2::NumericTableWidgetItem<t_real>(sy));
		m_sitestab->setItem(row, COL_SITE_SPIN_Z,
			new tl2::NumericTableWidgetItem<t_real>(sz));
		m_sitestab->setItem(row, COL_SITE_SPIN_MAG,
			new tl2::NumericTableWidgetItem<t_real>(S));
		m_sitestab->setItem(row, COL_SITE_RGB, new QTableWidgetItem(rgb.c_str()));
		if(m_allow_ortho_spin)
		{
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_X,
				new tl2::NumericTableWidgetItem<t_real>(sox));
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_Y,
				new tl2::NumericTableWidgetItem<t_real>(soy));
			m_sitestab->setItem(row, COL_SITE_SPIN_ORTHO_Z,
				new tl2::NumericTableWidgetItem<t_real>(soz));
		}
	}

	set_unique_tab_item_name(m_sitestab, m_sitestab->item(row, COL_SITE_NAME),
		COL_SITE_NAME, "site");

	m_sitestab->scrollToItem(m_sitestab->item(row, 0));
	m_sitestab->setCurrentCell(row, 0);
	m_sitestab->setSortingEnabled(true);

	UpdateVerticalHeader(m_sitestab);
	SyncSiteComboBoxes();
}



/**
 * update the contents of all site selection combo boxes to match the sites table
 */
void MagDynDlg::SyncSiteComboBoxes()
{
	if(m_ignoreSitesCalc)
		return;

	// iterate couplings and update their site selection combo boxes
	for(int row=0; row<m_termstab->rowCount(); ++row)
	{
		SitesComboBox *site_1 = reinterpret_cast<SitesComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM1_IDX));
		SitesComboBox *site_2 = reinterpret_cast<SitesComboBox*>(
			m_termstab->cellWidget(row, COL_XCH_ATOM2_IDX));

		SyncSiteComboBox(site_1, site_1->currentText().toStdString());
		SyncSiteComboBox(site_2, site_2->currentText().toStdString());
	}
}



/**
 * update the contents of a site selection combo box to match the sites table
 */
void MagDynDlg::SyncSiteComboBox(SitesComboBox* combo, const std::string& selected_site)
{
	BOOST_SCOPE_EXIT(combo)
	{
		combo->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	combo->blockSignals(true);

	combo->clear();

	int selected_idx = -1;
	int selected_idx_alt = -1;  // alternate selection in case of renamed site

	// iterate sites to get their names
	std::unordered_set<std::string> seen_names;
	for(int row = 0; row < m_sitestab->rowCount(); ++row)
	{
		auto *name = m_sitestab->item(row, COL_SITE_NAME);
		if(!name)
			continue;

		std::string site_name = name->text().toStdString();
		if(seen_names.find(site_name) != seen_names.end())
			continue;
		seen_names.insert(site_name);

		// add the item
		combo->addItem(name->text());

		// check if this item has to be selected
		if(name->text().toStdString() == selected_site)
			selected_idx = row;

		// check if the selection corresponds to the site's previous name
		std::string old_name = name->data(Qt::UserRole).toString().toStdString();
		if(old_name != "" && old_name == selected_site)
			selected_idx_alt = row;
	}

	combo->addItem("<invalid>");

	if(selected_idx >= 0)
	{
		// select the site
		combo->setCurrentIndex(selected_idx);
	}
	else if(selected_idx_alt >= 0)
	{
		// use alternate selection in case of a renamed site
		combo->setCurrentIndex(selected_idx_alt);
	}
	else
	{
		// set selection to invalid
		combo->setCurrentIndex(combo->count() - 1);
	}
}



/**
 * create a combo box with the site names
 */
SitesComboBox* MagDynDlg::CreateSitesComboBox(const std::string& selected_site)
{
	/**
	 * filter out wheel events (used for the combo boxes in the table)
	 */
	struct WheelEventFilter : public QObject
	{
		WheelEventFilter(QObject* comp) : m_comp{comp} {}
		WheelEventFilter(const WheelEventFilter&) = delete;
		WheelEventFilter& operator=(const WheelEventFilter&) = delete;

		bool eventFilter(QObject *comp, QEvent *evt) override
		{
			// ignore wheel events on the given component
			return comp == m_comp && evt->type() == QEvent::Wheel;
		}

		QObject* m_comp = nullptr;
	};


	SitesComboBox* combo = new SitesComboBox();
	combo->setParent(m_termstab);
	combo->installEventFilter(new WheelEventFilter(combo));
	SyncSiteComboBox(combo, selected_site);

	// connections
	connect(combo, static_cast<void (QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
		[this]()
		{
			TermsTableItemChanged(nullptr);
		});

	return combo;
}



/**
 * add an exchange term
 */
void MagDynDlg::AddTermTabItem(
	int row, const std::string& name,
	const std::string& atom_1, const std::string& atom_2,
	const std::string& dist_x, const std::string& dist_y, const std::string& dist_z,
	const std::string& J,
	const std::string& dmi_x, const std::string& dmi_y, const std::string& dmi_z,
	const std::string& gen_xx, const std::string& gen_xy, const std::string& gen_xz,
	const std::string& gen_yx, const std::string& gen_yy, const std::string& gen_yz,
	const std::string& gen_zx, const std::string& gen_zy, const std::string& gen_zz,
	const std::string& rgb)
{
	bool bclone = false;
	m_termstab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_termstab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_termstab->rowCount();
	else if(row == -2 && m_terms_cursor_row >= 0)	// use row from member variable
		row = m_terms_cursor_row;
	else if(row == -3 && m_terms_cursor_row >= 0)	// use row from member variable +1
		row = m_terms_cursor_row + 1;
	else if(row == -4 && m_terms_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_terms_cursor_row + 1;
		bclone = true;
	}

	m_termstab->setSortingEnabled(false);
	m_termstab->insertRow(row);

	if(bclone)
	{
		for(int thecol = 0; thecol < NUM_XCH_COLS; ++thecol)
		{
			m_termstab->setItem(row, thecol,
				m_termstab->item(m_terms_cursor_row, thecol)->clone());

			// also clone site selection combo boxes
			if(thecol == COL_XCH_ATOM1_IDX || thecol == COL_XCH_ATOM2_IDX)
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					m_termstab->cellWidget(m_terms_cursor_row, thecol));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				m_termstab->setCellWidget(row, thecol, combo);
				m_termstab->setItem(row, thecol, combo);
			}
		}
	}
	else
	{
		m_termstab->setItem(row, COL_XCH_NAME,
			new QTableWidgetItem(name.c_str()));

		SitesComboBox* combo1 = CreateSitesComboBox(atom_1);
		SitesComboBox* combo2 = CreateSitesComboBox(atom_2);
		m_termstab->setCellWidget(row, COL_XCH_ATOM1_IDX, combo1);
		m_termstab->setCellWidget(row, COL_XCH_ATOM2_IDX, combo2);
		m_termstab->setItem(row, COL_XCH_ATOM1_IDX, combo1);
		m_termstab->setItem(row, COL_XCH_ATOM2_IDX, combo2);

		m_termstab->setItem(row, COL_XCH_DIST_X,
			new tl2::NumericTableWidgetItem<t_real>(dist_x));
		m_termstab->setItem(row, COL_XCH_DIST_Y,
			new tl2::NumericTableWidgetItem<t_real>(dist_y));
		m_termstab->setItem(row, COL_XCH_DIST_Z,
			new tl2::NumericTableWidgetItem<t_real>(dist_z));
		m_termstab->setItem(row, COL_XCH_INTERACTION,
			new tl2::NumericTableWidgetItem<t_real>(J));
		m_termstab->setItem(row, COL_XCH_DMI_X,
			new tl2::NumericTableWidgetItem<t_real>(dmi_x));
		m_termstab->setItem(row, COL_XCH_DMI_Y,
			new tl2::NumericTableWidgetItem<t_real>(dmi_y));
		m_termstab->setItem(row, COL_XCH_DMI_Z,
			new tl2::NumericTableWidgetItem<t_real>(dmi_z));
		m_termstab->setItem(row, COL_XCH_RGB, new QTableWidgetItem(rgb.c_str()));
		if(m_allow_general_J)
		{
			m_termstab->setItem(row, COL_XCH_GEN_XX,
				new tl2::NumericTableWidgetItem<t_real>(gen_xx));
			m_termstab->setItem(row, COL_XCH_GEN_XY,
				new tl2::NumericTableWidgetItem<t_real>(gen_xy));
			m_termstab->setItem(row, COL_XCH_GEN_XZ,
				new tl2::NumericTableWidgetItem<t_real>(gen_xz));
			m_termstab->setItem(row, COL_XCH_GEN_YX,
				new tl2::NumericTableWidgetItem<t_real>(gen_yx));
			m_termstab->setItem(row, COL_XCH_GEN_YY,
				new tl2::NumericTableWidgetItem<t_real>(gen_yy));
			m_termstab->setItem(row, COL_XCH_GEN_YZ,
				new tl2::NumericTableWidgetItem<t_real>(gen_yz));
			m_termstab->setItem(row, COL_XCH_GEN_ZX,
				new tl2::NumericTableWidgetItem<t_real>(gen_zx));
			m_termstab->setItem(row, COL_XCH_GEN_ZY,
				new tl2::NumericTableWidgetItem<t_real>(gen_zy));
			m_termstab->setItem(row, COL_XCH_GEN_ZZ,
				new tl2::NumericTableWidgetItem<t_real>(gen_zz));
		}
	}

	set_unique_tab_item_name(m_termstab, m_termstab->item(row, COL_XCH_NAME),
		COL_XCH_NAME, "coupling");

	m_termstab->scrollToItem(m_termstab->item(row, 0));
	m_termstab->setCurrentCell(row, 0);
	m_termstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_termstab);
}



/**
 * add a variable
 */
void MagDynDlg::AddVariableTabItem(int row, const std::string& name, const t_cplx& value)
{
	bool bclone = false;
	m_varstab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_varstab->blockSignals(false);
		if(this_->m_autocalc->isChecked())
			this_->CalcAll();
	} BOOST_SCOPE_EXIT_END

	if(row == -1)	// append to end of table
		row = m_varstab->rowCount();
	else if(row == -2 && m_variables_cursor_row >= 0)	// use row from member variable
		row = m_variables_cursor_row;
	else if(row == -3 && m_variables_cursor_row >= 0)	// use row from member variable +1
		row = m_variables_cursor_row + 1;
	else if(row == -4 && m_variables_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_variables_cursor_row + 1;
		bclone = true;
	}

	m_varstab->setSortingEnabled(false);
	m_varstab->insertRow(row);

	if(bclone)
	{
		for(int thecol=0; thecol<NUM_VARS_COLS; ++thecol)
		{
			m_varstab->setItem(row, thecol,
				m_varstab->item(m_variables_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_varstab->setItem(row, COL_VARS_NAME,
			new QTableWidgetItem(name.c_str()));
		m_varstab->setItem(row, COL_VARS_VALUE_REAL,
			new tl2::NumericTableWidgetItem<t_real>(value.real()));
		m_varstab->setItem(row, COL_VARS_VALUE_IMAG,
			new tl2::NumericTableWidgetItem<t_real>(value.imag()));
	}

	set_unique_tab_item_name(m_varstab, m_varstab->item(row, COL_VARS_NAME),
		COL_VARS_NAME, "var");

	m_varstab->scrollToItem(m_varstab->item(row, 0));
	m_varstab->setCurrentCell(row, 0);
	m_varstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_varstab);
}



/**
 * add a magnetic field
 */
void MagDynDlg::AddFieldTabItem(
	int row,
	t_real Bh, t_real Bk, t_real Bl,
	t_real Bmag)
{
	bool bclone = false;

	if(row == -1)	// append to end of table
		row = m_fieldstab->rowCount();
	else if(row == -2 && m_fields_cursor_row >= 0)	// use row from member variable
		row = m_fields_cursor_row;
	else if(row == -3 && m_fields_cursor_row >= 0)	// use row from member variable +1
		row = m_fields_cursor_row + 1;
	else if(row == -4 && m_fields_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_fields_cursor_row + 1;
		bclone = true;
	}

	m_fieldstab->setSortingEnabled(false);
	m_fieldstab->insertRow(row);

	if(bclone)
	{
		for(int thecol=0; thecol<NUM_FIELD_COLS; ++thecol)
		{
			m_fieldstab->setItem(row, thecol,
				m_fieldstab->item(m_fields_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_fieldstab->setItem(row, COL_FIELD_H,
			new tl2::NumericTableWidgetItem<t_real>(Bh));
		m_fieldstab->setItem(row, COL_FIELD_K,
			new tl2::NumericTableWidgetItem<t_real>(Bk));
		m_fieldstab->setItem(row, COL_FIELD_L,
			new tl2::NumericTableWidgetItem<t_real>(Bl));
		m_fieldstab->setItem(row, COL_FIELD_MAG,
			new tl2::NumericTableWidgetItem<t_real>(Bmag));
	}

	m_fieldstab->scrollToItem(m_fieldstab->item(row, 0));
	m_fieldstab->setCurrentCell(row, 0);
	m_fieldstab->setSortingEnabled(true);

	UpdateVerticalHeader(m_fieldstab);
}



/**
 * add a coordinate path
 */
void MagDynDlg::AddCoordinateTabItem(
	int row,
	t_real hi, t_real ki, t_real li,
	t_real hf, t_real kf, t_real lf)
{
	bool bclone = false;

	if(row == -1)	// append to end of table
		row = m_coordinatestab->rowCount();
	else if(row == -2 && m_coordinates_cursor_row >= 0)	// use row from member variable
		row = m_coordinates_cursor_row;
	else if(row == -3 && m_coordinates_cursor_row >= 0)	// use row from member variable +1
		row = m_coordinates_cursor_row + 1;
	else if(row == -4 && m_coordinates_cursor_row >= 0)	// use row from member variable +1
	{
		row = m_coordinates_cursor_row + 1;
		bclone = true;
	}

	m_coordinatestab->setSortingEnabled(false);
	m_coordinatestab->insertRow(row);

	if(bclone)
	{
		for(int thecol=0; thecol<NUM_COORD_COLS; ++thecol)
		{
			m_coordinatestab->setItem(row, thecol,
				m_coordinatestab->item(m_coordinates_cursor_row, thecol)->clone());
		}
	}
	else
	{
		m_coordinatestab->setItem(row, COL_COORD_HI,
			new tl2::NumericTableWidgetItem<t_real>(hi));
		m_coordinatestab->setItem(row, COL_COORD_KI,
			new tl2::NumericTableWidgetItem<t_real>(ki));
		m_coordinatestab->setItem(row, COL_COORD_LI,
			new tl2::NumericTableWidgetItem<t_real>(li));
		m_coordinatestab->setItem(row, COL_COORD_HF,
			new tl2::NumericTableWidgetItem<t_real>(hf));
		m_coordinatestab->setItem(row, COL_COORD_KF,
			new tl2::NumericTableWidgetItem<t_real>(kf));
		m_coordinatestab->setItem(row, COL_COORD_LF,
			new tl2::NumericTableWidgetItem<t_real>(lf));
	}

	m_coordinatestab->scrollToItem(m_coordinatestab->item(row, 0));
	m_coordinatestab->setCurrentCell(row, 0);
	m_coordinatestab->setSortingEnabled(true);

	UpdateVerticalHeader(m_coordinatestab);
}



/**
 * delete table widget items
 */
void MagDynDlg::DelTabItem(QTableWidget *pTab, int begin, int end)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	// if nothing is selected, clear all items
	if(begin == -1 || pTab->selectedItems().count() == 0)
	{
		pTab->clearContents();
		pTab->setRowCount(0);
	}
	else if(begin == -2)	// clear selected
	{
		for(int row : GetSelectedRows(pTab, true))
		{
			pTab->removeRow(row);
		}
	}
	else if(begin >= 0 && end >= 0)		// clear given range
	{
		for(int row=end-1; row>=begin; --row)
		{
			pTab->removeRow(row);
		}
	}

	UpdateVerticalHeader(pTab);
	if(pTab == m_sitestab)
		SyncSiteComboBoxes();
}



void MagDynDlg::MoveTabItemUp(QTableWidget *pTab)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	pTab->setSortingEnabled(false);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	auto selected = GetSelectedRows(pTab, false);
	for(int row : selected)
	{
		if(row == 0)
			continue;

		auto *item = pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		pTab->insertRow(row-1);
		for(int col=0; col<pTab->columnCount(); ++col)
		{
			pTab->setItem(row-1, col, pTab->item(row+1, col)->clone());

			// also clone site selection combo boxes
			if(pTab == m_termstab && (col == COL_XCH_ATOM1_IDX || col == COL_XCH_ATOM2_IDX))
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					pTab->cellWidget(row+1, col));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				pTab->setCellWidget(row-1, col, combo);
				pTab->setItem(row-1, col, combo);
			}
		}
		pTab->removeRow(row+1);
	}

	for(int row=0; row<pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row+1) != selected.end())
		{
			for(int col=0; col<pTab->columnCount(); ++col)
				pTab->item(row, col)->setSelected(true);
		}
	}

	UpdateVerticalHeader(pTab);
}



void MagDynDlg::MoveTabItemDown(QTableWidget *pTab)
{
	bool needs_recalc = true;
	if(pTab == m_fieldstab || pTab == m_coordinatestab)
		needs_recalc = false;

	if(needs_recalc)
		pTab->blockSignals(true);
	pTab->setSortingEnabled(false);
	BOOST_SCOPE_EXIT(this_, pTab, needs_recalc)
	{
		if(needs_recalc)
		{
			pTab->blockSignals(false);
			if(this_->m_autocalc->isChecked())
				this_->CalcAll();
		}
	} BOOST_SCOPE_EXIT_END

	auto selected = GetSelectedRows(pTab, true);
	for(int row : selected)
	{
		if(row == pTab->rowCount()-1)
			continue;

		auto *item = pTab->item(row, 0);
		if(!item || !item->isSelected())
			continue;

		pTab->insertRow(row+2);
		for(int col=0; col<pTab->columnCount(); ++col)
		{
			pTab->setItem(row+2, col, pTab->item(row, col)->clone());

			// also clone site selection combo boxes
			if(pTab == m_termstab && (col == COL_XCH_ATOM1_IDX || col == COL_XCH_ATOM2_IDX))
			{
				SitesComboBox *combo_old = reinterpret_cast<SitesComboBox*>(
					pTab->cellWidget(row, col));
				SitesComboBox *combo = CreateSitesComboBox(
					combo_old->currentText().toStdString());
				pTab->setCellWidget(row+2, col, combo);
				pTab->setItem(row+2, col, combo);
			}
		}
		pTab->removeRow(row);
	}

	for(int row=0; row<pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0);
			item && std::find(selected.begin(), selected.end(), row-1) != selected.end())
		{
			for(int col=0; col<pTab->columnCount(); ++col)
				pTab->item(row, col)->setSelected(true);
		}
	}

	UpdateVerticalHeader(pTab);
}



/**
 * insert a vertical header column showing the row index
 */
void MagDynDlg::UpdateVerticalHeader(QTableWidget *pTab)
{
	for(int row=0; row<pTab->rowCount(); ++row)
	{
		QTableWidgetItem *item = pTab->verticalHeaderItem(row);
		if(!item)
			item = new QTableWidgetItem{};
		item->setText(QString::number(row));
		pTab->setVerticalHeaderItem(row, item);
	}
}



std::vector<int> MagDynDlg::GetSelectedRows(QTableWidget *pTab, bool sort_reversed) const
{
	std::vector<int> vec;
	vec.reserve(pTab->selectedItems().size());

	for(int row=0; row<pTab->rowCount(); ++row)
	{
		if(auto *item = pTab->item(row, 0); item && item->isSelected())
			vec.push_back(row);
	}

	if(sort_reversed)
	{
		std::stable_sort(vec.begin(), vec.end(), [](int row1, int row2)
		{ return row1 > row2; });
	}

	return vec;
}



/**
 * sites table item contents changed
 */
void MagDynDlg::SitesTableItemChanged(QTableWidgetItem *item)
{
	if(!item)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_sitestab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_sitestab->blockSignals(true);

	// was the site renamed?
	if(m_sitestab->column(item) == COL_SITE_NAME)
	{
		int row = m_sitestab->row(item);
		if(row >= 0 && row < m_sitestab->rowCount() && row < (int)m_dyn.GetMagneticSitesCount())
		{
			// get the previous name of the site
			std::string old_name = m_dyn.GetMagneticSite(row).name;
			// and set it as temporary, alternate site name
			item->setData(Qt::UserRole, QString(old_name.c_str()));

			set_unique_tab_item_name(m_sitestab, item, COL_SITE_NAME, "site");
		}
	}

	SyncSiteComboBoxes();

	// clear alternate name
	item->setData(Qt::UserRole, QVariant());

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * coupling table item contents changed
 */
void MagDynDlg::TermsTableItemChanged(QTableWidgetItem *item)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_termstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_termstab->blockSignals(true);

	set_unique_tab_item_name(m_termstab, item, COL_XCH_NAME, "coupling");

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * variable table item contents changed
 */
void MagDynDlg::VariablesTableItemChanged(QTableWidgetItem *item)
{
	BOOST_SCOPE_EXIT(this_)
	{
		this_->m_varstab->blockSignals(false);
	} BOOST_SCOPE_EXIT_END
	m_varstab->blockSignals(true);

	set_unique_tab_item_name(m_varstab, item, COL_VARS_NAME, "var");

	if(m_autocalc->isChecked())
		CalcAll();
}



/**
 * show the popup menu for the tables
 */
void MagDynDlg::ShowTableContextMenu(
	QTableWidget *pTab, QMenu *pMenu, QMenu *pMenuNoItem, const QPoint& ptLocal)
{
	QPoint pt = ptLocal;
	// transform the point to widget coordinates first if it has a viewport
	if(pTab->viewport())
		pt = pTab->viewport()->mapToParent(pt);

	// transform the point to global coordinates
	auto ptGlob = pTab->mapToGlobal(pt);

	if(const auto* item = pTab->itemAt(ptLocal); item)
		pMenu->popup(ptGlob);
	else
		pMenuNoItem->popup(ptGlob);
}
