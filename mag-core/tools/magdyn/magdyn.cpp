/**
 * magnetic dynamics -- main dialog handler functions
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
#include "libs/loadcif.h"

#include <QtCore/QMimeData>
#include <iostream>



MagDynDlg::MagDynDlg(QWidget* pParent) : QDialog{pParent},
	m_sett{new QSettings{"takin", "magdyn", this}}
{
	// restore settings done from takin main settings dialog
	get_settings_from_takin_core();
	if(g_font != "")
	{
		QFont font = this->font();
		if(font.fromString(g_font))
			setFont(font);
	}

	InitSettingsDlg();

	// read settings that require a restart
	m_allow_ortho_spin = (g_allow_ortho_spin != 0);
	m_allow_general_J = (g_allow_general_J != 0);

	// create gui
	CreateMainWindow();
	CreateMenuBar();

	// create dialogs
	CreateInfoDlg();
	CreateNotesDlg();

	// create input panels
	CreateSitesPanel();
	CreateExchangeTermsPanel();
	CreateSamplePanel();
	CreateSampleEnvPanel();
	CreateVariablesPanel();

	// create output panels
	CreateDispersionPanel();
	CreateHamiltonPanel();
	CreateCoordinatesPanel();
	CreateExportPanel();


	// get space groups and symops
	auto spacegroups = get_sgs<t_mat_real>();
	m_SGops.clear();
	m_SGops.reserve(spacegroups.size());

	for(QComboBox* combo : {m_comboSG, m_comboSGSites, m_comboSGTerms})
	{
		if(!combo)
			continue;
		combo->clear();
	}

	for(auto [sgnum, descr, ops] : spacegroups)
	{
		for(QComboBox* combo : {m_comboSG, m_comboSGSites, m_comboSGTerms})
		{
			if(!combo)
				continue;
			combo->addItem(descr.c_str(), combo->count());
		}

		m_SGops.emplace_back(std::move(ops));
	}


	InitSettings();

	// restore settings
	if(m_sett)
	{
		// restore window size and position
		if(m_sett->contains("geo"))
			restoreGeometry(m_sett->value("geo").toByteArray());
		else
			resize(800, 600);

		if(m_sett->contains("recent_files"))
			m_recent.SetRecentFiles(m_sett->value("recent_files").toStringList());

		if(m_sett->contains("recent_struct_files"))
			m_recent_struct.SetRecentFiles(m_sett->value("recent_struct_files").toStringList());

		if(m_sett->contains("splitter"))
			m_split_inout->restoreState(m_sett->value("splitter").toByteArray());
	}

	setAcceptDrops(true);
}



MagDynDlg::~MagDynDlg()
{
	Clear();

	if(m_settings_dlg)
	{
		delete m_settings_dlg;
		m_settings_dlg = nullptr;
	}

	if(m_structplot_dlg)
	{
		delete m_structplot_dlg;
		m_structplot_dlg = nullptr;
	}

	if(m_table_import_dlg)
	{
		delete m_table_import_dlg;
		m_table_import_dlg = nullptr;
	}

	if(m_notes_dlg)
	{
		delete m_notes_dlg;
		m_notes_dlg = nullptr;
	}

	if(m_info_dlg)
	{
		delete m_info_dlg;
		m_info_dlg = nullptr;
	}
}



void MagDynDlg::mousePressEvent(QMouseEvent *evt)
{
	QDialog::mousePressEvent(evt);
}



/**
 * dialog is closing
 */
void MagDynDlg::closeEvent(QCloseEvent *)
{
	if(!m_sett)
		return;

	m_recent.TrimEntries();
	m_sett->setValue("recent_files", m_recent.GetRecentFiles());

	m_recent_struct.TrimEntries();
	m_sett->setValue("recent_struct_files", m_recent_struct.GetRecentFiles());

	m_sett->setValue("geo", saveGeometry());

	if(m_split_inout)
		m_sett->setValue("splitter", m_split_inout->saveState());

	if(m_structplot_dlg)
		m_sett->setValue("geo_struct_view", m_structplot_dlg->saveGeometry());
}



/**
 * a file is being dragged over the window
 */
void MagDynDlg::dragEnterEvent(QDragEnterEvent *evt)
{
	if(evt)
		evt->accept();
}



/**
 * a file is being dropped onto the window
 */
void MagDynDlg::dropEvent(QDropEvent *evt)
{
	const QMimeData *mime = evt->mimeData();
	if(!mime)
		return;

	for(const QUrl& url : mime->urls())
	{
		if(!url.isLocalFile())
			continue;

		Load(url.toLocalFile());
		evt->accept();
		break;
	}
}



/**
 * refresh and calculate everything
 */
void MagDynDlg::CalcAll()
{
	// calculate structure
	SyncToKernel();
	StructPlotSync();

	// calculate dynamics
	CalcDispersion();
	CalcHamiltonian();
}



/**
 * enable GUI inputs after calculation threads have finished
 */
void MagDynDlg::EnableInput()
{
	m_tabs_in->setEnabled(true);
	m_tabs_out->setEnabled(true);
	m_menu->setEnabled(true);
	m_btnStart->setEnabled(true);
}



/**
 * disable GUI inputs for calculation threads
 */
void MagDynDlg::DisableInput()
{
	m_menu->setEnabled(false);
	m_tabs_out->setEnabled(false);
	m_tabs_in->setEnabled(false);
	m_btnStart->setEnabled(false);
}
