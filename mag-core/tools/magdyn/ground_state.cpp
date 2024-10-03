/**
 * magnetic dynamics -- minimise ground state energy
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

#include "ground_state.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QTableWidgetItem>
#include <QtWidgets/QDialogButtonBox>

#include "tlibs2/libs/qt/numerictablewidgetitem.h"



/**
 * columns of the sites table
 */
enum : int
{
	COL_GSSITE_NAME = 0,                // label
	COL_GSSITE_PHI, COL_GSSITE_THETA,   // spin direction

	NUM_GSSITE_COLS
};



/**
 * shows the 3d view of the magnetic structure
 */
GroundStateDlg::GroundStateDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Ground State Energy Minimisation");
	setSizeGripEnabled(true);

	m_sitestab = new QTableWidget(this);
	m_sitestab->setShowGrid(true);
	m_sitestab->setAlternatingRowColors(true);
	m_sitestab->setSortingEnabled(true);
	m_sitestab->setMouseTracking(true);
	m_sitestab->setSelectionBehavior(QTableWidget::SelectRows);
	m_sitestab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_sitestab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});
	//m_sitestab->setContextMenuPolicy(Qt::CustomContextMenu);

	m_sitestab->verticalHeader()->setDefaultSectionSize(
		fontMetrics().lineSpacing() + 4);
	m_sitestab->verticalHeader()->setVisible(true);

	m_sitestab->setColumnCount(NUM_GSSITE_COLS);

	m_sitestab->setHorizontalHeaderItem(COL_GSSITE_NAME, new QTableWidgetItem{"Name"});
	m_sitestab->setHorizontalHeaderItem(COL_GSSITE_PHI, new QTableWidgetItem{"Spin u(φ)"});
	m_sitestab->setHorizontalHeaderItem(COL_GSSITE_THETA, new QTableWidgetItem{"Spin v(θ)"});

	m_sitestab->setColumnWidth(COL_GSSITE_NAME, 90);
	m_sitestab->setColumnWidth(COL_GSSITE_PHI, 90);
	m_sitestab->setColumnWidth(COL_GSSITE_THETA, 90);

	QPushButton *btnFromKernel = new QPushButton("Get Spins", this);
	QPushButton *btnToKernel = new QPushButton("Set Spins", this);
	QPushButton *btnOk = new QPushButton("OK", this);

	int y = 0;
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);
	grid->addWidget(m_sitestab, y++, 0, 1, 4);
	grid->addWidget(btnFromKernel, y, 0, 1, 1);
	grid->addWidget(btnToKernel, y, 1, 1, 1);
	grid->addWidget(btnOk, y++, 3, 1, 1);

	if(m_sett && m_sett->contains("ground_state/geo"))
		restoreGeometry(m_sett->value("ground_state/geo").toByteArray());
	else
		resize(800, 800);

	connect(btnFromKernel, &QAbstractButton::clicked, this, &GroundStateDlg::SyncFromKernel);
	connect(btnToKernel, &QAbstractButton::clicked, this, &GroundStateDlg::SyncToKernel);
	connect(btnOk, &QAbstractButton::clicked, this, &QDialog::accept);
}



/**
 * set a pointer to the main magdyn kernel
 */
void GroundStateDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn_kern = dyn;
	if(!m_dyn)
		SyncFromKernel();
}



/**
 * get spin configuration from kernel
 */
void GroundStateDlg::SyncFromKernel()
{
	if(!m_dyn_kern)
		return;

	m_dyn = *m_dyn_kern;

	m_sitestab->setSortingEnabled(false);

	for(const t_site& site : m_dyn->GetMagneticSites())
	{
		const t_vec_real& S = site.spin_dir_calc;
		const auto [ rho, phi, theta ] =  tl2::cart_to_sph<t_real>(S[0], S[1], S[2]);
		const auto [ u, v ] = tl2::sph_to_uv<t_real>(phi, theta);

		int row = m_sitestab->rowCount();
		m_sitestab->insertRow(row);

		m_sitestab->setItem(row, COL_GSSITE_NAME, new QTableWidgetItem(site.name.c_str()));
		m_sitestab->setItem(row, COL_GSSITE_PHI, new tl2::NumericTableWidgetItem<t_real>(u));
		m_sitestab->setItem(row, COL_GSSITE_THETA, new tl2::NumericTableWidgetItem<t_real>(v));
	}

	m_sitestab->setSortingEnabled(true);
}



/**
 * send spin configuration to kernel
 */
void GroundStateDlg::SyncToKernel()
{
	emit SpinsUpdated(&*m_dyn);
}



/**
 * dialog is closing
 */
void GroundStateDlg::closeEvent(QCloseEvent *)
{
	if(m_sett)
		m_sett->setValue("ground_state/geo", saveGeometry());
}
