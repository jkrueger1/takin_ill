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
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QTableWidgetItem>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QMessageBox>

#include <string>
#include <unordered_set>

#include "tlibs2/libs/qt/numerictablewidgetitem.h"



/**
 * columns of the sites table
 */
enum : int
{
	COL_SPIN_NAME = 0,             // label
	COL_SPIN_PHI, COL_SPIN_THETA,  // spin direction
	COL_SPIN_U, COL_SPIN_V,        // spin direction
	COL_SPIN_U_FIXED,              // don't minimise u?
	COL_SPIN_V_FIXED,              // don't minimise v?

	NUM_SPIN_COLS
};



/**
 * shows the 3d view of the magnetic structure
 */
GroundStateDlg::GroundStateDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Ground State Energy Minimisation");
	setSizeGripEnabled(true);

	m_spinstab = new QTableWidget(this);
	m_spinstab->setShowGrid(true);
	m_spinstab->setAlternatingRowColors(true);
	m_spinstab->setSortingEnabled(true);
	m_spinstab->setMouseTracking(true);
	m_spinstab->setSelectionBehavior(QTableWidget::SelectRows);
	m_spinstab->setSelectionMode(QTableWidget::ContiguousSelection);
	m_spinstab->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	m_spinstab->verticalHeader()->setDefaultSectionSize(fontMetrics().lineSpacing() + 4);
	m_spinstab->verticalHeader()->setVisible(true);

	m_spinstab->setColumnCount(NUM_SPIN_COLS);

	m_spinstab->setHorizontalHeaderItem(COL_SPIN_NAME, new QTableWidgetItem{"Name"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_PHI, new QTableWidgetItem{"Spin φ"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_THETA, new QTableWidgetItem{"Spin θ"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_U, new QTableWidgetItem{"Spin u"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_V, new QTableWidgetItem{"Spin v"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_U_FIXED, new QTableWidgetItem{"Fix u"});
	m_spinstab->setHorizontalHeaderItem(COL_SPIN_V_FIXED, new QTableWidgetItem{"Fix v"});

	m_spinstab->setColumnWidth(COL_SPIN_NAME, 125);
	m_spinstab->setColumnWidth(COL_SPIN_PHI, 90);
	m_spinstab->setColumnWidth(COL_SPIN_THETA, 90);
	m_spinstab->setColumnWidth(COL_SPIN_U, 90);
	m_spinstab->setColumnWidth(COL_SPIN_V, 90);
	m_spinstab->setColumnWidth(COL_SPIN_U_FIXED, 50);
	m_spinstab->setColumnWidth(COL_SPIN_V_FIXED, 50);

	QPushButton *btnFromKernel = new QPushButton("Get Spins", this);
	QPushButton *btnToKernel = new QPushButton("Set Spins", this);
	QPushButton *btnMinimise = new QPushButton("Minimise", this);
	QPushButton *btnOk = new QPushButton("OK", this);

	int y = 0;
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(6, 6, 6, 6);
	grid->addWidget(m_spinstab, y++, 0, 1, 4);
	grid->addWidget(btnFromKernel, y, 0, 1, 1);
	grid->addWidget(btnToKernel, y, 1, 1, 1);
	grid->addWidget(btnMinimise, y, 2, 1, 1);
	grid->addWidget(btnOk, y++, 3, 1, 1);

	if(m_sett && m_sett->contains("ground_state/geo"))
		restoreGeometry(m_sett->value("ground_state/geo").toByteArray());
	else
		resize(800, 800);

	connect(btnFromKernel, &QAbstractButton::clicked, [this](){ GroundStateDlg::SyncFromKernel(); });
	connect(btnToKernel, &QAbstractButton::clicked, this, &GroundStateDlg::SyncToKernel);
	connect(btnMinimise, &QAbstractButton::clicked, this, &GroundStateDlg::Minimise);
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
void GroundStateDlg::SyncFromKernel(const t_magdyn *dyn,
	const std::unordered_set<std::string> *fixed_spins)
{
	// if no kernel is given, use the main one
	if(!dyn)
		dyn = m_dyn_kern;
	if(!dyn)
		return;

	// copy the data (if it's not the same structure anyway)
	if(dyn != &*m_dyn)
		m_dyn = *dyn;

	m_spinstab->setSortingEnabled(false);
	m_spinstab->clearContents();
	m_spinstab->setRowCount(0);

	for(const t_site& site : m_dyn->GetMagneticSites())
	{
		const t_vec_real& S = site.spin_dir_calc;
		const auto [ rho, phi, theta ] =  tl2::cart_to_sph<t_real>(S[0], S[1], S[2]);
		const auto [ u, v ] = tl2::sph_to_uv<t_real>(phi, theta);

		int row = m_spinstab->rowCount();
		m_spinstab->insertRow(row);

		QTableWidgetItem *item_name = new QTableWidgetItem(site.name.c_str());
		item_name->setFlags(item_name->flags() & ~Qt::ItemIsEditable);

		QCheckBox* u_fixed = new QCheckBox(this);
		QCheckBox* v_fixed = new QCheckBox(this);
		u_fixed->setChecked(false);
		v_fixed->setChecked(false);

		m_spinstab->setItem(row, COL_SPIN_NAME, item_name);
		m_spinstab->setItem(row, COL_SPIN_PHI, new tl2::NumericTableWidgetItem<t_real>(tl2::r2d(phi)));
		m_spinstab->setItem(row, COL_SPIN_THETA, new tl2::NumericTableWidgetItem<t_real>(tl2::r2d(theta)));
		m_spinstab->setItem(row, COL_SPIN_U, new tl2::NumericTableWidgetItem<t_real>(u));
		m_spinstab->setItem(row, COL_SPIN_V, new tl2::NumericTableWidgetItem<t_real>(v));
		m_spinstab->setCellWidget(row, COL_SPIN_U_FIXED, u_fixed);
		m_spinstab->setCellWidget(row, COL_SPIN_V_FIXED, v_fixed);

		// keep fixed spins
		if(fixed_spins)
		{
			if(fixed_spins->find(site.name + "_phi") != fixed_spins->end())
				u_fixed->setChecked(true);
			if(fixed_spins->find(site.name + "_theta") != fixed_spins->end())
				v_fixed->setChecked(true);
		}
	}

	m_spinstab->setSortingEnabled(true);
}



/**
 * send spin configuration to kernel
 */
void GroundStateDlg::SyncToKernel()
{
	emit SpinsUpdated(&*m_dyn);
}



/**
 * minimise the ground state energy
 */
void GroundStateDlg::Minimise()
{
	// set fixed spin parameters
	std::unordered_set<std::string> fixed_spins;

	for(int row = 0; row < m_spinstab->rowCount(); ++row)
	{
		QCheckBox* u_fixed = reinterpret_cast<QCheckBox*>(m_spinstab->cellWidget(row, COL_SPIN_U_FIXED));
		if(u_fixed && u_fixed->isChecked())
		{
			std::string fixed_name = m_spinstab->item(row, COL_SPIN_NAME)->text().toStdString() + "_phi";
			std::cout << fixed_name << std::endl;
			fixed_spins.emplace(std::move(fixed_name));
		}

		QCheckBox* v_fixed = reinterpret_cast<QCheckBox*>(m_spinstab->cellWidget(row, COL_SPIN_V_FIXED));
		if(v_fixed && v_fixed->isChecked())
		{
			std::string fixed_name = m_spinstab->item(row, COL_SPIN_NAME)->text().toStdString() + "_theta";
			std::cout << fixed_name << std::endl;
			fixed_spins.emplace(std::move(fixed_name));
		}
	}

	// minimise
	if(!m_dyn->CalcGroundState(&fixed_spins, true))
		QMessageBox::critical(this, "Error", "Could not find minimum ground state energy.");

	// set new parameters
	SyncFromKernel(&*m_dyn, &fixed_spins);
}



/**
 * dialog is closing
 */
void GroundStateDlg::closeEvent(QCloseEvent *)
{
	if(m_sett)
		m_sett->setValue("ground_state/geo", saveGeometry());
}
