/**
 * Elastic Positions Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 5-jan-2024
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2024  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * Copyright (C) 2013-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; version 2 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 * ----------------------------------------------------------------------------
 */

#include "ElasticDlg.h"
#include "tlibs/phys/neutrons.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"


// position table columns
#define POSTAB_H     0
#define POSTAB_K     1
#define POSTAB_L     2
#define POSTAB_KI    3
#define POSTAB_KF    4

#define POSTAB_COLS  5


using t_real = t_real_glob;
static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


ElasticDlg::ElasticDlg(QWidget* pParent, QSettings* pSett) : QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") &&
			font.fromString(m_pSettings->value("main/font_gen", "").toString()))
		{
			setFont(font);
		}
	}

	btnAddPosition->setIcon(load_icon("res/icons/list-add.svg"));
	btnDelPosition->setIcon(load_icon("res/icons/list-remove.svg"));

	// connections
	QObject::connect(btnAddPosition, &QAbstractButton::clicked,
		this, &ElasticDlg::AddPosition);
	QObject::connect(btnDelPosition, &QAbstractButton::clicked,
		this, &ElasticDlg::DelPosition);
	//QObject::connect(tablePositions, &QTableWidget::itemChanged,
	//	this, &ElasticDlg::CalcElasticPositions);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &ElasticDlg::ButtonBoxClicked);

	if(m_pSettings)
	{
		if(m_pSettings->contains("elastic_pos/geo"))
			restoreGeometry(m_pSettings->value("elastic_pos/geo").toByteArray());
	}
}


ElasticDlg::~ElasticDlg()
{
}


/**
 * calculate the elastic positions corresponding to the given inelastic ones
 */
void ElasticDlg::CalcElasticPositions()
{
	std::wostringstream ostrResults;

	// iterate positions
	for(int row = 0; row < tablePositions->rowCount(); ++row)
	{
		QTableWidgetItem *item_h = tablePositions->item(row, POSTAB_H);
		QTableWidgetItem *item_k = tablePositions->item(row, POSTAB_K);
		QTableWidgetItem *item_l = tablePositions->item(row, POSTAB_L);
		QTableWidgetItem *item_ki = tablePositions->item(row, POSTAB_KI);
		QTableWidgetItem *item_kf = tablePositions->item(row, POSTAB_KF);

		if(!item_h || !item_k || !item_l || !item_ki || !item_kf)
			continue;

		std::cout << row << ": " << item_kf->text().toStdString() << std::endl;

		// get coordinates
		t_real dH = tl::str_to_var_parse<t_real>(item_h->text().toStdString());
		t_real dK = tl::str_to_var_parse<t_real>(item_k->text().toStdString());
		t_real dL = tl::str_to_var_parse<t_real>(item_l->text().toStdString());
		t_real dKi = tl::str_to_var_parse<t_real>(item_ki->text().toStdString());
		t_real dKf = tl::str_to_var_parse<t_real>(item_kf->text().toStdString());

		t_real dSampleTheta, dSample2Theta;
		ublas::vector<t_real> vecQ;

		t_real dMono2Theta = tl::get_mono_twotheta(dKi / angs, m_dMono*angs, m_bSenses[0]) / rads;
		t_real dAna2Theta = tl::get_mono_twotheta(dKf / angs, m_dAna*angs, m_bSenses[2]) / rads;

		if(tl::is_nan_or_inf<t_real>(dMono2Theta) || tl::is_nan_or_inf<t_real>(dAna2Theta))
			continue;

		// get angles corresponding to given inelastic position
		tl::get_tas_angles(m_lattice,
			m_vec1, m_vec2,
			dKi, dKf,
			dH, dK, dL,
			m_bSenses[1],
			&dSampleTheta, &dSample2Theta,
			&vecQ);

		if(tl::is_nan_or_inf<t_real>(dSample2Theta)
			|| tl::is_nan_or_inf<t_real>(dSampleTheta))
			continue;

		t_real dE = (tl::k2E(dKi / angs) - tl::k2E(dKf / angs)) / meV;

		tl::set_eps_0(dE, g_dEps);
		tl::set_eps_0(vecQ, g_dEps);
		tl::set_eps_0(dSample2Theta, g_dEps);
		tl::set_eps_0(dSampleTheta, g_dEps);

		const std::wstring strAA = tl::get_spec_char_utf16("AA") +
			tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");

		// print inelastic position
		ostrResults << "Q_org = " << vecQ << " " << strAA << ", |Q| = "
			<< ublas::norm_2(vecQ) << " " << strAA
			<< ", E = " << dE << "\n";

		// get elastic analyser angle corresponding for kf' = ki
		t_real dAna2Theta_elast = tl::get_mono_twotheta(dKi/angs, m_dAna*angs, m_bSenses[2]) / rads;

		// find Q position for kf' = ki elastic position
		tl::get_hkl_from_tas_angles<t_real>(m_lattice,
			m_vec1, m_vec2,
			m_dMono, m_dAna,
			dMono2Theta*0.5, dAna2Theta_elast*0.5, dSampleTheta, dSample2Theta,
			m_bSenses[0], m_bSenses[2], m_bSenses[1],
			&dH, &dK, &dL,
			&dKi, &dKf, &dE, 0,
			&vecQ);

		if(tl::is_nan_or_inf<t_real>(dH)
			|| tl::is_nan_or_inf<t_real>(dK)
			|| tl::is_nan_or_inf<t_real>(dL))
			continue;

		for(t_real* d : { &dH, &dK, &dL, &dKi, &dKf, &dE })
			tl::set_eps_0(*d, g_dEps);

		// print results for elastic kf' = ki position
		ostrResults << "Q_new = " << vecQ << " " << strAA
			<< ", |Q| = " << ublas::norm_2(vecQ) << " " << strAA
			<< ", E = " << dE << "\n";
	}

	textResults->setText(QString::fromWCharArray(ostrResults.str().c_str()));
}


/**
 * add an inelastic position to the table
 */
void ElasticDlg::AddPosition()
{
	// add a row
	int row = tablePositions->rowCount();
	tablePositions->setRowCount(row + 1);

	// insert items with default values
	for(int col=0; col<POSTAB_COLS; ++col)
		tablePositions->setItem(row, col, new QTableWidgetItem("0"));
}


/**
 * delete the selected inelastic positions from the table
 */
void ElasticDlg::DelPosition()
{
	int row = tablePositions->currentRow();
	if(row >= 0)
		tablePositions->removeRow(row);
}


/**
 * save the configuration
 */
void ElasticDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	// TODO
}


/**
 * load the configuration
 */
void ElasticDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	// TODO
}


void ElasticDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("elastic_pos/geo", saveGeometry());

		QDialog::accept();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}
}


void ElasticDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


#include "moc_ElasticDlg.cpp"
