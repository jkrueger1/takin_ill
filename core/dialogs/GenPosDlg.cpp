/**
 * Generate Scan Positions
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 8-jan-2024
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

#include "GenPosDlg.h"

#include "tlibs/phys/neutrons.h"
#include "tlibs/string/string.h"


using t_real = t_real_glob;

static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();


GenPosDlg::GenPosDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSettings(pSett)
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

		if(m_pSettings->contains("gen_pos/geo"))
			restoreGeometry(m_pSettings->value("gen_pos/geo").toByteArray());
	}

	// connections
	QObject::connect(btnGenerate, &QAbstractButton::clicked,
		this, &GenPosDlg::GeneratePositions);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &GenPosDlg::ButtonBoxClicked);
	QObject::connect(checkCKf, &QCheckBox::stateChanged, [this](int state)
	{
		if(state == Qt::Unchecked)
			labelFixedK->setText("ki (Å⁻¹):");
		else
			labelFixedK->setText("kf (Å⁻¹):");
	});
}


GenPosDlg::~GenPosDlg()
{
}


/**
 * generate the scan positions and send them to the listeners
 */
void GenPosDlg::GeneratePositions()
{
	t_real hi = tl::str_to_var<t_real>(editHi->text().toStdString());
	t_real ki = tl::str_to_var<t_real>(editKi->text().toStdString());
	t_real li = tl::str_to_var<t_real>(editLi->text().toStdString());
	t_real Ei = tl::str_to_var<t_real>(editEi->text().toStdString());
	t_real hf = tl::str_to_var<t_real>(editHf->text().toStdString());
	t_real kf = tl::str_to_var<t_real>(editKf->text().toStdString());
	t_real lf = tl::str_to_var<t_real>(editLf->text().toStdString());
	t_real Ef = tl::str_to_var<t_real>(editEf->text().toStdString());

	t_real kfix = tl::str_to_var<t_real>(editKfix->text().toStdString());
	bool kf_fixed = checkCKf->isChecked();

	std::vector<ScanPosition> positions;
	int steps = spinSteps->value();
	positions.reserve(steps);

	for(int step = 0; step < steps; ++step)
	{
		ScanPosition pos;

		pos.h = tl::lerp(hi, hf, t_real(step)/t_real(steps-1));
		pos.k = tl::lerp(ki, kf, t_real(step)/t_real(steps-1));
		pos.l = tl::lerp(li, lf, t_real(step)/t_real(steps-1));
		pos.E = tl::lerp(Ei, Ef, t_real(step)/t_real(steps-1));

		if(kf_fixed)
		{
			pos.kf = kfix;
			pos.ki = tl::get_other_k(pos.E*meV, kfix/angs, !kf_fixed) * angs;
		}
		else
		{
			pos.ki = kfix;
			pos.kf = tl::get_other_k(pos.E*meV, kfix/angs, !kf_fixed) * angs;
		}

		positions.emplace_back(std::move(pos));
	}

	emit GeneratedPositions(positions);
}


void GenPosDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("gen_pos/geo", saveGeometry());

		QDialog::accept();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}
}


#include "moc_GenPosDlg.cpp"
