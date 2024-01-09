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

#ifndef __GENPOS_DLG_H__
#define __GENPOS_DLG_H__

#include <vector>

#include <QDialog>
#include <QSettings>
#include "ui/ui_genpos.h"

#include "tlibs/math/linalg.h"
#include "libs/globals.h"


struct ScanPosition
{
	t_real_glob h, k, l, E, ki, kf;
};


class GenPosDlg : public QDialog, Ui::GenPosDlg
{ Q_OBJECT
	public:
		GenPosDlg(QWidget* pParent = nullptr, QSettings* pSett = nullptr);
		virtual ~GenPosDlg();

	protected:
		QSettings *m_pSettings = nullptr;

	protected slots:
		void GeneratePositions();
		void ButtonBoxClicked(QAbstractButton* pBtn);

	signals:
		void GeneratedPositions(const std::vector<ScanPosition>& pos);
};

#endif
