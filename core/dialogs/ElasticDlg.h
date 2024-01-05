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

#ifndef __ELASTIC_DLG_H__
#define __ELASTIC_DLG_H__

#include <QDialog>
#include <QSettings>
#include "ui/ui_elastic.h"

#include "tlibs/phys/lattice.h"
#include "tlibs/math/linalg.h"
#include "tlibs/file/prop.h"
#include "libs/globals.h"
#include "libs/globals_qt.h"
#include "tools/taz/tasoptions.h"


class ElasticDlg : public QDialog, Ui::ElasticDlg
{ Q_OBJECT
	public:
		ElasticDlg(QWidget* pParent = nullptr, QSettings* pSett = nullptr);
		virtual ~ElasticDlg();

	protected:
		QSettings *m_pSettings = nullptr;

		ublas::vector<t_real_glob> m_vec1, m_vec2;
		tl::Lattice<t_real_glob> m_lattice;

		t_real_glob m_dAna = 3.355, m_dMono = 3.355;
		bool m_bSenses [3] = { false, true, false };

	public:
		void SetLattice(const tl::Lattice<t_real_glob>& lattice)
		{ m_lattice = lattice; }
		void SetScatteringPlane(const ublas::vector<t_real_glob>& vec1, const ublas::vector<t_real_glob>& vec2)
		{ m_vec1 = vec1; m_vec2 = vec2; }

		void SetD(t_real_glob dMono, t_real_glob dAna)
		{ m_dMono = dMono; m_dAna = dAna; }

		void SetMonoSense(bool bSense) { m_bSenses[0] = bSense; }
		void SetSampleSense(bool bSense) { m_bSenses[1] = bSense; }
		void SetAnaSense(bool bSense) { m_bSenses[2] = bSense; }
		void SetSenses(bool bM, bool bS, bool bA)
		{ m_bSenses[0] = bM; m_bSenses[1] = bS; m_bSenses[2] = bA; }

	public slots:
		void CalcElasticPositions();

	protected slots:
		void ButtonBoxClicked(QAbstractButton* pBtn);

	protected:
		void AddPosition();
		void DelPosition();
		virtual void showEvent(QShowEvent *pEvt) override;

	public:
		void Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot);
		void Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot);
};

#endif
