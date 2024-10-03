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

#ifndef __MAG_DYN_GROUNDSTATE_DLG_H__
#define __MAG_DYN_GROUNDSTATE_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QTableWidget>

#include "gui_defs.h"




/**
 * ground state dialog
 */
class GroundStateDlg : public QDialog
{ Q_OBJECT
public:
	GroundStateDlg(QWidget *parent, QSettings *sett = nullptr);
	~GroundStateDlg() = default;
	GroundStateDlg(const GroundStateDlg&) = delete;
	GroundStateDlg& operator=(const GroundStateDlg&) = delete;

	void SetKernel(const t_magdyn* dyn);


protected:
	void SyncFromKernel();
	void SyncToKernel();

	virtual void closeEvent(QCloseEvent *) override;


private:
	const t_magdyn *m_dyn_kern{};     // main kernel
	std::optional<t_magdyn> m_dyn{};  // local copy to work on

	QSettings *m_sett{};
	QTableWidget *m_sitestab{};


signals:
	void SpinsUpdated(const t_magdyn* dyn);
};


#endif
