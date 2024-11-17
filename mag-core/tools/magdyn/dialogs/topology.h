/**
 * magnetic dynamics -- topological calculations
 * @author Tobias Weber <tweber@ill.fr>
 * @date November 2024
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

#ifndef __MAG_DYN_TOPO_DLG_H__
#define __MAG_DYN_TOPO_DLG_H__

#include <QtCore/QSettings>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenu>

#include <qcustomplot.h>

#include "gui_defs.h"




/**
 * topology dialog
 */
class TopologyDlg : public QDialog
{ Q_OBJECT
public:
	TopologyDlg(QWidget *parent, QSettings *sett = nullptr);
	virtual ~TopologyDlg();

	TopologyDlg(const TopologyDlg&) = delete;
	TopologyDlg& operator=(const TopologyDlg&) = delete;

	void SetKernel(const t_magdyn* dyn);


protected:
	virtual void accept() override;

	void RescalePlot();
	void PlotMouseMove(QMouseEvent* evt);
	void PlotMousePress(QMouseEvent* evt);


private:
	const t_magdyn *m_dyn{};  // main calculation kernel
	QCustomPlot *m_plot{};    // plotter

	QSettings *m_sett{};      // program settings
	QLabel *m_status{};       // status bar
	QMenu *m_menuPlot{};      // context menu for plot
};


#endif
