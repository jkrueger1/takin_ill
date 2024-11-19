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
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QProgressBar>

#include <qcustomplot.h>
#include <vector>

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
	void ClearPlot(bool replot = true);
	void Plot();
	void PlotMouseMove(QMouseEvent* evt);
	void PlotMousePress(QMouseEvent* evt);

	void EnableCalculation(bool enable = true);
	void Calculate();


private:
	const t_magdyn *m_dyn{};         // main calculation kernel

	QCustomPlot *m_plot{};                    // plotter
	std::vector<QCPCurve*> m_curves{};        // plot curves
	std::vector<QVector<qreal>> m_Qs_data{};  // momentum transfer per band
	std::vector<QVector<qreal>> m_Bs_data{};  // berry curvature per band
	t_size m_Q_idx{};                         // index of dominant Q component
	t_real m_Q_min{}, m_Q_max{};              // range of dominant Q component
	t_real m_B_min{}, m_B_max{};              // range of berry curvature

	QDoubleSpinBox *m_Q_start[3]{};  // Q start coordinate
	QDoubleSpinBox *m_Q_end[3]{};    // Q end coordinate
	QSpinBox *m_num_Q{};             // number of Q coordinates

	QDoubleSpinBox *m_B_filter{};    // maximum B value
	QSpinBox *m_coords[2]{};         // berry curvature component indices
	QCheckBox *m_imag{};             // imaginary or real components?

	QPushButton *m_btnStartStop{};   // start/stop calculation
	bool m_calcEnabled{};            // enable calculations
	bool m_stopRequested{};          // stop running calculations

	QSettings *m_sett{};             // program settings
	QProgressBar *m_progress{};      // progress bar
	QLabel *m_status{};              // status bar
	QMenu *m_menuPlot{};             // context menu for plot
};


#endif
