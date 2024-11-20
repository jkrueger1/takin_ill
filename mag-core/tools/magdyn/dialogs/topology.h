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
#include <QtWidgets/QTabWidget>
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

	// set kernel and Q path from main window
	void SetKernel(const t_magdyn* dyn);
	void SetDispersionQ(const t_vec_real& Qstart, const t_vec_real& Qend);


protected:
	virtual void accept() override;

	void ShowError(const char* msg);

	// ------------------------------------------------------------------------
	// berry curvature tab
	void RescaleBerryCurvaturePlot();
	void ClearBerryCurvaturePlot(bool replot = true);
	void PlotBerryCurvature();
	void SaveBerryCurvaturePlotFigure();
	void BerryCurvaturePlotMouseMove(QMouseEvent *evt);
	void BerryCurvaturePlotMousePress(QMouseEvent *evt);

	void EnableBerryCurvatureCalculation(bool enable = true);
	void CalculateBerryCurvature();
	void SaveBerryCurvatureData();

	void SetBerryCurvatureQ();
	// ------------------------------------------------------------------------


private:
	const t_magdyn *m_dyn{};            // main calculation kernel
	t_vec_real m_Qstart{}, m_Qend{};    // Qs from main window

	// ------------------------------------------------------------------------
	// main dialog
	QTabWidget *m_tabs{};               // tabs
	QSettings *m_sett{};                // program settings
	QLabel *m_status{};                 // status bar
	// ------------------------------------------------------------------------

	// ------------------------------------------------------------------------
	// berry curvature tab
	QCustomPlot *m_plot_bc{};                 // berry curvature plotter
	std::vector<QCPCurve*> m_curves_bc{};     // berry cyrvature plot curves
	std::vector<QVector<qreal>> m_Qs_data_bc{};  // momentum transfer per band
	std::vector<QVector<qreal>> m_Bs_data_bc{};  // berry curvature per band
	t_size m_Q_idx_bc{};                      // index of dominant Q component
	t_real m_Q_min_bc{}, m_Q_max_bc{};        // range of dominant Q component
	t_real m_B_min_bc{}, m_B_max_bc{};        // range of berry curvature

	QDoubleSpinBox *m_Q_start_bc[3]{};  // Q start coordinate
	QDoubleSpinBox *m_Q_end_bc[3]{};    // Q end coordinate
	QSpinBox *m_num_Q_bc{};             // number of Q coordinates

	QDoubleSpinBox *m_B_filter_bc{};    // maximum B value
	QSpinBox *m_coords_bc[2]{};         // berry curvature component indices
	QCheckBox *m_imag_bc{};             // imaginary or real components?

	QPushButton *m_btnStartStop_bc{};   // start/stop calculation
	bool m_calcEnabled_bc{};            // enable calculations
	bool m_stopRequested_bc{};          // stop running calculations

	QProgressBar *m_progress_bc{};      // progress bar
	QMenu *m_menuPlot_bc{};             // context menu for plot
	// ------------------------------------------------------------------------
};


#endif
