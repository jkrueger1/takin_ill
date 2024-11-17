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

#include "topology.h"
#include "defs.h"

#include <QtWidgets/QGridLayout>



/**
 * shows the 3d view of the magnetic structure
 */
TopologyDlg::TopologyDlg(QWidget *parent, QSettings *sett)
	: QDialog{parent}, m_sett{sett}
{
	setWindowTitle("Topology");
	setSizeGripEnabled(true);

	// plotter
	m_plot = new QCustomPlot(this);
	m_plot->xAxis->setLabel("Momentum Transfer Q (rlu)");
	m_plot->yAxis->setLabel("Berry Curvature B");
	m_plot->setInteraction(QCP::iRangeDrag, true);
	m_plot->setInteraction(QCP::iRangeZoom, true);
	m_plot->setSelectionRectMode(QCP::srmZoom);
	m_plot->setSizePolicy(QSizePolicy{
		QSizePolicy::Expanding, QSizePolicy::Expanding});

	// context menu for plotter
	m_menuPlot = new QMenu("Plotter", this);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot);
	m_menuPlot->addAction(acRescalePlot);

	QPushButton *btnOk = new QPushButton("OK", this);

	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);


	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	grid->addWidget(m_plot, y++, 0, 1, 4);
	grid->addWidget(btnOk, y++, 3, 1, 1);
	grid->addWidget(m_status, y, 0, 1, 4);

	if(m_sett && m_sett->contains("topology/geo"))
		restoreGeometry(m_sett->value("topology/geo").toByteArray());
	else
		resize(640, 480);

	// connections
	connect(m_plot, &QCustomPlot::mouseMove, this, &TopologyDlg::PlotMouseMove);
	connect(m_plot, &QCustomPlot::mousePress, this, &TopologyDlg::PlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &TopologyDlg::RescalePlot);
	connect(btnOk, &QAbstractButton::clicked, this, &QDialog::accept);
}



TopologyDlg::~TopologyDlg()
{
}



void TopologyDlg::PlotMouseMove(QMouseEvent* evt)
{
	if(!m_status)
		return;

	t_real Q = m_plot->xAxis->pixelToCoord(evt->pos().x());
	t_real berry = m_plot->yAxis->pixelToCoord(evt->pos().y());

	QString status("Q = %1 rlu, B = %2.");
	status = status.arg(Q, 0, 'g', g_prec_gui).arg(berry, 0, 'g', g_prec_gui);
	m_status->setText(status);
}



void TopologyDlg::PlotMousePress(QMouseEvent* evt)
{
	// show context menu
	if(evt->buttons() & Qt::RightButton)
	{
		if(!m_menuPlot)
			return;
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
		QPoint pos = evt->globalPos();
#else
		QPoint pos = evt->globalPosition().toPoint();
#endif
		m_menuPlot->popup(pos);
		evt->accept();
        }
}



/**
 * rescale plot axes to fit the content
 */
void TopologyDlg::RescalePlot()
{
	if(!m_plot)
		return;

	m_plot->rescaleAxes();
	m_plot->replot();
}



/**
 * set a pointer to the main magdyn kernel
 */
void TopologyDlg::SetKernel(const t_magdyn* dyn)
{
	m_dyn = dyn;
}



/**
 * dialog is closing
 */
void TopologyDlg::accept()
{
	if(m_sett)
		m_sett->setValue("topology/geo", saveGeometry());

	QDialog::accept();
}
