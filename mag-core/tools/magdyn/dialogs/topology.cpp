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
#include <boost/scope_exit.hpp>
#include <boost/asio.hpp>
namespace asio = boost::asio;

#include <limits>
#include <mutex>
#include <memory>
#include <sstream>

#include <QtWidgets/QGridLayout>

#include "topology.h"
#include "helper.h"

#include "tlibs2/libs/algos.h"



/**
 * shows the topology dialog
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
	m_plot->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Expanding});

	// context menu for plotter
	m_menuPlot = new QMenu("Plotter", this);
	QAction *acRescalePlot = new QAction("Rescale Axes", m_menuPlot);
	m_menuPlot->addAction(acRescalePlot);

	// start and stop coordinates
	m_Q_start[0] = new QDoubleSpinBox(this);
	m_Q_start[1] = new QDoubleSpinBox(this);
	m_Q_start[2] = new QDoubleSpinBox(this);
	m_Q_end[0] = new QDoubleSpinBox(this);
	m_Q_end[1] = new QDoubleSpinBox(this);
	m_Q_end[2] = new QDoubleSpinBox(this);

	m_Q_start[0]->setToolTip("Dispersion initial momentum transfer, h_i (rlu).");
	m_Q_start[1]->setToolTip("Dispersion initial momentum transfer, k_i (rlu).");
	m_Q_start[2]->setToolTip("Dispersion initial momentum transfer, l_i (rlu).");
	m_Q_end[0]->setToolTip("Dispersion final momentum transfer, h_f (rlu).");
	m_Q_end[1]->setToolTip("Dispersion final momentum transfer, k_f (rlu).");
	m_Q_end[2]->setToolTip("Dispersion final momentum transfer, l_f (rlu).");

	static const char* hklPrefix[] = { "h = ", "k = ","l = ", };
	for(int i = 0; i < 3; ++i)
	{
		m_Q_start[i]->setDecimals(4);
		m_Q_start[i]->setMinimum(-99.9999);
		m_Q_start[i]->setMaximum(+99.9999);
		m_Q_start[i]->setSingleStep(0.01);
		m_Q_start[i]->setValue(i == 0 ? -1. : 0.);
		m_Q_start[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_start[i]->setPrefix(hklPrefix[i]);

		m_Q_end[i]->setDecimals(4);
		m_Q_end[i]->setMinimum(-99.9999);
		m_Q_end[i]->setMaximum(+99.9999);
		m_Q_end[i]->setSingleStep(0.01);
		m_Q_end[i]->setValue(i == 0 ? 1. : 0.);
		m_Q_end[i]->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
		m_Q_end[i]->setPrefix(hklPrefix[i]);
	}

	// number of Q points in the plot
	m_num_Q = new QSpinBox(this);
	m_num_Q->setMinimum(1);
	m_num_Q->setMaximum(99999);
	m_num_Q->setValue(64);
	m_num_Q->setSizePolicy(QSizePolicy{QSizePolicy::Expanding, QSizePolicy::Fixed});
	m_num_Q->setToolTip("Number of Q points in the plot.");

	// start/stop button
	m_btnStartStop = new QPushButton("Calculate", this);

	// close button
	QPushButton *btnOk = new QPushButton("OK", this);

	// progress bar
	m_progress = new QProgressBar(this);
	m_progress->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// status bar
	m_status = new QLabel(this);
	m_status->setAlignment(Qt::AlignVCenter | Qt::AlignLeft);
	m_status->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

	// component grid
	auto grid = new QGridLayout(this);
	grid->setSpacing(4);
	grid->setContentsMargins(8, 8, 8, 8);

	int y = 0;
	grid->addWidget(m_plot, y++, 0, 1, 4);
	grid->addWidget(new QLabel("Start Q (rlu):", this), y, 0, 1, 1);
	grid->addWidget(m_Q_start[0], y, 1, 1, 1);
	grid->addWidget(m_Q_start[1], y, 2, 1, 1);
	grid->addWidget(m_Q_start[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("End Q (rlu):", this), y, 0, 1, 1);
	grid->addWidget(m_Q_end[0], y, 1, 1, 1);
	grid->addWidget(m_Q_end[1], y, 2, 1, 1);
	grid->addWidget(m_Q_end[2], y++, 3, 1, 1);
	grid->addWidget(new QLabel("Q Count:", this), y, 0, 1, 1);
	grid->addWidget(m_num_Q, y,1,1,1);
	grid->addWidget(m_btnStartStop, y, 2, 1, 1);
	grid->addWidget(btnOk, y++, 3, 1, 1);
	grid->addWidget(m_progress, y++, 0, 1, 4);
	grid->addWidget(m_status, y++, 0, 1, 4);

	// restore settings
	if(m_sett && m_sett->contains("topology/geo"))
		restoreGeometry(m_sett->value("topology/geo").toByteArray());
	else
		resize(640, 480);

	// connections
	connect(m_plot, &QCustomPlot::mouseMove, this, &TopologyDlg::PlotMouseMove);
	connect(m_plot, &QCustomPlot::mousePress, this, &TopologyDlg::PlotMousePress);
	connect(acRescalePlot, &QAction::triggered, this, &TopologyDlg::RescalePlot);
	connect(m_btnStartStop, &QAbstractButton::clicked, [this]()
	{
		// behaves as start or stop button?
		if(m_calcEnabled)
			Calculate();
		else
			m_stopRequested = true;
	});

	connect(btnOk, &QAbstractButton::clicked, this, &QDialog::accept);

	EnableCalculation();
}



TopologyDlg::~TopologyDlg()
{
}



/**
 * plot the berry curvature
 */
void TopologyDlg::Plot()
{
	if(!m_plot)
		return;

	m_plot->clearPlottables();
	m_curves.clear();

	// plot berry curvatures per band
	const t_size num_bands = m_Bs_data.size();
	for(t_size band = 0; band < num_bands; ++band)
	{
		QCPCurve *curve = new QCPCurve(m_plot->xAxis, m_plot->yAxis);

		// colour
		QPen pen = curve->pen();
		int col[3] = {
			int(std::lerp(1., 0., t_real(band) / t_real(num_bands - 1)) * 255.),
			0x00,
			int(std::lerp(0., 1., t_real(band) / t_real(num_bands - 1)) * 255.),
		};

		//get_colour<int>(g_colPlot, col);
		const QColor colFull(col[0], col[1], col[2]);
		pen.setColor(colFull);
		pen.setWidthF(2.);

		curve->setPen(pen);
		curve->setLineStyle(QCPCurve::lsLine);
		curve->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssNone, 1));
		curve->setAntialiased(true);
		curve->setData(m_Qs_data[band], m_Bs_data[band]);

		m_curves.push_back(curve);
	}

	// set labels
	const char* Q_label[]{ "h (rlu)", "k (rlu)", "l (rlu)" };
	m_plot->xAxis->setLabel(QString("Momentum Transfer ") + Q_label[m_Q_idx]);

	// set ranges
	m_plot->xAxis->setRange(m_Q_min, m_Q_max);
	m_plot->yAxis->setRange(m_B_min, m_B_max);

	m_plot->replot();
}



/**
 * calculate the berry curvature
 */
void TopologyDlg::Calculate()
{
	if(!m_dyn)
		return;

	BOOST_SCOPE_EXIT(this_)
	{
		this_->EnableCalculation(true);
	} BOOST_SCOPE_EXIT_END
	EnableCalculation(false);

	ClearPlot(false);

	// get coordinates
	t_vec_real Q_start = tl2::create<t_vec_real>(
	{
		m_Q_start[0]->value(),
		m_Q_start[1]->value(),
		m_Q_start[2]->value(),
	});

	t_vec_real Q_end = tl2::create<t_vec_real>(
	{
		m_Q_end[0]->value(),
		m_Q_end[1]->value(),
		m_Q_end[2]->value(),
	});

	// get Q component with maximum range
	t_vec_real Q_range = Q_end - Q_start;
	m_Q_idx = 0;
	if(std::abs(Q_range[1]) > std::abs(Q_range[m_Q_idx]))
		m_Q_idx = 1;
	if(std::abs(Q_range[2]) > std::abs(Q_range[m_Q_idx]))
		m_Q_idx = 2;

	// keep the scanned Q component in ascending order
	if(Q_start[m_Q_idx] > Q_end[m_Q_idx])
		std::swap(Q_start, Q_end);

	// get settings
	t_size Q_count = m_num_Q->value();
	t_real delta = 1e-12; //g_eps;
	std::vector<t_size> *perm = nullptr;
	t_size dim1 = 0, dim2 = 1;
	bool evecs_ortho = true;
	bool show_imag_comp = false;
	t_real max_curv = 100.;

	// calculate berry curvature
	t_magdyn dyn = *m_dyn;
	dyn.SetUniteDegenerateEnergies(false);

	// tread pool and mutex to protect m_Qs_data and m_Bs_data
	asio::thread_pool pool{g_num_threads};
	std::mutex mtx;

	m_stopRequested = false;
	m_progress->setMinimum(0);
	m_progress->setMaximum(Q_count);
	m_progress->setValue(0);
	m_status->setText(QString("Starting calculation using %1 thread(s).").arg(g_num_threads));
	tl2::Stopwatch<t_real> stopwatch;
	stopwatch.start();

	// create calculation tasks
	using t_task = std::packaged_task<void()>;
	using t_taskptr = std::shared_ptr<t_task>;
	std::vector<t_taskptr> tasks;
	tasks.reserve(Q_count);

	for(t_size i = 0; i < Q_count; ++i)
	{
		auto task = [this, &mtx, &dyn, &Q_start, &Q_end, i, Q_count, delta,
			perm, dim1, dim2, evecs_ortho, show_imag_comp, max_curv]()
		{
			const t_vec_real Q = Q_count > 1
				? tl2::lerp(Q_start, Q_end, t_real(i) / t_real(Q_count - 1))
				: Q_start;

			std::vector<t_cplx> curvs = dyn.CalcBerryCurvatures(
				Q, delta, perm, dim1, dim2, evecs_ortho);

			std::lock_guard<std::mutex> _lck{mtx};

			if(curvs.size() > m_Bs_data.size())
				m_Bs_data.resize(curvs.size());
			if(curvs.size() > m_Qs_data.size())
				m_Qs_data.resize(curvs.size());

			for(t_size band = 0; band < curvs.size(); ++band)
			{
				t_real berry_comp = show_imag_comp
					? curvs[band].imag()
					: curvs[band].real();

				if(std::abs(berry_comp) > max_curv)
					continue;

				m_Qs_data[band].push_back(Q[m_Q_idx]);
				m_Bs_data[band].push_back(berry_comp);
			}
		};

		t_taskptr taskptr = std::make_shared<t_task>(task);
		tasks.push_back(taskptr);
		asio::post(pool, [taskptr]() { (*taskptr)(); });
	}

	m_status->setText(QString("Calculating in %1 thread(s)...").arg(g_num_threads));

	// get results from tasks
	for(std::size_t task_idx = 0; task_idx < tasks.size(); ++task_idx)
	{
		t_taskptr task = tasks[task_idx];

		qApp->processEvents();  // process events to see if the stop button was clicked
		if(m_stopRequested)
		{
			pool.stop();
			break;
		}

		task->get_future().get();
		m_progress->setValue(task_idx + 1);
	}

	pool.join();
	stopwatch.stop();

	// show elapsed time
	std::ostringstream ostrMsg;
	ostrMsg.precision(g_prec_gui);
	ostrMsg << "Calculation";
	if(m_stopRequested)
		ostrMsg << " stopped ";
	else
		ostrMsg << " finished ";
	ostrMsg << "after " << stopwatch.GetDur() << " s.";
	m_status->setText(ostrMsg.str().c_str());

	// sort data by Q
	auto sort_data = [](QVector<qreal>& Qvec, QVector<qreal>& Bvec)
	{
		// sort vectors by Q component
		std::vector<std::size_t> perm = tl2::get_perm(Qvec.size(),
			[&Qvec, &Bvec](std::size_t idx1, std::size_t idx2) -> bool
		{
			return Qvec[idx1] < Qvec[idx2];
		});

		Qvec = tl2::reorder(Qvec, perm);
		Bvec = tl2::reorder(Bvec, perm);
	};

	for(t_size band = 0; band < m_Bs_data.size(); ++band)
		sort_data(m_Qs_data[band], m_Bs_data[band]);

	// data ranges
	m_Q_min = Q_start[m_Q_idx];
	m_Q_max = Q_end[m_Q_idx];
	m_B_min = std::numeric_limits<t_real>::max();
	m_B_max = -m_B_min;

	// get berry curvature range
	for(const QVector<t_real>& Bs_data : m_Bs_data)
	{
		auto [min_B_iter, max_B_iter] = std::minmax_element(Bs_data.begin(), Bs_data.end());
		if(min_B_iter != Bs_data.end() && max_B_iter != Bs_data.end())
		{
			t_real B_range = *max_B_iter - *min_B_iter;

			m_B_max = std::max(m_B_max, *max_B_iter + B_range*0.05);
			m_B_min = std::min(m_B_min, *min_B_iter - B_range*0.05);
		}
	}

	Plot();
}



/**
 * clears the dispersion graph
 */
void TopologyDlg::ClearPlot(bool replot)
{
	m_curves.clear();

	if(m_plot)
	{
		m_plot->clearPlottables();
		if(replot)
			m_plot->replot();
	}

	m_Qs_data.clear();
	m_Bs_data.clear();
	m_Q_idx = 0;
}



/**
 * show current cursor coordinates
 */
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



/**
 * show plot context menu
 */
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
 * toggle between "calculate" and "stop" button
 */
void TopologyDlg::EnableCalculation(bool enable)
{
	m_calcEnabled = enable;

	if(enable)
	{
		m_btnStartStop->setText("Calculate");
		m_btnStartStop->setToolTip("Start calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-start"));
	}
	else
	{
		m_btnStartStop->setText("Stop");
		m_btnStartStop->setToolTip("Stop running calculation.");
		m_btnStartStop->setIcon(QIcon::fromTheme("media-playback-stop"));
	}
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
