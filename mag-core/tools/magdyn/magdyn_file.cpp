/**
 * magnetic dynamics -- saving of dispersion data
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

#include "magdyn.h"

#include <QtCore/QString>

#include <vector>
#include <array>



// --------------------------------------------------------------------------------
/**
 * save the plot as pdf
 */
void MagDynDlg::SavePlotFigure()
{
	if(!m_plot)
		return;

	QString dirLast = m_sett->value("dir", "").toString();
	QString filename = QFileDialog::getSaveFileName(
		this, "Save Figure", dirLast, "PDF Files (*.pdf)");
	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	m_plot->savePdf(filename);
}



/**
 * save the data for a dispersion direction
 */
void MagDynDlg::SaveDispersion(bool as_scr)
{
	QString dirLast = m_sett->value("dir", "").toString();

	QString filename;
	if(as_scr)
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data As Script", dirLast, "Py Files (*.py)");
	else
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data", dirLast, "Data Files (*.dat)");

	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	const t_real Q_start[]
	{
		(t_real)m_Q_start[0]->value(),
		(t_real)m_Q_start[1]->value(),
		(t_real)m_Q_start[2]->value(),
	};

	const t_real Q_end[]
	{
		(t_real)m_Q_end[0]->value(),
		(t_real)m_Q_end[1]->value(),
		(t_real)m_Q_end[2]->value(),
	};

	const t_size num_pts = m_num_points->value();

	// TODO
	bool stop_request = false;
	m_statusFixed->setText("Calculating dispersion.");
	m_dyn.SaveDispersion(filename.toStdString(),
		Q_start[0], Q_start[1], Q_start[2],
		Q_end[0], Q_end[1], Q_end[2],
		num_pts, g_num_threads, as_scr,
		&stop_request);
	m_statusFixed->setText("Ready.");
}



/**
 * save the data for multiple dispersion directions
 * given by the saved coordinates
 */
void MagDynDlg::SaveMultiDispersion(bool as_scr)
{
	if(!m_coordinatestab->rowCount())
	{
		QMessageBox::critical(this, "Magnetic Dynamics", "No Q coordinates available, "
			"please define them in the \"Coordinates\" tab.");
		return;
	}

	QString dirLast = m_sett->value("dir", "").toString();
	QString filename;
	if(as_scr)
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data As Script", dirLast, "Py Files (*.py)");
	else
		filename = QFileDialog::getSaveFileName(
			this, "Save Dispersion Data", dirLast, "Data Files (*.dat)");
	if(filename == "")
		return;
	m_sett->setValue("dir", QFileInfo(filename).path());

	const t_size num_pts = m_num_points->value();

	// get all Qs from the coordinates table
	std::vector<std::array<t_real, 3>> Qs;
	std::vector<std::string> Q_names;
	Qs.reserve(m_coordinatestab->rowCount());
	Q_names.reserve(m_coordinatestab->rowCount());

	for(int coord_row = 0; coord_row < m_coordinatestab->rowCount(); ++coord_row)
	{
		std::string name = m_coordinatestab->item(
			coord_row, COL_COORD_NAME)->text().toStdString();
		t_real h = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_coordinatestab->item(coord_row, COL_COORD_H))->GetValue();
		t_real k = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_coordinatestab->item(coord_row, COL_COORD_K))->GetValue();
		t_real l = static_cast<tl2::NumericTableWidgetItem<t_real>*>(
			m_coordinatestab->item(coord_row, COL_COORD_L))->GetValue();

		Q_names.emplace_back(std::move(name));
		Qs.emplace_back(std::array<t_real, 3>{{ h, k, l }});
	}

	// TODO
	bool stop_request = false;
	m_statusFixed->setText("Calculating dispersion.");
	m_dyn.SaveMultiDispersion(filename.toStdString(),
		Qs, num_pts, g_num_threads, as_scr,
		&stop_request, &Q_names);
	m_statusFixed->setText("Ready.");
}
// --------------------------------------------------------------------------------
