/**
 * Log Dialog
 * @author Tobias Weber <tobias.weber@tum.de>
 * @date 14-oct-2017
 * @license GPLv2
 *
 * ----------------------------------------------------------------------------
 * Takin (inelastic neutron scattering software package)
 * Copyright (C) 2017-2023  Tobias WEBER (Institut Laue-Langevin (ILL),
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

#ifndef __LOG_DLG_H__
#define __LOG_DLG_H__

#include <QDialog>
#include <QSettings>
#include <QFileSystemWatcher>
#include "ui/ui_log.h"
#include <memory>


class LogDlg : public QDialog, Ui::LogDlg
{ Q_OBJECT
	protected:
		QSettings *m_pSettings = 0;
		std::unique_ptr<QFileSystemWatcher> m_pFileWatcher;

	public:
		LogDlg(QWidget* pParent=0, QSettings *pSett=0,
			const std::string& strLogFile="");
		virtual ~LogDlg() = default;

	protected:
		virtual void accept() override;

	protected slots:
		void LogFileChanged(const QString&);
};

#endif
