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

#include "ElasticDlg.h"

#include "dialogs/FilePreviewDlg.h"

#include "tlibs/phys/neutrons.h"
#include "tlibs/file/loadinstr.h"
#include "tlibs/string/string.h"
#include "tlibs/string/spec_char.h"


// position table columns
#define POSTAB_H     0
#define POSTAB_K     1
#define POSTAB_L     2
#define POSTAB_KI    3
#define POSTAB_KF    4

#define POSTAB_COLS  5


using t_real = t_real_glob;
using t_mat = ublas::matrix<t_real>;
using t_vec = ublas::vector<t_real>;

static const tl::t_length_si<t_real> angs = tl::get_one_angstrom<t_real>();
static const tl::t_energy_si<t_real> meV = tl::get_one_meV<t_real>();
static const tl::t_angle_si<t_real> rads = tl::get_one_radian<t_real>();


ElasticDlg::ElasticDlg(QWidget* pParent, QSettings* pSett)
	: QDialog(pParent), m_pSettings(pSett)
{
	this->setupUi(this);

	if(m_pSettings)
	{
		QFont font;
		if(m_pSettings->contains("main/font_gen") &&
			font.fromString(m_pSettings->value("main/font_gen", "").toString()))
		{
			setFont(font);
		}

		if(m_pSettings->contains("elastic_pos/geo"))
			restoreGeometry(m_pSettings->value("elastic_pos/geo").toByteArray());
	}

	btnAddPosition->setIcon(load_icon("res/icons/list-add.svg"));
	btnDelPosition->setIcon(load_icon("res/icons/list-remove.svg"));

	// connections
	QObject::connect(btnAddPosition, &QAbstractButton::clicked,
		this, &ElasticDlg::AddPosition);
	QObject::connect(btnDelPosition, &QAbstractButton::clicked,
		this, &ElasticDlg::DelPosition);
	QObject::connect(btnGenPositions, &QAbstractButton::clicked,
		this, &ElasticDlg::GeneratePositions);
	// TODO: import scan files
	//QObject::connect(btnImportPositions, &QAbstractButton::clicked,
	//	this, &ElasticDlg::ImportPositions);
	QObject::connect(tablePositions, &QTableWidget::itemChanged,
		this, &ElasticDlg::CalcElasticPositions);
	QObject::connect(buttonBox, &QDialogButtonBox::clicked,
		this, &ElasticDlg::ButtonBoxClicked);
}


ElasticDlg::~ElasticDlg()
{
}


/**
 * calculate the elastic positions corresponding to the given inelastic ones
 */
void ElasticDlg::CalcElasticPositions()
{
	if(!m_bAllowCalculation)
		return;

	t_mat matB = tl::get_B(m_lattice, 1);
	t_mat matBinv;
	bool B_ok = tl::inverse(matB, matBinv);

	const std::wstring strAA = tl::get_spec_char_utf16("AA") +
		tl::get_spec_char_utf16("sup-") + tl::get_spec_char_utf16("sup1");

	std::wostringstream ostrResults1, ostrResults2;
	ostrResults1.precision(g_iPrec);
	ostrResults2.precision(g_iPrec);

	ostrResults1 << "<b>Elastic Positions Q' Corresponding to (Q, E) with kf' := ki:</b>";
	ostrResults1 << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResults1 << "<tr>";
	ostrResults1 << "<th><b>No.</b></th>";
	if(B_ok)
		ostrResults1 << "<th><b>Q (rlu)</b></th>";
	else
		ostrResults1 << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResults1 << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResults1 << "<th><b>E (meV)</b></th>";
	if(B_ok)
		ostrResults1 << "<th><b>Q' (rlu)</b></th>";
	else
		ostrResults1 << "<th><b>Q' (" << strAA << ")</b></th>";
	ostrResults1 << "<th><b>|Q'| (" << strAA << ")</b></th>";
	ostrResults1 << "</tr>";

	ostrResults2 << "<b>Elastic Positions Q'' Corresponding to (Q, E) with ki'' := kf:</b>";
	ostrResults2 << "<center><table border=\"1\" cellpadding=\"0\" width=\"95%\">";
	ostrResults2 << "<tr>";
	ostrResults2 << "<th><b>No.</b></th>";
	if(B_ok)
		ostrResults2 << "<th><b>Q (rlu)</b></th>";
	else
		ostrResults2 << "<th><b>Q (" << strAA << ")</b></th>";
	ostrResults2 << "<th><b>|Q| (" << strAA << ")</b></th>";
	ostrResults2 << "<th><b>E (meV)</b></th>";
	if(B_ok)
		ostrResults2 << "<th><b>Q'' (rlu)</b></th>";
	else
		ostrResults2 << "<th><b>Q'' (" << strAA << ")</b></th>";
	ostrResults2 << "<th><b>|Q''| (" << strAA << ")</b></th>";
	ostrResults2 << "</tr>";

	// iterate positions
	for(int row = 0; row < tablePositions->rowCount(); ++row)
	{
		try
		{
			QTableWidgetItem *item_h = tablePositions->item(row, POSTAB_H);
			QTableWidgetItem *item_k = tablePositions->item(row, POSTAB_K);
			QTableWidgetItem *item_l = tablePositions->item(row, POSTAB_L);
			QTableWidgetItem *item_ki = tablePositions->item(row, POSTAB_KI);
			QTableWidgetItem *item_kf = tablePositions->item(row, POSTAB_KF);

			if(!item_h || !item_k || !item_l || !item_ki || !item_kf)
				throw tl::Err("Invalid hkl, ki or kf.");

			// get coordinates
			t_real dH = tl::str_to_var_parse<t_real>(item_h->text().toStdString());
			t_real dK = tl::str_to_var_parse<t_real>(item_k->text().toStdString());
			t_real dL = tl::str_to_var_parse<t_real>(item_l->text().toStdString());
			t_real dKi = tl::str_to_var_parse<t_real>(item_ki->text().toStdString());
			t_real dKf = tl::str_to_var_parse<t_real>(item_kf->text().toStdString());

			t_real dMono2Theta = tl::get_mono_twotheta(dKi / angs, m_dMono*angs, m_bSenses[0]) / rads;
			t_real dAna2Theta = tl::get_mono_twotheta(dKf / angs, m_dAna*angs, m_bSenses[2]) / rads;
			t_real dSampleTheta, dSample2Theta;
			t_vec vecQ;

			if(tl::is_nan_or_inf<t_real>(dMono2Theta) || tl::is_nan_or_inf<t_real>(dAna2Theta))
				throw tl::Err("Invalid monochromator or analyser angle.");

			// get angles corresponding to given inelastic position
			tl::get_tas_angles(m_lattice,
				m_vec1, m_vec2,
				dKi, dKf,
				dH, dK, dL,
				m_bSenses[1],
				&dSampleTheta, &dSample2Theta,
				&vecQ);

			if(tl::is_nan_or_inf<t_real>(dSample2Theta)
				|| tl::is_nan_or_inf<t_real>(dSampleTheta))
				throw tl::Err("Invalid sample 2theta angle.");

			t_real dE = (tl::k2E(dKi / angs) - tl::k2E(dKf / angs)) / meV;
			t_vec vecQrlu = tl::prod_mv(matBinv, vecQ);

			tl::set_eps_0(dE, g_dEps);
			tl::set_eps_0(vecQ, g_dEps);
			tl::set_eps_0(vecQrlu, g_dEps);
			tl::set_eps_0(dSample2Theta, g_dEps);
			tl::set_eps_0(dSampleTheta, g_dEps);


			try
			{
				// get corresponding elastic analyser angle for kf' = ki
				t_real dAna2Theta_elast = tl::get_mono_twotheta(dKi/angs, m_dAna*angs, m_bSenses[2]) / rads;

				// find Q position for kf' = ki elastic position
				t_vec vecQ1;
				t_real dH1 = dH, dK1 = dK, dL1 = dL, dE1 = 0.;
				t_real dKi1 = dKi, dKf1 = dKi;
				tl::get_hkl_from_tas_angles<t_real>(m_lattice,
					m_vec1, m_vec2,
					m_dMono, m_dAna,
					dMono2Theta*0.5, dAna2Theta_elast*0.5, dSampleTheta, dSample2Theta,
					m_bSenses[0], m_bSenses[2], m_bSenses[1],
					&dH1, &dK1, &dL1,
					&dKi1, &dKf1, &dE1, 0,
					&vecQ1);

				if(tl::is_nan_or_inf<t_real>(dH1)
					|| tl::is_nan_or_inf<t_real>(dK1)
					|| tl::is_nan_or_inf<t_real>(dL1))
					throw tl::Err("Invalid h' k' l'.");

				t_vec vecQ1rlu = tl::prod_mv(matBinv, vecQ1);

				tl::set_eps_0(vecQ1, g_dEps);
				tl::set_eps_0(vecQ1rlu, g_dEps);

				for(t_real* d : { &dH1, &dK1, &dL1, &dKi1, &dKf1, &dE1 })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResults1 << "<tr>";
				ostrResults1 << "<td>" << row+1 << "</td>";
				if(B_ok)
					ostrResults1 << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResults1 << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResults1 << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResults1 << "<td>" << dE << "</td>";

				// print results for elastic kf' = ki position
				if(B_ok)
					ostrResults1 << "<td>" << vecQ1rlu[0] << ", " << vecQ1rlu[1] << ", " << vecQ1rlu[2] << "</td>";
				else
					ostrResults1 << "<td>" << vecQ1[0] << ", " << vecQ1[1] << ", " << vecQ1[2] << "</td>";
				ostrResults1 << "<td>" << tl::veclen(vecQ1) << "</td>";
				ostrResults1 << "</tr>";
			}
			catch(const std::exception& ex)
			{
				ostrResults1 << "<tr>";
				ostrResults1 << "<td>" << row+1 << "</td>";
				ostrResults1 << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResults1 << "</tr>";
			}


			try
			{
				// get corresponding elastic monochromator angle for ki'' = kf
				t_real dMono2Theta_elast = tl::get_mono_twotheta(dKf/angs, m_dMono*angs, m_bSenses[0]) / rads;

				// find Q position for kf' = ki elastic position
				t_vec vecQ2;
				t_real dH2 = dH, dK2 = dK, dL2 = dL, dE2 = 0.;
				t_real dKi2 = dKf, dKf2 = dKf;
				tl::get_hkl_from_tas_angles<t_real>(m_lattice,
					m_vec1, m_vec2,
					m_dMono, m_dAna,
					dMono2Theta_elast*0.5, dAna2Theta*0.5, dSampleTheta, dSample2Theta,
					m_bSenses[0], m_bSenses[2], m_bSenses[1],
					&dH2, &dK2, &dL2,
					&dKi2, &dKf2, &dE2, 0,
					&vecQ2);

				if(tl::is_nan_or_inf<t_real>(dH2)
					|| tl::is_nan_or_inf<t_real>(dK2)
					|| tl::is_nan_or_inf<t_real>(dL2))
					throw tl::Err("Invalid h'' k'' l''.");

				t_vec vecQ2rlu = tl::prod_mv(matBinv, vecQ2);

				tl::set_eps_0(vecQ2, g_dEps);
				tl::set_eps_0(vecQ2rlu, g_dEps);

				for(t_real* d : { &dH2, &dK2, &dL2, &dKi2, &dKf2, &dE2 })
					tl::set_eps_0(*d, g_dEps);


				// print inelastic position
				ostrResults2 << "<tr>";
				ostrResults2 << "<td>" << row+1 << "</td>";
				if(B_ok)
					ostrResults2 << "<td>" << vecQrlu[0] << ", " << vecQrlu[1] << ", " << vecQrlu[2] << "</td>";
				else
					ostrResults2 << "<td>" << vecQ[0] << ", " << vecQ[1] << ", " << vecQ[2] << "</td>";
				ostrResults2 << "<td>" << tl::veclen(vecQ) << "</td>";
				ostrResults2 << "<td>" << dE << "</td>";

				// print results for elastic ki'' = kf position
				if(B_ok)
					ostrResults2 << "<td>" << vecQ2rlu[0] << ", " << vecQ2rlu[1] << ", " << vecQ2rlu[2] << "</td>";
				else
					ostrResults2 << "<td>" << vecQ2[0] << ", " << vecQ2[1] << ", " << vecQ2[2] << "</td>";
				ostrResults2 << "<td>" << tl::veclen(vecQ2) << "</td>";
				ostrResults2 << "</tr>";
			}
			catch(const std::exception& ex)
			{
				ostrResults2 << "<tr>";
				ostrResults2 << "<td>" << row+1 << "</td>";
				ostrResults2 << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				ostrResults2 << "</tr>";
			}
		}
		catch(const std::exception& ex)
		{
			for(std::wostream* ostr : { &ostrResults1, &ostrResults2 })
			{
				(*ostr) << "<tr>";
				(*ostr) << "<td>" << row+1 << "</td>";
				(*ostr) << "<td><font color=\"#ff0000\"><b>" << ex.what() << "</b></font></td>";
				(*ostr) << "</tr>";
			}
		}
	}

	ostrResults1 << "</table></center>";
	ostrResults2 << "</table></center>";

	std::wostringstream ostrResults;
	ostrResults << "<html><body>";
	ostrResults << "<p>" << ostrResults1.str() << "</p>";
	ostrResults << "<br>";
	ostrResults << "<p>" << ostrResults2.str() << "</p>";
	ostrResults << "</body></html>";
	textResults->setHtml(QString::fromWCharArray(ostrResults.str().c_str()));
}


// ----------------------------------------------------------------------------
// positions table
// ----------------------------------------------------------------------------
/**
 * add an inelastic position to the table
 */
void ElasticDlg::AddPosition()
{
	// add a row
	int row = tablePositions->rowCount();
	tablePositions->setRowCount(row + 1);

	// insert items with default values
	for(int col=0; col<POSTAB_COLS; ++col)
		tablePositions->setItem(row, col, new QTableWidgetItem("0"));
}


/**
 * delete the selected inelastic positions from the table
 */
void ElasticDlg::DelPosition()
{
	int row = tablePositions->currentRow();
	if(row >= 0)
	{
		tablePositions->removeRow(row);
	}
	else
	{
		tablePositions->clearContents();
		tablePositions->setRowCount(0);
	}

	CalcElasticPositions();
}


/**
 * show dialog to generate positions
 */
void ElasticDlg::GeneratePositions()
{
	if(!m_pGenPosDlg)
	{
		m_pGenPosDlg = new GenPosDlg(this, m_pSettings);
		QObject::connect(m_pGenPosDlg, &GenPosDlg::GeneratedPositions,
			this, &ElasticDlg::GeneratedPositions);
	}

	focus_dlg(m_pGenPosDlg);
}


/**
 * add received generated positions
 */
void ElasticDlg::GeneratedPositions(const std::vector<ScanPosition>& positions)
{
	m_bAllowCalculation = false;

	// clear old positions
	tablePositions->clearContents();
	tablePositions->setRowCount(0);

	for(const ScanPosition& pos : positions)
	{
		AddPosition();

		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_H)->setText(
			tl::var_to_str(pos.h, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_K)->setText(
			tl::var_to_str(pos.k, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_L)->setText(
			tl::var_to_str(pos.l, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KI)->setText(
			tl::var_to_str(pos.ki, g_iPrec).c_str());
		tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KF)->setText(
			tl::var_to_str(pos.kf, g_iPrec).c_str());
	}

	// recalculate
	m_bAllowCalculation = true;
	CalcElasticPositions();
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// TODO: importing positions from scan files
// ----------------------------------------------------------------------------
/**
 * import positions from a scan file
 * TODO: for importing positions, set the lattice and scattering plane from the scan file
 *       and disable the setting of these properties from the main window
 */
void ElasticDlg::ImportPositions()
{
	std::vector<std::string> vecFiles = GetFiles();
	if(!vecFiles.size() || vecFiles[0]=="")
		return;

	m_bAllowCalculation = false;

	// clear old positions
	tablePositions->clearContents();
	tablePositions->setRowCount(0);

	// iterate scan files
	for(const std::string& strFile : vecFiles)
	{
		std::unique_ptr<tl::FileInstrBase<t_real>> ptrScan(
			tl::FileInstrBase<t_real>::LoadInstr(strFile.c_str()));
		if(!ptrScan)
		{
			tl::log_err("Invalid scan file: \"", strFile, "\".");
			continue;
		}

		// iterate scan positions and add them to the table
		for(std::size_t idx = 0; idx < ptrScan->GetScanCount(); ++idx)
		{
			AddPosition();

			auto pos = ptrScan->GetScanHKLKiKf(idx);
			tablePositions->item(tablePositions->rowCount() - 1, POSTAB_H)->setText(
				tl::var_to_str(std::get<0>(pos), g_iPrec).c_str());
			tablePositions->item(tablePositions->rowCount() - 1, POSTAB_K)->setText(
				tl::var_to_str(std::get<1>(pos), g_iPrec).c_str());
			tablePositions->item(tablePositions->rowCount() - 1, POSTAB_L)->setText(
				tl::var_to_str(std::get<2>(pos), g_iPrec).c_str());
			tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KI)->setText(
				tl::var_to_str(std::get<3>(pos), g_iPrec).c_str());
			tablePositions->item(tablePositions->rowCount() - 1, POSTAB_KF)->setText(
				tl::var_to_str(std::get<4>(pos), g_iPrec).c_str());
		}
	}

	// recalculate
	m_bAllowCalculation = true;
	CalcElasticPositions();
}


/**
 * opens a file selection dialog to select scan files
 */
std::vector<std::string> ElasticDlg::GetFiles()
{
	std::vector<std::string> vecFiles;

	QFileDialog::Option fileopt = QFileDialog::Option(0);
	bool bShowPreview = true;
	QString strDirLast;
	if(m_pSettings)
	{
		bShowPreview = m_pSettings->value("main/dlg_previews", true).toBool();
		strDirLast = m_pSettings->value("elastic_pos/last_import_dir", ".").toString();
		if(!m_pSettings->value("main/native_dialogs", true).toBool())
			fileopt = QFileDialog::DontUseNativeDialog;
	}

	QFileDialog *pdlg = nullptr;
	if(bShowPreview)
	{
		pdlg = new FilePreviewDlg(this, "Import Data File...", m_pSettings);
		pdlg->setOptions(QFileDialog::DontUseNativeDialog);
	}
	else
	{
		pdlg = new QFileDialog(this, "Import Data File...");
		pdlg->setOptions(fileopt);
	}
	std::unique_ptr<QFileDialog> ptrdlg(pdlg);

	pdlg->setDirectory(strDirLast);
	pdlg->setFileMode(QFileDialog::ExistingFiles);
	pdlg->setViewMode(QFileDialog::Detail);
#if !defined NO_IOSTR
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.ng0 *.NG0 *.log *.LOG *.scn.gz *.SCN.GZ *.dat.gz *.DAT.GZ *.ng0.gz *.NG0.GZ *.log.gz *.LOG.GZ *.scn.bz2 *.SCN.BZ2 *.dat.bz2 *.DAT.BZ2 *.ng0.bz2 *.NG0.BZ2 *.log.bz2 *.LOG.BZ2);;All files (*.*)";
#else
	QString strFilter = "Data files (*.dat *.scn *.DAT *.SCN *.NG0 *.ng0 *.log *.LOG);;All files (*.*)";
#endif
	pdlg->setNameFilter(strFilter);
	if(!pdlg->exec())
		return vecFiles;
	if(!pdlg->selectedFiles().size())
		return vecFiles;

	vecFiles.reserve(pdlg->selectedFiles().size());
	for(int iFile=0; iFile<pdlg->selectedFiles().size(); ++iFile)
		vecFiles.push_back(pdlg->selectedFiles()[iFile].toStdString());

	if(m_pSettings && vecFiles[0] != "")
	{
		std::string strDir = tl::get_dir(vecFiles[0]);
		m_pSettings->setValue("elastic_pos/last_import_dir", QString(strDir.c_str()));
	}

	return vecFiles;
}
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// loading / saving
// ----------------------------------------------------------------------------
/**
 * save the configuration
 */
void ElasticDlg::Save(std::map<std::string, std::string>& mapConf, const std::string& strXmlRoot)
{
	// TODO
}


/**
 * load the configuration
 */
void ElasticDlg::Load(tl::Prop<std::string>& xml, const std::string& strXmlRoot)
{
	// TODO
}
// ----------------------------------------------------------------------------


void ElasticDlg::ButtonBoxClicked(QAbstractButton* pBtn)
{
	if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::AcceptRole)
	{
		if(m_pSettings)
			m_pSettings->setValue("elastic_pos/geo", saveGeometry());

		QDialog::accept();
	}
	else if(buttonBox->buttonRole(pBtn) == QDialogButtonBox::RejectRole)
	{
		reject();
	}
}


void ElasticDlg::showEvent(QShowEvent *pEvt)
{
	QDialog::showEvent(pEvt);
}


// ----------------------------------------------------------------------------
// setter
// ----------------------------------------------------------------------------
void ElasticDlg::SetLattice(const tl::Lattice<t_real_glob>& lattice)
{
	m_lattice = lattice;
}


void ElasticDlg::SetScatteringPlane(const ublas::vector<t_real_glob>& vec1, const ublas::vector<t_real_glob>& vec2)
{
	m_vec1 = vec1;
	m_vec2 = vec2;
}


void ElasticDlg::SetD(t_real_glob dMono, t_real_glob dAna)
{
	m_dMono = dMono;
	m_dAna = dAna;
}


void ElasticDlg::SetMonoSense(bool bSense)
{
	m_bSenses[0] = bSense;
}


void ElasticDlg::SetSampleSense(bool bSense)
{
	m_bSenses[1] = bSense;
}


void ElasticDlg::SetAnaSense(bool bSense)
{
	m_bSenses[2] = bSense;
}


void ElasticDlg::SetSenses(bool bM, bool bS, bool bA)
{
	m_bSenses[0] = bM;
	m_bSenses[1] = bS;
	m_bSenses[2] = bA;
}
// ----------------------------------------------------------------------------


#include "moc_ElasticDlg.cpp"
