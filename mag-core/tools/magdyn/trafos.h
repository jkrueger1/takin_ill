/**
 * transformation calculator
 * @author Tobias Weber <tweber@ill.fr>
 * @date 29-dec-2022
 * @license GPLv3, see 'LICENSE' file
 * @note Forked on 7-sep-2023 from my privately developed "gl" project: https://github.com/t-weber/gl .
 */

#ifndef __MAGDYN_TRAFOCALC_H__
#define __MAGDYN_TRAFOCALC_H__

#include <QtCore/QSettings>
#include <QtWidgets/QDialog>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QDoubleSpinBox>


class TrafoCalculator : public QDialog
{ Q_OBJECT
public:
	TrafoCalculator(QWidget* pParent = nullptr, QSettings *sett = nullptr);
	virtual ~TrafoCalculator() = default;

	TrafoCalculator(const TrafoCalculator&) = delete;
	const TrafoCalculator& operator=(const TrafoCalculator&) = delete;


private:
	QSettings *m_sett{};

	QTextEdit *m_textRotation{};
	QDoubleSpinBox *m_spinAxis[3]{nullptr, nullptr, nullptr};
	QDoubleSpinBox *m_spinAngle{};
	QDoubleSpinBox *m_spinVecToRotate[3]{nullptr, nullptr, nullptr};


protected slots:
	virtual void accept() override;

	void CalculateRotation();
};


#endif
