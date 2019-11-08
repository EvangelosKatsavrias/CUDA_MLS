//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#include<QApplication>
#include <QtGui>
#include <QtWidgets>


class OpenFileAction: public QAction
{

protected:
	QString m_fileName;


public:
	using QAction::QAction;

	QString get_fileName(  ) { return m_fileName; }

	QString openFile() {
		m_fileName = QFileDialog::getOpenFileName(0, "Open File", "", "Mesh Files(*.msh)"); 
		if ( m_fileName != "") {
			QFile file(m_fileName);
			if (!file.open(QIODevice::ReadOnly) ) {
				QMessageBox::critical(0, "Error", "Could not open the file.");
				return "";
			}
			file.close();
			return m_fileName;
		}
		return "";
	}

};


class SaveFileAction: public QAction
{

protected:
	QString m_fileName;


public:
	using QAction::QAction;

	QString get_fileName(  ) { return m_fileName; }

	QString saveFile() {
		m_fileName = QFileDialog::getSaveFileName(0, "Save File", "", "Mesh Files(*.msh)"); 
		if ( m_fileName != "") {
			QFile file(m_fileName);
			if (!file.open(QIODevice::WriteOnly) ) {
				QMessageBox::critical(0, "Error", "Could not write the file.");
				return "";
			}
			file.close();
			return m_fileName;
		}
		return "";
	}

};



