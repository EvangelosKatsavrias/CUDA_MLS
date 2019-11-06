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



