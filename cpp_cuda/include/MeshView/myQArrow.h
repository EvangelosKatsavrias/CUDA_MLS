#ifndef MYQARROWHEADER
#define MYQARROWHEADER

#include <QApplication>
#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DCore/QTransform>

class myQArrow: public QObject
{
	private:
		Qt3DCore::QEntity*	m_rootEntity;
		Qt3DCore::QEntity*	m_coneEntity;
		Qt3DCore::QEntity*	m_cylinderEntity;
		QVector3D		m_start;
		QVector3D		m_end;
		float			m_scale;
		QColor			m_color;


	public:
		explicit myQArrow( Qt3DCore::QEntity *rootEntity, QVector3D start, QVector3D end, float scale = 1, QColor color = QColor(Qt::red) );
		~myQArrow() {}

};

#endif
