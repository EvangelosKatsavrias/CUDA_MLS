#ifndef UCSAXESLINESHEADER
#define UCSAXESLINESHEADER
#include"myQLine.h"
#include"myQArrow.h"


class UCSAxesLinesScene : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	myQLine*		m_axis_x;
	myQLine*		m_axis_y;
	myQLine*		m_axis_z;


public:
	explicit UCSAxesLinesScene(Qt3DCore::QEntity *rootEntity): m_rootEntity(rootEntity)
	{
		m_axis_x = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(1.f, 0.f, 0.f), QColor(Qt::red) );
		m_axis_y = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 1.f, 0.f), QColor(Qt::green) );
		m_axis_z = new myQLine( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 0.f, 1.f), QColor(Qt::blue) );
	}
	~UCSAxesLinesScene() {}


public slots:
void enableDraw(bool enabled)
{
    //m_axis_x_->setEnabled(enabled);
}

};



class UCSAxesScene : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	myQArrow*		m_axis_x;
	myQArrow*		m_axis_y;
	myQArrow*		m_axis_z;


public:
	explicit UCSAxesScene(Qt3DCore::QEntity *rootEntity): m_rootEntity(rootEntity)
	{
		m_axis_x = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(1.f, 0.f, 0.f), 1.f, QColor(Qt::red) );
		m_axis_y = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 1.f, 0.f), 1.f, QColor(Qt::green) );
		m_axis_z = new myQArrow( rootEntity, QVector3D(0.f, 0.f, 0.f),  QVector3D(0.f, 0.f, 1.f), 1.f, QColor(Qt::blue) );
	}
	~UCSAxesScene() {}


public slots:
void enableDraw(bool enabled)
{
    //m_axis_x_->setEnabled(enabled);
}

};



#endif
