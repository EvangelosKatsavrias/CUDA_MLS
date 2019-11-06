#ifndef MYQMESHRENDEREDHEADER
#define MYQMESHRENDEREDHEADER

#include <QApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QMesh>
#include <Qt3DExtras/QPhongMaterial>
#include <Qt3DExtras/QPerVertexColorMaterial>
#include"myAlgorithms.h"
#include<set>
#include"Mesh.h"




template<class T>
class myQMeshRendered : public QObject
{

private:
	Qt3DCore::QEntity*	m_rootEntity;
	Qt3DCore::QEntity*	m_lineEntity;
	QByteArray		m_vertexBufferData;
	QByteArray		m_colorBufferData;
	T*			m_positions;
	QByteArray		m_indexBufferData;
	unsigned int*		m_indices;
	QColor			m_materialColor;
	Qt3DExtras::QPerVertexColorMaterial* m_material;
	Mesh<T, 3, int>*		m_mesh;


public:
	explicit myQMeshRendered(Qt3DCore::QEntity *rootEntity, Mesh<T, 3, int>* mesh, QColor color = QColor(Qt::red) ); 
	~myQMeshRendered() {}


public slots:
	void enableDraw(bool enabled) { m_lineEntity->setEnabled(enabled); }

};

template class myQMeshRendered<float>;
template class myQMeshRendered<double>;

#endif
