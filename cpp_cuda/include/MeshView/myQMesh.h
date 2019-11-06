#ifndef MYQMESHHEADER
#define MYQMESHHEADER

#include <QApplication>
#include <Qt3DCore/QEntity>
#include <Qt3DRender/QMesh>
#include <Qt3DExtras/QPhongMaterial>
#include"myAlgorithms.h"
#include<set>
#include"Mesh.h"

template<class T>
class myQMesh : public QObject
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
	Qt3DExtras::QPhongMaterial* m_material;
	size_t			m_numOfNodes;
	size_t			m_numOfConnectivities;
	std::set< std::pair<size_t, size_t >, unorderLess<size_t> > 	m_uniqueEdges;
	Mesh<T, 3, int>*		m_mesh;


public:
	explicit myQMesh(Qt3DCore::QEntity *rootEntity, Mesh<T, 3, int>* mesh, QColor color = QColor(Qt::red) ); 
	~myQMesh() {}


public slots:
	void enableDraw(bool enabled) { m_lineEntity->setEnabled(enabled); }

};

template class myQMesh<float>;
template class myQMesh<double>;

#endif
