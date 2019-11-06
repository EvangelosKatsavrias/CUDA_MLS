#ifndef NODESHEADER
#define NODESHEADER

#include<vector>
#include<fstream>
#include"fileFunctions.h"
#include<iostream>
#include<algorithm>
#include<sstream>
#include<iterator>
#include<string>


template<class T=float>
class Nodes
{
protected:
	std::vector<T*>			m_data;
	size_t				m_numOfNodes;
	std::vector<int>		m_indexOfStackedComponents;
	std::vector<std::vector<T> >	m_stackOfComponents;


	size_t	create_stackedComponent 	(size_t componentIndex);
	bool	remove_stackedComponent 	(size_t componentIndex);
	int 	checkForAllocatedMemorySpace 	(size_t* componentsIndices, size_t numOfComponents);
	int 	set_defaultInput 		(size_t* componentsIndices, size_t* numOfNodes, size_t* numOfComponents);

	std::vector<T*>			m_nodeReference;

public:
	Nodes (size_t numOfComponents = 2, size_t numOfNodes = 0, T** data = 0);
	Nodes (std::string filePath);
	Nodes (size_t numOfComponents, std::vector<T>** components);
	Nodes (Nodes<T>* referenceNodes, size_t numOfComponents = 0, size_t numOfNodes = 0);
	Nodes (Nodes<T>& referenceNodes, size_t numOfComponents = 0, size_t numOfNodes = 0);
	~Nodes ();


	int	set_numOfComponents (size_t);
	int	set_numOfNodes (size_t);

	// Setting is done by assigning a memory pointer to the values, except the case where a reference to a vector is provided
	int	set_components    (				   size_t* componentsIndices=0, size_t numOfComponents=0);
	int	set_components    (T** newData,			   size_t* componentsIndices=0, size_t numOfComponents=0);
	int	set_components    (std::vector<T>** newData,	   size_t* componentsIndices=0, size_t numOfComponents=0);
	int	set_components    (std::vector<T>& newData,	   size_t  componentIndex);
	int	set_components 	  (std::istream& inStream,	   size_t* componentsIndices=0, size_t numOfComponents=0);
	int	set_components 	  (std::string filePath,	   size_t* componentsIndices=0, size_t numOfComponents=0);
	int	copySet_components(T** newData,			   size_t* componentsIndices=0, size_t numOfComponents=0);


	// Setting is done by assigning a memory pointer to the values, except the case where a reference to a vector is provided
	int	insert_components (				   size_t numOfComponents=1, int insertInPosition=-1);
	int	insert_components (T** newComponents,		   size_t numOfComponents=1, int insertInPosition=-1);
	int	insert_components (std::vector<T>** newComponents, size_t numOfComponents=1, int insertInPosition=-1);
	int	insert_components (std::vector<T>& newComponent,   			     int insertInPosition=-1);
	int	insert_components (std::istream& inStream, 	   size_t numOfComponents=1, int insertInPosition=-1);
	int	insert_components (std::string filePath, 	   size_t numOfComponents=1, int insertInPosition=-1);
	int	remove_components (size_t* componentsIndices, 	   size_t numOfComponents=1);
	int	remove_allComponents();
	int	copyInsert_components(T** newComponents,	   size_t numOfComponents=1, int insertInPosition=-1);


	// Setting is done by copying the values
	int	set_dataValues 	  (T** newData, 		   size_t* componentsIndices=0, size_t firstDataPoint=0, size_t numOfNodes=0, size_t numOfComponents=0);
	int	set_dataValues    (std::vector<T>** newData,	   size_t* componentsIndices=0, size_t firstDataPoint=0, size_t numOfNodes=0, size_t numOfComponents=0);
	int	set_dataValues    (std::vector<T>& newData,	   size_t  componentIndex=0,    size_t firstDataPoint=0, size_t numOfNodes=0);
	int	set_dataValues    (std::istream& inStream, 	   size_t* componentsIndices=0, size_t firstDataPoint=0, size_t numOfNodes=0, size_t numOfComponents=0);
	int	set_dataValues    (std::string filePath, 	   size_t* componentsIndices=0, size_t firstDataPoint=0, size_t numOfNodes=0, size_t numOfComponents=0);


	// Setting is done by copying the values
	int	insert_data 	  (T** newData, 		   size_t numOfNodes, int insertInPosition=-1);
	int	insert_data 	  (std::vector<T>** newData,	   size_t numOfNodes, int insertInPosition=-1);
	int	insert_data 	  (std::vector<T>&  newData,	   size_t numOfNodes, int insertInPosition=-1);
	int	insert_data 	  (std::istream& inStream, 	   size_t numOfNodes, int insertInPosition=-1);
	int	insert_data 	  (std::string filePath, 	   size_t numOfNodes, int insertInPosition=-1);
	int	remove_data 	  (size_t fromPosition,    	   size_t numOfNodes);



	// Data accessors
	std::string virtual 	get_typename() { return "Nodes";}
	size_t			get_numOfComponents () {return m_data.size();}
	size_t			get_numOfNodes () {return m_numOfNodes;}


	T* 	operator[]		(size_t componentIndex) { return get_component(componentIndex); }
	void	get_node		(size_t nodeIndex, T** array) { T** allComponents = get_allComponents(); for (int i = 0; i < m_data.size(); i++) array[i] = allComponents[i]+nodeIndex; }

	T*	get_component		(size_t componentIndex) { return m_data.at(componentIndex); }
	int	get_component		(std::ostream& outStream, size_t componentIndex =0);
	int	get_component		(std::string file, size_t componentIndex =0);

	int	get_components		(std::ostream& outStream, size_t* componentsIndices=0, size_t numOfComponents=0);
	int	get_components		(std::string file, size_t* componentsIndices=0, size_t numOfComponents=0);
	int	get_components		(T** componentsData, size_t* componentsIndices=0, size_t numOfComponents=0) { for (int index=0; index<numOfComponents; index++) componentsData[index] = m_data[componentsIndices[index]]; }

	T**	get_allComponents	( ) { return m_data.data(); }
	int	get_allComponents	(std::ostream& outStream);
	int	get_allComponents	(std::string file);

	int	get_dataValues		(std::ostream& outStream, size_t* componentsIndices=0, size_t firstNode=0, size_t numOfNodes=0, size_t numOfComponents=0);
	int	get_dataValues		(std::string file, size_t* componentsIndices=0, size_t firstNode=0, size_t numOfNodes=0, size_t numOfComponents=0);
};


//template class Nodes<float>;
//template class Nodes<double>;
//template class Nodes<Nodes<T>>;


template<typename T>
std::ostream& operator<< (std::ostream& out, Nodes<T>& data)
{
	T** allComps = data.get_allComponents();
	for (size_t pointIndex = 0; pointIndex < data.get_numOfNodes(); pointIndex++ )
	{
		for (size_t compIndex = 0; compIndex < data.get_numOfComponents(); compIndex++ )
			out << allComps[compIndex][pointIndex] << "\t\t";
		out << std::endl;
	}

	return out;
}


template<class T>
Nodes<T>::Nodes(size_t numOfComponents, size_t numOfNodes, T** data)
{
	set_numOfNodes(numOfNodes);
	insert_components(data, numOfComponents, 0);
	m_nodeReference.resize(numOfComponents);
}


template<class T>
Nodes<T>::Nodes(std::string filePath)
{
	size_t numOfComponents(0), numOfNodes(0);
	getFromTextFile_numOfColumnsAndLines(filePath, &numOfComponents, &numOfNodes);
	set_numOfNodes(numOfNodes);
	insert_components(filePath, numOfComponents, 0);
	m_nodeReference.resize( numOfComponents );
}


template<class T>
Nodes<T>::Nodes(size_t numOfComponents, std::vector<T>** components)
{
	set_numOfNodes(components[0]->size());
	insert_components(components, numOfComponents);
	m_nodeReference.resize( numOfComponents );
}


template<class T>
Nodes<T>::Nodes( Nodes<T>* referenceNodes, size_t numOfComponents, size_t numOfNodes)
{
	set_numOfNodes(referenceNodes->get_numOfNodes());
	if (numOfNodes) set_numOfNodes(numOfNodes);
	if (!numOfComponents) numOfComponents = referenceNodes->get_numOfComponents();
	insert_components ( referenceNodes->get_allComponents(), numOfComponents );
	m_nodeReference.resize(numOfComponents);
}


template<class T>
Nodes<T>::Nodes( Nodes<T>& referenceNodes, size_t numOfComponents, size_t numOfNodes)
{
	set_numOfNodes(referenceNodes.get_numOfNodes());
	if (numOfNodes) set_numOfNodes(numOfNodes);
	if (!numOfComponents) numOfComponents = referenceNodes.get_numOfComponents();
	copyInsert_components( referenceNodes.get_allComponents(), numOfComponents );
	m_nodeReference.resize(numOfComponents);
}


template<class T>
Nodes<T>::~Nodes()
{}


// ===========================================================================================

template<class T>
int Nodes<T>::set_numOfComponents(size_t numOfComponents)
{
	m_data.resize(numOfComponents, 0);
	m_indexOfStackedComponents.resize(numOfComponents, -1);
	m_nodeReference.resize(numOfComponents);
}


template<class T>
int Nodes<T>::set_numOfNodes(size_t newNumOfNodes)
{
	m_numOfNodes = newNumOfNodes;
/*
	auto minel = std::min_element(m_indexOfStackedComponents.begin(), m_indexOfStackedComponents.end());
	if ( *minel == -1 )
	{
		try { throw std::exception(); } catch (std::exception e) {std::cout << "The number of data points for some of the components is managed externally." << std::endl;}
		return 0;
	}
*/

	for (int componentIndex = 0; componentIndex < m_data.size(); componentIndex++)
	{
		size_t stackIndex = m_indexOfStackedComponents[componentIndex];
		if (stackIndex !=-1)
		{
			m_stackOfComponents.at(stackIndex).reserve(m_numOfNodes);
			m_data.at(componentIndex) = m_stackOfComponents.at(stackIndex).data();
		}
	}
	return 1;
}


// =======================================================================================

template<class T>
size_t Nodes<T>::create_stackedComponent(size_t componentIndex)
{
	int index = m_indexOfStackedComponents.at(componentIndex);
	if ( index != -1 ) { m_stackOfComponents.at(index).clear(); return index; }

	m_stackOfComponents.push_back(std::vector<T>());
	m_stackOfComponents.back().reserve(m_numOfNodes);
	m_indexOfStackedComponents.at(componentIndex) = m_stackOfComponents.size()-1;

	m_data.at(componentIndex) = m_stackOfComponents.back().data();
	return m_stackOfComponents.size()-1;
}


template<class T>
bool Nodes<T>::remove_stackedComponent(size_t componentIndex)
{
	if ( m_indexOfStackedComponents.at(componentIndex) != -1 )
	{
		m_stackOfComponents.erase( m_stackOfComponents.begin() + m_indexOfStackedComponents.at(componentIndex));
		m_indexOfStackedComponents.at(componentIndex) = -1;
		return 1;
	}

	return 0;
}


// ===================================================================================

template<class T>
int Nodes<T>::set_components (size_t* componentsIndices, size_t numOfComponents)
{
	size_t allComponents[get_numOfComponents()];
	if ( !componentsIndices || !numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) allComponents[i] = i; componentsIndices =allComponents; numOfComponents = get_numOfComponents(); }

	size_t componentIndex, stackIndex;
	for (size_t index=0; index < numOfComponents; index++) {
		componentIndex = componentsIndices[index];
		stackIndex = create_stackedComponent(componentIndex); 
		m_stackOfComponents[stackIndex].resize(get_numOfNodes()); 
	}

	return 1;
}


template<class T>
int Nodes<T>::set_components (T** newData, size_t* componentsIndices, size_t numOfComponents)
{
	size_t allComponents[get_numOfComponents()];
	if ( !componentsIndices || !numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) allComponents[i] = i; componentsIndices =allComponents; numOfComponents = get_numOfComponents(); }

	size_t componentIndex;
	for (size_t index=0; index < numOfComponents; index++){
		componentIndex = componentsIndices[index];
		remove_stackedComponent(componentIndex);
		m_data[componentIndex] = newData[index];}

	return 1;
}


template<class T>
int Nodes<T>::set_components (std::vector<T>** newData, size_t* componentsIndices, size_t numOfComponents)
{
	size_t allComponents[get_numOfComponents()];
	if ( !componentsIndices || !numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) allComponents[i] = i; componentsIndices =allComponents; numOfComponents = get_numOfComponents(); }

	T* passedData[numOfComponents];
	for (size_t componentIndex=0; componentIndex < numOfComponents; componentIndex++) 
		passedData[componentIndex] = newData[componentIndex]->data();

	return set_components (passedData, componentsIndices, numOfComponents);
}


template<class T>
int Nodes<T>::set_components(std::vector<T>& newData, size_t componentIndex)
{ set_dataValues( newData, componentIndex, size_t(0), get_numOfNodes() ); }


template<class T>
int Nodes<T>::set_components (std::istream& fileStream, size_t* componentsIndices, size_t numOfComponents)
{
	size_t allComponents[get_numOfComponents()];
	if ( !componentsIndices || !numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) allComponents[i] = i; componentsIndices =allComponents; numOfComponents = get_numOfComponents(); }

	size_t stackIndex;
	for (size_t index = 0; index < numOfComponents; index++)
	{	stackIndex = create_stackedComponent(componentsIndices[index]); m_stackOfComponents[stackIndex].resize(m_numOfNodes); }
	set_dataValues(fileStream, componentsIndices, 0, m_numOfNodes, numOfComponents);
}


template<class T>
int Nodes<T>::set_components (std::string filePath, size_t* componentsIndices, size_t numOfComponents)
{
	size_t allComponents[get_numOfComponents()];
	if ( !componentsIndices || !numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) allComponents[i] = i; componentsIndices =allComponents; numOfComponents = get_numOfComponents(); }

	size_t stackIndex;
	for (size_t index = 0; index < numOfComponents; index++)
	{	stackIndex = create_stackedComponent(componentsIndices[index]); m_stackOfComponents[stackIndex].resize(m_numOfNodes); }
	set_dataValues(filePath, componentsIndices, 0, m_numOfNodes, numOfComponents);
}



template<class T>
int Nodes<T>::copySet_components (T** newData, size_t* componentsIndices, size_t numOfComponents)
{ set_dataValues( newData, componentsIndices, 0, get_numOfNodes(), numOfComponents); }



// ====================================================================================================


template<class T>
int Nodes<T>::insert_components(size_t numOfComponents, int insertInPosition)
{
	if (insertInPosition==-1) insertInPosition = get_numOfComponents();

	m_data.insert(m_data.begin()+insertInPosition, numOfComponents, 0);
	m_indexOfStackedComponents.insert(m_indexOfStackedComponents.begin()+insertInPosition, numOfComponents, -1);
	for (int index = 0; index<numOfComponents; index++) create_stackedComponent(insertInPosition+index);

	m_nodeReference.resize(m_data.size());

}


template<class T>
int Nodes<T>::insert_components(T** newComponents, size_t numOfComponents, int insertInPosition)
{
	T* allEmpty[numOfComponents];
	if (!newComponents) { for (int index=0; index<numOfComponents; index++) allEmpty[index] = 0; newComponents =allEmpty; }
	if (insertInPosition==-1) insertInPosition = get_numOfComponents(); // to insert at the end

	m_data.insert(m_data.begin()+insertInPosition, newComponents, newComponents+numOfComponents);
	m_indexOfStackedComponents.insert(m_indexOfStackedComponents.begin()+insertInPosition, numOfComponents, -1);

	m_nodeReference.resize(m_data.size());
}


template<class T>
int Nodes<T>::insert_components(std::vector<T>** newComponents, size_t numOfComponents, int insertInPosition)
{
	T* passedData[numOfComponents];
	for (size_t componentIndex=0; componentIndex < numOfComponents; componentIndex++) 
		passedData[componentIndex] = newComponents[componentIndex]->data();

	insert_components(passedData, numOfComponents, insertInPosition);
}


template<class T>
int Nodes<T>::insert_components(std::vector<T>& newComponent, int insertInPosition)
{
	if (insertInPosition==-1) insertInPosition = get_numOfComponents();

	m_data.insert(m_data.begin()+insertInPosition, 1, 0);
	m_indexOfStackedComponents.insert(m_indexOfStackedComponents.begin()+insertInPosition, 1, -1);
	create_stackedComponent(insertInPosition);

	set_dataValues( newComponent, insertInPosition, 0, get_numOfNodes() );


	m_nodeReference.resize(m_data.size());

}


template<class T>
int Nodes<T>::insert_components(std::istream& fileStream, size_t numOfComponents, int insertInPosition)
{
	if (insertInPosition==-1) insertInPosition = get_numOfComponents();
	insert_components(numOfComponents, insertInPosition);
	size_t compIndices[numOfComponents]; for ( size_t index =0; index<numOfComponents; index++ ) compIndices[index] = insertInPosition +index;
	set_components(fileStream, compIndices, numOfComponents); 
}


template<class T>
int Nodes<T>::insert_components(std::string filePath, size_t numOfComponents, int insertInPosition)
{
	if (insertInPosition==-1) insertInPosition = get_numOfComponents();
	insert_components(numOfComponents, insertInPosition);
	size_t compIndices[numOfComponents]; for ( size_t index =0; index<numOfComponents; index++ ) compIndices[index] = insertInPosition +index;
	set_components(filePath, compIndices, numOfComponents); 
}


class max_value
{
public:
	bool operator() (int a, int b) {
		if ( a > b ) return true;
		else return false;
	}

};

template<class T>
int Nodes<T>::remove_components(size_t* componentsIndices, size_t numOfComponents)
{
	std::sort(componentsIndices, componentsIndices+numOfComponents, max_value() );

	for (int index = 0; index<numOfComponents; index++ ) {
		remove_stackedComponent(componentsIndices[index]);
		m_data.erase( m_data.begin() +componentsIndices[index] );
		m_indexOfStackedComponents.erase(m_indexOfStackedComponents.begin()+ componentsIndices[index] ); 
	}


	m_nodeReference.resize(m_data.size());

}


template<class T>
int Nodes<T>::remove_allComponents()
{
	size_t components[get_numOfComponents()];
	for (size_t i = 0; i< get_numOfComponents(); i++) components[i] = i;

	remove_components(components, get_numOfComponents());
}



template<class T>
int Nodes<T>::copyInsert_components(T** newComponents, size_t numOfComponents, int insertInPosition)
{
	if (insertInPosition==-1) insertInPosition = get_numOfComponents();
	size_t componentsIndices[numOfComponents];
	for (int index=0; index<numOfComponents; index++) componentsIndices[index] = insertInPosition+index;

	insert_components  ( numOfComponents, insertInPosition);
	copySet_components ( newComponents, componentsIndices, numOfComponents);
}



// =============================================================================================

template<class T>
int Nodes<T>::checkForAllocatedMemorySpace(size_t* componentsIndices, size_t numOfComponents)
{
	for (size_t index=0; index < numOfComponents; index++) {
		size_t componentIndex = componentsIndices[index];
		if ( !m_data[componentIndex] )
		{
			try{ throw std::exception();} catch(std::exception e) {std::cout << "Before setting data, set the allocated memory location or container for the corresponding components." << std::endl;} return 0; } 
	}
	return 1;
}


template<class T>
int Nodes<T>::set_defaultInput(size_t* componentsIndices, size_t* numOfNodes, size_t* numOfComponents)
{
	if ( !*numOfNodes) *numOfNodes = get_numOfNodes();
	if ( !componentsIndices || !*numOfComponents ) { for (size_t i = 0; i < get_numOfComponents(); i++) componentsIndices[i] = i; *numOfComponents = get_numOfComponents(); }
}


template<class T>
int Nodes<T>::set_dataValues(T** newData, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
	size_t tempPointer[1]; if (!componentsIndices) componentsIndices = tempPointer;
	set_defaultInput(componentsIndices, &numOfNodes, &numOfComponents);
	if ( !checkForAllocatedMemorySpace(componentsIndices, numOfComponents) ) return 0;


	size_t componentIndex, stackIndex;
	for (size_t index=0; index < numOfComponents; index++){
		componentIndex = componentsIndices[index]; stackIndex = m_indexOfStackedComponents[componentIndex];
		if ( stackIndex !=-1 ) m_stackOfComponents[stackIndex].resize(m_numOfNodes);
		std::copy(newData[index], newData[index]+numOfNodes, m_data[componentIndex]+firstDataPoint );
	}

	return 1;
}


template<class T>
int Nodes<T>::set_dataValues(std::vector<T>** newData, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
	size_t tempPointer[1]; if (!componentsIndices) componentsIndices = tempPointer;
	set_defaultInput(componentsIndices, &numOfNodes, &numOfComponents);
	if ( !checkForAllocatedMemorySpace(componentsIndices, numOfComponents) ) return 0;

	T* passedData[numOfComponents];
	for (size_t componentIndex=0; componentIndex < numOfComponents; componentIndex++) 
		passedData[componentIndex] = newData[componentIndex]->data();

	return set_dataValues(passedData, componentsIndices, firstDataPoint, numOfNodes, numOfComponents);
}


template<class T>
int Nodes<T>::set_dataValues(std::vector<T>& newData, size_t componentIndex, size_t firstDataPoint, size_t numOfNodes)
{ T* component = newData.data(); set_dataValues(&component, &componentIndex, firstDataPoint, numOfNodes, 1); }


template<class T>
int Nodes<T>::set_dataValues(std::istream& fileStream, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
	size_t tempPointer[1]; if (!componentsIndices) componentsIndices = tempPointer;
	set_defaultInput(componentsIndices, &numOfNodes, &numOfComponents);
	if ( !checkForAllocatedMemorySpace(componentsIndices, numOfComponents) ) return 0;

	size_t componentIndex, stackIndex;
	for (size_t index=0; index < numOfComponents; index++){
		componentIndex = componentsIndices[index]; stackIndex = m_indexOfStackedComponents[componentIndex];
		if ( stackIndex !=-1 ) m_stackOfComponents[stackIndex].resize(m_numOfNodes);
	}


	size_t pointIndex(firstDataPoint);
	while (	fileStream >> m_data[componentsIndices[0]][pointIndex] && pointIndex<numOfNodes-1 )
	{
		for (size_t componentIndex = 1; componentIndex < numOfComponents; componentIndex++)
			fileStream >> m_data[componentsIndices[componentIndex]][pointIndex];
		pointIndex++;
	}
	return 1;
}


template<class T>
int Nodes<T>::set_dataValues(std::string filePath, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
try 
{
	size_t tempPointer[1]; if (!componentsIndices) componentsIndices = tempPointer;
	set_defaultInput(componentsIndices, &numOfNodes, &numOfComponents);
	if ( !checkForAllocatedMemorySpace(componentsIndices, numOfComponents) ) return 0;

	size_t reqNumOfComponents(numOfComponents), reqNumOfNodes(numOfNodes);
	std::ifstream newFileStream = openFile(filePath, &reqNumOfComponents, &reqNumOfNodes);

	set_dataValues(newFileStream, componentsIndices, firstDataPoint, numOfNodes, numOfComponents);

	newFileStream.close();

	return 1;
}
catch (const std::ios_base::failure& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }
catch (const exception_badNumOfColumns& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }
catch (const exception_badNumOfLines& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }

}


// =========================================================================================


template<class T>
int Nodes<T>::insert_data(T** newData, size_t numOfNewNodes, int insertInPosition)
{
	if ( insertInPosition == -1) insertInPosition = get_numOfNodes();
	set_numOfNodes(m_numOfNodes+numOfNewNodes);
	std::vector<T*> currentData(m_data); size_t stackIndex;
	for (size_t componentIndex = 0; componentIndex < get_numOfComponents(); componentIndex++) {
		stackIndex = m_indexOfStackedComponents[componentIndex];
		if (stackIndex ==-1) {
			size_t pos = create_stackedComponent(componentIndex); m_stackOfComponents[pos].resize(m_numOfNodes);
			std::copy( currentData[componentIndex], currentData[componentIndex]+insertInPosition, m_stackOfComponents.at(pos).begin() );
			std::copy( currentData[componentIndex] +insertInPosition, currentData[componentIndex]+m_numOfNodes-numOfNewNodes, m_stackOfComponents.at(pos).begin()+insertInPosition+1+numOfNewNodes);
			std::copy(newData[componentIndex], newData[componentIndex]+numOfNewNodes, m_stackOfComponents.at(pos).begin()+insertInPosition);
		}
		else m_stackOfComponents[stackIndex].insert(m_stackOfComponents[stackIndex].begin()+insertInPosition, newData[componentIndex], newData[componentIndex]+numOfNewNodes); 
	}
}

template<class T>
int Nodes<T>::insert_data(std::vector<T>** newData, size_t numOfNewNodes, int insertInPosition)
{
	T* passedData[get_numOfComponents()];
	for (size_t componentIndex=0; componentIndex < get_numOfComponents(); componentIndex++) 
		passedData[componentIndex] = newData[componentIndex]->data();

	insert_data(passedData, numOfNewNodes, insertInPosition);
}


template<class T>
int Nodes<T>::insert_data(std::vector<T>& newData, size_t numOfNewNodes, int insertInPosition)
{
	T* passedData[1] = {newData.data()};
	insert_data(passedData, numOfNewNodes, insertInPosition);
}


template<class T>
int Nodes<T>::insert_data(std::istream& fileStream, size_t numOfNewNodes, int insertInPosition)
{

	if ( insertInPosition == -1) insertInPosition = get_numOfNodes();
	set_numOfNodes(m_numOfNodes+numOfNewNodes);

	std::vector<T*> currentData(m_data); size_t stackIndex;
	for (size_t componentIndex = 0; componentIndex < get_numOfComponents(); componentIndex++) {
		stackIndex = m_indexOfStackedComponents[componentIndex];
		if (stackIndex ==-1) {
			size_t pos = create_stackedComponent(componentIndex); m_stackOfComponents[pos].resize(m_numOfNodes);
			std::copy( currentData[componentIndex], currentData[componentIndex]+insertInPosition, m_stackOfComponents.at(pos).begin() );
			std::copy( currentData[componentIndex] +insertInPosition, currentData[componentIndex]+m_numOfNodes-numOfNewNodes, m_stackOfComponents.at(pos).begin()+insertInPosition+1+numOfNewNodes);
		}
		else m_stackOfComponents[stackIndex].insert(m_stackOfComponents[stackIndex].begin()+insertInPosition, numOfNewNodes, 0);}


	size_t pointIndex(insertInPosition);
	while (	fileStream >> m_data[0][pointIndex] && pointIndex < numOfNewNodes-1 ) {
		for (size_t componentIndex = 1; componentIndex < get_numOfComponents(); componentIndex++)
			fileStream >> m_data[componentIndex][pointIndex];
		pointIndex++;
	}
}


template<class T>
int Nodes<T>::insert_data(std::string filePath, size_t numOfNewNodes, int insertInPosition)
{
try
{
	if ( insertInPosition == -1) insertInPosition = get_numOfNodes();
	set_numOfNodes(m_numOfNodes+numOfNewNodes);

	size_t reqNumOfComponents(get_numOfComponents()), reqNumOfNodes(numOfNewNodes);
	std::ifstream fileStream = openFile(filePath, &reqNumOfComponents, &reqNumOfNodes);

	insert_data(fileStream, numOfNewNodes, insertInPosition);

	fileStream.close();
}
catch (const std::ios_base::failure& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }
catch (const exception_badNumOfColumns& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }
catch (const exception_badNumOfLines& exception) { std::cout << exception.what() << " Would you like to redefine the filepath? " ; if (0) throw; }

}



template<class T>
int Nodes<T>::remove_data(size_t fromPosition, size_t numOfNodes)
{
	std::vector<T*> currentData(m_data); size_t stackIndex;
	for (size_t componentIndex = 0; componentIndex < get_numOfComponents(); componentIndex++) {
		stackIndex = m_indexOfStackedComponents[componentIndex];
		if (stackIndex ==-1) {
			size_t pos = create_stackedComponent(componentIndex); m_stackOfComponents[pos].resize(m_numOfNodes-numOfNodes);
			std::copy(currentData[componentIndex], currentData[componentIndex]+fromPosition, m_stackOfComponents.at(pos).begin() );
			std::copy(currentData[componentIndex]+fromPosition+numOfNodes, currentData[componentIndex]+m_numOfNodes, m_stackOfComponents.at(pos).begin()+fromPosition+1);}
		else m_stackOfComponents[stackIndex].erase( m_stackOfComponents[stackIndex].begin()+fromPosition, m_stackOfComponents[stackIndex].begin()+fromPosition+numOfNodes ); }

	set_numOfNodes(m_numOfNodes-numOfNodes);
}



//===========================================================================================


template<typename T>
int Nodes<T>::get_component (std::ostream& outStream, size_t componentIndex )
{
	for (size_t pointIndex = 0; pointIndex < get_numOfNodes(); pointIndex++)
		outStream << m_data[componentIndex][pointIndex] << std::endl;
}

template<typename T>
int Nodes<T>::get_component (std::string file, size_t componentIndex )
{
	std::ofstream out; out.exceptions(std::ios_base::badbit);
	out.open(file); get_component (out); out.close();
}


template<typename T>
int Nodes<T>::get_allComponents (std::ostream& outStream)
{
	for (size_t componentIndex = 0; componentIndex < get_numOfComponents(); componentIndex++)
		get_component (outStream, componentIndex );
}


template<typename T>
int Nodes<T>::get_allComponents (std::string file )
{
	std::ofstream out; out.exceptions(std::ios_base::badbit);
	out.open(file); get_allComponents (out); out.close();
}


template<class T>
int Nodes<T>::get_components (std::ostream& outStream, size_t* componentsIndices, size_t numOfComponents)
{ get_dataValues(outStream, componentsIndices, 0, m_numOfNodes, numOfComponents); }


template<typename T>
int Nodes<T>::get_components (std::string file, size_t* componentsIndices, size_t numOfComponents)
{
	std::ofstream out; out.exceptions(std::ios_base::badbit);
	out.open(file);
	get_components ( out, componentsIndices, numOfComponents);
	out.close();
}


template<typename T>
int Nodes<T>::get_dataValues (std::ostream& outStream, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
	size_t tempPointer[1]; if (!componentsIndices) componentsIndices = tempPointer;
	set_defaultInput(componentsIndices, &numOfNodes, &numOfComponents);

	for (size_t pointIndex = firstDataPoint; pointIndex < numOfNodes; pointIndex++)
	{
		for (size_t componentIndex = 0; componentIndex < numOfComponents; componentIndex++)
			outStream << m_data[componentsIndices[componentIndex]][pointIndex] << "\t\t";
		outStream << std::endl;
	}
}


template<typename T>
int Nodes<T>::get_dataValues (std::string file, size_t* componentsIndices, size_t firstDataPoint, size_t numOfNodes, size_t numOfComponents)
{
	std::ofstream out; out.exceptions(std::ios_base::badbit);
	out.open(file);
	get_dataValues ( out, componentsIndices, firstDataPoint, numOfNodes, numOfComponents);
	out.close();
}


#endif
