#ifndef SCATTEREDDATAHEADER
#define SCATTEREDDATAHEADER

#include<vector>
#include<fstream>
#include"Nodes.h"

template<class T=float>
class ScatteredData
{
protected:
	Nodes<T>*		m_domain;
	Nodes<T>*		m_field;
	size_t			m_numOfNodes;
	std::vector<bool>	m_memAllocatedMembers;

public:
	ScatteredData(size_t numOfDomainDimensions=1, size_t numOfFieldsPerPoint=1, size_t numOfNodes=0);
	ScatteredData(Nodes<T>* domain, Nodes<T>* field = 0);
	~ScatteredData();


	void		set_numOfDomainDimensions	(size_t);
	void		set_numOfFieldsPerPoint		(size_t);
	void		set_numOfNodes			(size_t);


	int		set_domain 			(Nodes<T>*);
	int		set_field  			(Nodes<T>*);

	void		set_AllComponents 		();
	void		set_AllComponents 		(T** newData);
	void		set_AllComponents 		(std::string filePath);
	void		set_AllDomainComponents		();
	void		set_AllDomainComponents 	(T** newData);
	void		set_AllDomainComponents 	(std::string filePath);
	void		set_AllFieldComponents  	();
	void		set_AllFieldComponents 		(T** newData);
	void		set_AllFieldComponents  	(std::string filePath);
	void		set_domainComponent		(		size_t domainComponent);
	void		set_domainComponent		(T* newData, 	size_t domainComponent);
	void		set_fieldComponent		(		size_t fieldComponent);
	void		set_fieldComponent		(T* newData, 	size_t fieldComponent);

	int		insert_domainComponents		(T** newData, 	size_t numOfComponents);
	int		insert_fieldComponents		(T** newData, 	size_t numOfComponents);
	int		copyInsert_domainComponents 	(T** newData, size_t numOfComponents);
	int		copyInsert_fieldComponents 	(T** newData, size_t numOfComponents);


	Nodes<T>*	get_domain 			() 			{return m_domain;}
	Nodes<T>* 	get_field  			() 			{return m_field;}

	size_t		get_numOfDomainDimensions	() 			{return m_domain->get_numOfComponents();}
	size_t		get_numOfFieldsPerPoint		() 			{return m_field->get_numOfComponents();}
	size_t		get_numOfNodes			() 			{return m_numOfNodes;}
	T** 		get_AllDomainComponents 	() 			{return m_domain->get_allComponents();}
	T** 		get_AllFieldComponents  	() 			{return m_field->get_allComponents();}
	T*		get_domainComponent		(size_t componentIndex) {return m_domain->get_component(componentIndex);}
	T*		get_fieldComponent		(size_t componentIndex) {return m_field->get_component(componentIndex);}

};

template class ScatteredData<float>;
template class ScatteredData<double>;



template <typename T>
std::ostream& operator<< (std::ostream& out, ScatteredData<T>& scatteredData)
{

	for (int index = 0; index < scatteredData.get_numOfNodes(); index++)
	{
		for (int domainCompIndex = 0; domainCompIndex < scatteredData.get_numOfDomainDimensions(); domainCompIndex++) out << scatteredData.get_domainComponent(domainCompIndex)[index] << "\t";

		for (int fieldIndex = 0; fieldIndex < scatteredData.get_numOfFieldsPerPoint(); fieldIndex++) out << scatteredData.get_fieldComponent(fieldIndex)[index] << "\t";

		out << std::endl;
	}

	return out;
}


#endif
