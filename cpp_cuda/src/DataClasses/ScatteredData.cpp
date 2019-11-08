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

#include"ScatteredData.h"
#include<stdexcept>
#include<iostream>


template<class T>
ScatteredData<T>::ScatteredData(size_t numOfDomainDimensions, size_t numOfFieldsPerPoint, size_t numOfNodes): m_numOfNodes(numOfNodes), m_memAllocatedMembers(2, 1)
{
	m_domain = new Nodes<T>(numOfDomainDimensions, numOfNodes);
	m_field  = new Nodes<T>(numOfFieldsPerPoint, numOfNodes);
}


template<class T>
ScatteredData<T>::ScatteredData(Nodes<T>* domain, Nodes<T>* field): m_domain(domain), m_field(field), m_numOfNodes(domain->get_numOfNodes()), m_memAllocatedMembers(2, 0)
{ if(!field) { new Nodes<T>(1, domain->get_numOfNodes()); m_memAllocatedMembers[1] = 1;} }


template<class T>
ScatteredData<T>::~ScatteredData()
{
	if(m_memAllocatedMembers[0] == 1) delete(m_domain);
	if(m_memAllocatedMembers[1] == 1) delete(m_field);
}


template<class T>
void ScatteredData<T>::set_numOfDomainDimensions(size_t domainDimensions)
{ m_domain->set_numOfComponents(domainDimensions); }


template<class T>
void ScatteredData<T>::set_numOfFieldsPerPoint(size_t numOfFields)
{ m_field->set_numOfComponents(numOfFields); }


template<class T>
void ScatteredData<T>::set_numOfNodes(size_t numOfNodes)
{
	m_domain->set_numOfNodes(numOfNodes);
	m_field->set_numOfNodes(numOfNodes);
	m_numOfNodes = numOfNodes;
}

template<class T>
void ScatteredData<T>::set_AllComponents()
{ m_domain->set_components(); m_field->set_components();}


template<class T>
void ScatteredData<T>::set_AllComponents(T** newData)
{ m_domain->set_components(newData); m_field->set_components(newData+get_numOfDomainDimensions()); }

#include"vectorPlotOutStreamer.h"
template<class T>
void ScatteredData<T>::set_AllComponents(std::string filePath)
{
	std::ifstream newFileStream; newFileStream.open(filePath);

	std::vector<T> *tempDomain = new std::vector<T>[get_numOfDomainDimensions()];
	std::vector<T> *tempFields = new std::vector<T>[get_numOfFieldsPerPoint()];

	for (size_t index=0; index<get_numOfDomainDimensions(); index++)
		tempDomain[index].resize(get_numOfNodes());

	for (size_t index=0; index<get_numOfFieldsPerPoint(); index++)
		tempFields[index].resize(get_numOfNodes());

	size_t pointIndex(0);
	while (	newFileStream >> tempDomain[0][pointIndex] )
	{
		for (size_t index=1;index<get_numOfDomainDimensions();index++) 
			newFileStream >> tempDomain[index][pointIndex];
		for (size_t index=0;index<get_numOfFieldsPerPoint();index++) 
			newFileStream >> tempFields[index][pointIndex];
		pointIndex++;
	}

	newFileStream.close();

	m_domain->set_components();
	m_field ->set_components();
	m_domain->set_dataValues(&tempDomain);
	m_field ->set_dataValues(&tempFields);

	delete[] tempDomain;
	delete[] tempFields;
}

template<class T>
void ScatteredData<T>::set_AllDomainComponents()
{ m_domain->set_components(); }


template<class T>
void ScatteredData<T>::set_AllDomainComponents(T** newData)
{ m_domain->set_components(newData); }


template<class T>
void ScatteredData<T>::set_AllDomainComponents(std::string filePath)
{ m_domain->set_components(filePath); }

template<class T>
void ScatteredData<T>::set_AllFieldComponents()
{ m_field->set_components(); }


template<class T>
void ScatteredData<T>::set_AllFieldComponents(T** newData)
{ m_field->set_components(newData); }

template<class T>
void ScatteredData<T>::set_domainComponent(size_t domainComponent)
{ m_domain->set_components(&domainComponent, 1); }


template<class T>
void ScatteredData<T>::set_domainComponent(T* newData, size_t domainComponent)
{ m_domain->set_components(&newData, &domainComponent, 1); }


template<class T>
void ScatteredData<T>::set_fieldComponent(size_t fieldComponent)
{ m_field->set_components(&fieldComponent, 1); }


template<class T>
void ScatteredData<T>::set_fieldComponent(T* newData, size_t fieldComponent)
{ m_field->set_components(&newData, &fieldComponent, 1); }


template<class T>
int ScatteredData<T>::set_domain(Nodes<T>* newData)
{ m_domain = newData; }

template<class T>
int ScatteredData<T>::set_field(Nodes<T>* newData)
{ m_field = newData; }


template<class T>
int ScatteredData<T>::insert_domainComponents(T** newData, size_t numOfComponents)
{ return m_domain->insert_components( newData, numOfComponents); }

template<class T>
int ScatteredData<T>::insert_fieldComponents(T** newData, size_t numOfComponents)
{ return m_field->insert_components( newData, numOfComponents); }

template<class T>
int ScatteredData<T>::copyInsert_domainComponents(T** newData, size_t numOfComponents)
{ return m_domain->copyInsert_components( newData, numOfComponents); }

template<class T>
int ScatteredData<T>::copyInsert_fieldComponents(T** newData, size_t numOfComponents)
{ return m_field->copyInsert_components( newData, numOfComponents); }


