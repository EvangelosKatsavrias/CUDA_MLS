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

#ifndef MATRIXMEMORYSTORAGEHEADER
#define MATRIXMEMORYSTORAGEHEADER

#include<vector>
#include"ContainerBase.h"
#include"LUDecomposers.h"


template<class T>
int find_numOfBands(T* A, int n);

int map2BandedStorage(int i, int j, int n);
int mapFromBandedStorage(int i_band, int n, int *i, int *j);

template<class T>
int find_skylineColumns(T* A, int n, int* n_max);

int map2SkylineStorage();
int mapFromSkylineStorage();


template<class T>
void convert_matrixStorage_2Banded(T* A, int n, int numOfBands, T* A_band);

template<class T>
void convert_matrixStorage_2Skyline(T* A, int n, int* n_max, T* A_skyline);


template<class T>
class Matrix: public ContainerBase<T, std::vector<T> >
{
	T m_numOfRows;
	T m_numOfCols;

	Matrix() {}
	~Matrix() {}

};


template<class T>
class SquareMatrix: Matrix<T>
{

	SquareMatrix() {}
	~SquareMatrix() {}

};


template<class T>
class BandedMatrix: public ContainerBase<T, std::vector<T> >
{

typedef ContainerBase<T, std::vector<T>> baseCont_t ;

protected:
	int		m_numOfColumns;
	int		m_numOfBands;
	T*		m_referenceFullMatrix;


public:
	BandedMatrix ( int numOfColumns=3, int numOfBands=1, T* bandedMatrix=0 ): baseCont_t(numOfColumns*numOfBands), m_numOfColumns(numOfColumns), m_numOfBands(numOfBands) { if (bandedMatrix) this->set_content( bandedMatrix ); }
	BandedMatrix ( T* referenceFullMatrix, int numOfColumns ): m_numOfColumns(numOfColumns)
	{ set_fromFullMatrix(referenceFullMatrix, m_numOfColumns); }


	int set_numOfBands( int numOfBands ) { m_numOfBands = numOfBands; this->resize(m_numOfColumns*m_numOfBands); }
	int get_numOfBands() { return m_numOfBands; }

	int set_numOfColumns( int numOfColumns ) { m_numOfColumns = numOfColumns; this->resize(m_numOfColumns*m_numOfBands);}
	int get_numOfColumns() { return m_numOfColumns; }

	T* get_bandedMatrix() { return this->data(); }
	int set_bandedMatrix( T* bandedMatrix, int numOfBands=0, int numOfColumns=0 ) { if (numOfBands) m_numOfBands =numOfBands; if (numOfColumns ) m_numOfColumns= numOfColumns; this->set_content(bandedMatrix, numOfBands*numOfColumns); }
	int set_fromFullMatrix( T* fullMatrix, int numOfColumns =0, int numOfBands=0 ) 
	{ 
		m_referenceFullMatrix = fullMatrix; 
		if ( numOfColumns ) m_numOfColumns = numOfColumns;
		if ( numOfBands ) set_numOfBands(numOfBands);
		else set_numOfBands( find_NumOfBands( fullMatrix, m_numOfColumns ) );
		convert_fullMatrix( fullMatrix, m_numOfColumns, m_numOfBands, this->data() ); 
	}


	static int find_NumOfBands ( T* fullMatrix, int numOfColumns ) { return find_numOfBands(fullMatrix, numOfColumns);  }
	static void convert_fullMatrix(T* fullMatrix, int numOfColumns, int numOfBands, T* bandedMatrix) { convert_matrixStorage_2Banded(fullMatrix, numOfColumns, numOfBands, bandedMatrix); }

};
template class BandedMatrix <float>;
template class BandedMatrix <double>;


template<class T>
class SkylineMatrix: public ContainerBase<T, std::vector<T> >
{

typedef ContainerBase<T, std::vector<T> > baseCont_t ;

protected:
	int					m_numOfColumns;
	ContainerBase<int, std::vector<int> >	m_columnPointers;
	T*					m_referenceFullMatrix;


public:
	SkylineMatrix ( int numOfColumns=3, int* columnPointers=0, T* skylineMatrix=0 ): baseCont_t(0, skylineMatrix), m_numOfColumns(numOfColumns), m_columnPointers(numOfColumns+1) { m_columnPointers.set_content(); if (columnPointers) set_columnPointers(columnPointers);  if (skylineMatrix) this->set_content(skylineMatrix, m_columnPointers.back()); }

	SkylineMatrix ( T* referenceFullMatrix, int numOfColumns ): m_numOfColumns(numOfColumns)
	{ set_fromFullMatrix(referenceFullMatrix, numOfColumns ); }


	int set_numOfColumns( int numOfColumns ) { m_numOfColumns = numOfColumns; set_columnPointers(); this->resize(0); }
	int get_numOfColumns( ) { return m_numOfColumns; }


	int* get_columnPointers( ) { return m_columnPointers.data(); }
	int set_columnPointers( ) { m_columnPointers.resize(m_numOfColumns+1); }
	int set_columnPointers( int* columnPointers, int numOfColumns=0 ) { if (numOfColumns) set_numOfColumns(numOfColumns); m_columnPointers.set_content(columnPointers); this->resize(m_columnPointers[m_numOfColumns]); }


	T* get_skylineMatrix() { return this->data(); }
	int set_skylineMatrix( T* skylineMatrix, int* columnPointers=0, int numOfColumns=0 ) { if (columnPointers) set_columnPointers(columnPointers, numOfColumns); this->set_content(skylineMatrix); }

	int set_fromFullMatrix( T* fullMatrix, int numOfColumns=0 )
	{
		m_referenceFullMatrix = fullMatrix; if( numOfColumns ) set_numOfColumns(numOfColumns); 
		find_SkylineColumns( fullMatrix, m_numOfColumns, m_columnPointers.data() );
		this->resize(m_columnPointers.back()); 
		convert_fullMatrix( fullMatrix, m_numOfColumns, m_columnPointers.data(), this->data() );
	}


	static int find_SkylineColumns ( T* fullMatrix, int numOfColumns, int* columnPointers ) { return find_skylineColumns(fullMatrix, numOfColumns, columnPointers); }
	static void convert_fullMatrix(T* fullMatrix, int numOfColumns, int* columnPointers, T* skylineMatrix) { convert_matrixStorage_2Skyline(fullMatrix, numOfColumns, columnPointers, skylineMatrix ); }

};
template class SkylineMatrix <float>;
template class SkylineMatrix <double>;


#endif
