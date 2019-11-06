#ifndef CONTAINERBASEHEADERFILE
#define CONTAINERBASEHEADERFILE


#include<stdexcept>
#include<exception>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

template<class contentType, class containerType>
class ContainerBase
{
protected:
	contentType*	m_data;
	size_t		m_size;
	bool		m_internalContainerFlag;
	containerType	m_internalContainer;


public:
	ContainerBase( size_t size=0 );
	ContainerBase( size_t size, contentType* data );
	ContainerBase( ContainerBase* referenceContainerBase );
	ContainerBase( ContainerBase& referenceContainerBase );
	~ContainerBase();


	contentType*	data 	    ( );
	int		set_content ( );
	int		set_content ( contentType* data, int size=0 );
	int		set_content ( std::vector<contentType>& data );
	int		set_content ( std::istream& streamIn );
	int		set_content ( std::string fileIn );
	
	contentType*	get_content (  );
	int		get_content ( std::ostream& streamOut );
	int		get_content ( std::string fileOut );
	containerType&	get_internalContainer ( ) { return m_internalContainer; }


	typedef		contentType* iterator;
	typedef		const contentType* const_iterator;
	contentType*	begin ();
	contentType*	end ();
	contentType	front ();
	contentType	back ();
	size_t		size ();
	contentType&	operator[] (size_t index);
	int		resize(size_t newSize, contentType value = 0 );
	int		reserve(size_t newSize );
	contentType&	at ( size_t index );
	void		clear ( );
	contentType*	insert( contentType* position, contentType* begin, contentType* end );
};


template<class T, class T1>
ContainerBase<T, T1>::ContainerBase( size_t size ):
m_size(size)
{ m_internalContainerFlag = 1; m_internalContainer.reserve(size); }


template<class T, class T1>
ContainerBase<T, T1>::ContainerBase( size_t size, T* data ):
m_size(size), m_data(data)
{ m_internalContainer.clear();	m_internalContainerFlag = 0; }


template<class T, class T1>
ContainerBase<T, T1>::ContainerBase( ContainerBase<T, T1>* referenceContainerBase ):
ContainerBase( referenceContainerBase->size(), referenceContainerBase->data() ) 
{}

template<class T, class T1>
ContainerBase<T, T1>::ContainerBase( ContainerBase& referenceContainerBase ):
ContainerBase(referenceContainerBase.size(), referenceContainerBase.data() )
{ resize( referenceContainerBase.size() ); }


template<class T, class T1>
ContainerBase<T, T1>::~ContainerBase()
{ }


template<class T, class T2>
std::istream& operator >> (std::istream& streamIn, ContainerBase<T, T2>& container)
{
	for ( auto& element: container ) streamIn >> element;
	return streamIn;
}


template<class T, class T2>
std::ostream& operator << (std::ostream& streamOut, ContainerBase<T, T2>& container)
{
	for ( auto& element: container ) streamOut << element << std::endl;
	return streamOut;
}


template<class T, class T1>
T* ContainerBase<T, T1>::data ( ) { return m_data; }

template<class T, class T1>
int ContainerBase<T, T1>::set_content ( )
{ m_internalContainerFlag=1; m_internalContainer.resize(m_size, 0); m_data = m_internalContainer.data(); }

template<class T, class T1>
int ContainerBase<T, T1>::set_content ( T* data, int size ) { if (size) m_size=size; m_data = data; m_internalContainerFlag=0; m_internalContainer.clear(); }

template<class T, class T1>
int ContainerBase<T, T1>::set_content ( std::vector<T>& data ) { resize(data.size()); std::copy( data.begin(), data.end(), begin() ); }


template<class T, class T1>
int ContainerBase<T, T1>::set_content ( std::istream& streamIn ) { streamIn >> *this; }

template<class T, class T1>
int ContainerBase<T, T1>::set_content ( std::string fileIn ) {

	std::ifstream newFileStream; newFileStream.exceptions(std::ios_base::badbit);
	newFileStream.open( fileIn );

	set_content ( newFileStream );
	newFileStream.close();
}


template<class T, class T1>
T* ContainerBase<T, T1>::get_content ( ) { return m_data; }

template<class T, class T1>
int ContainerBase<T, T1>::get_content ( std::ostream& streamOut ) { streamOut << *this; }

template<class T, class T1>
int ContainerBase<T, T1>::get_content ( std::string fileOut ) { 
	
	std::ofstream newFileStream; newFileStream.exceptions(std::ios_base::badbit);
	newFileStream.open( fileOut );

	get_content ( newFileStream );
	newFileStream.close();
}

template<class T, class T1>
T* ContainerBase<T, T1>::begin (){ return m_data; }

template<class T, class T1>
T* ContainerBase<T, T1>::end (){ return m_data+m_size; }

template<class T, class T1>
T ContainerBase<T, T1>::front (){ return *m_data; }

template<class T, class T1>
T ContainerBase<T, T1>::back (){ return *(end()-1); }


template<class T, class T1>
size_t ContainerBase<T, T1>::size () { return m_size; }


template<class T, class T1>
T& ContainerBase<T, T1>::at ( size_t index ) { return m_data[index]; }


template<class T, class T1>
T& ContainerBase<T, T1>::operator[] (size_t index) { return m_data[index]; }


template<class T, class T1>
int ContainerBase<T, T1>::resize(size_t newSize, T value )
{
	m_internalContainer.resize(newSize, value); 
	if ( !m_internalContainerFlag ) 
	{ std::copy(begin(), end(), m_internalContainer.begin() );
	  m_internalContainerFlag = 1; }

	m_data = m_internalContainer.data();
	m_size = newSize;
	return 1;
}


template<class T, class T1>
int ContainerBase<T, T1>::reserve( size_t newSize )
{ 
	m_internalContainer.reserve(newSize);

	if ( !m_internalContainerFlag )
	{ std::copy(begin(), end(), m_internalContainer.begin() );
	  m_internalContainerFlag = 1; }

	m_data = m_internalContainer.data();
	m_size = newSize;
	return 1;
}


template<class T, class T1>
void ContainerBase<T, T1>::clear ( )
{
	m_internalContainer.resize( 0 ); m_size = 0;	
	m_internalContainerFlag = 1;
	m_data = m_internalContainer.data();
}

template<class T, class T1>
T* ContainerBase<T, T1>::insert( T* position, T* data_begin, T* data_end )
{
	size_t pos = position -begin();
	resize( m_size );

	m_internalContainer.insert( m_internalContainer.begin()+pos, data_begin, data_end );
	m_data = m_internalContainer.data();
	m_size = m_internalContainer.size();
}


#endif
