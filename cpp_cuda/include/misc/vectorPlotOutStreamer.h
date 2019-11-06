#ifndef VECTORPLOTERHEADER
#define VECTORPLOTERHEADER

#include<vector>
#include<set>
#include<iostream>

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& givenVector)
{
	for (auto& element: givenVector) out << element << std::endl;

	return out;
}


//#include<iterator>
// if ( !v.empty() ) std::copy (v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));


template <typename T>
std::istream& operator>> (std::istream& streamIn, std::vector<T>& givenVector)
{
	for (auto& element: givenVector) streamIn >> element;

	return streamIn;
}


template <typename T>
std::ostream& operator<< (std::ostream& out, const std::set<T>& givenVector)
{
	for (auto& element: givenVector) out << element << std::endl;

	return out;
}


template <typename T>
std::istream& operator>> (std::istream& streamIn, std::set<T>& givenVector)
{
	for (auto& element: givenVector) streamIn >> element;

	return streamIn;
}



#endif
