#ifndef STATIONARYPOINTS
#define STATIONARYPOINTS

#include"Nodes.h"

template<class T=float>
class StationaryPoints: public Nodes<T>
{
public:
	T* get_spatialData(int component){ return this->get_component(component);}
	using Nodes<T>::Nodes;
};


template class StationaryPoints<float>;
template class StationaryPoints<double>;


#endif
