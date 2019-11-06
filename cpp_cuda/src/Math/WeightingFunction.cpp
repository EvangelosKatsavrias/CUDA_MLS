#include"WeightingFunction.h"


std::string get_DistributionType_name( DistributionType& type )
{
	std::string name;
	switch ( type ) {
		case 0:	name = "Wendland"; break;
		case 1: name = "Gaussian"; break;
		case 2: name = "Inverse Distance"; break;
		case 3: name = "Non predefined"; break;
       	}
	return name;
}


template<class T>
std::string WeightingFunction<T>::get_distributionTypeName()
{
	DistributionType type = get_distributionType();
	return get_DistributionType_name( type );
}


std::istream& operator >> (std::istream& stream, DistributionType& var )
{
	int temp; stream >> temp; 
	switch (temp) {
		case 0: { var = Wendland; break; }
		case 1: { var = Gaussian; break; }
		case 2: { var = InverseDistance; break; }
		case 3: { var = NonPredefined; break; }
	}
	return stream;
}

std::ostream& operator << (std::ostream& stream, DistributionType& var )
{
	stream << get_DistributionType_name( var );
	return stream;
}

