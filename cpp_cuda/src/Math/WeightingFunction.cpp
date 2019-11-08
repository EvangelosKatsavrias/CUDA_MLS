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

