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

#include"MeshAdapters.h"


std::istream& operator >> (std::istream& stream, displacementsSource& var )
{
	int temp; stream >> temp; 
	var = static_cast<displacementsSource> ( temp );
	return stream;
}


std::istream& operator >> (std::istream& stream, solutionMethod& var )
{
	int temp; stream >> temp; 
	var = static_cast<solutionMethod> ( temp );
	return stream;
}


std::string get_displacementsSource_name( displacementsSource& source )
{
	std::string name;
	switch ( source ) {
		case 0:
			name = "ProvidedInFile"; break;
		case 1:
			name = "AffineTransformation"; break;
		case 2:
			name = "Translation"; break;
		case 3:
			name = "Rotation_EulerAngles"; break;
		case 4:
			name = "Rotation_AxisAngle"; break;
		case 5:
			name = "Rotation_Quaternion"; break;
	}
	return name;
}

std::string get_solutionMethod_name( solutionMethod& method )
{
	std::string name;
	switch ( method ) {
		case 0:
			name = "MLS_LCS"; break;
	}
	return name;
}


std::ostream& operator << (std::ostream& stream, displacementsSource& var )
{
	stream << get_displacementsSource_name( var );
	return stream;
}


std::ostream& operator << (std::ostream& stream, solutionMethod& var )
{
	stream << get_solutionMethod_name( var );
	return stream;
}

