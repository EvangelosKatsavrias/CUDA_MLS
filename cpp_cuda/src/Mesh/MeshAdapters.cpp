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

