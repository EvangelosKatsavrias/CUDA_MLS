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

#include"App_MeshAdapter.h"
#include"fileFunctions.h"


typedef double data_t;


void find_qualityStats( ContainerBase<data_t, std::vector< data_t > >* qualityMetrics, data_t* min, data_t* max, data_t* mean, data_t* stdDev ) {

	data_t meanQualityValue(0);
	for ( auto& value: qualityMetrics[0] ) {
		meanQualityValue += value;
	}
	auto minmaxQuality = std::minmax_element( qualityMetrics->begin(), qualityMetrics->end() );
	meanQualityValue /= qualityMetrics->size();

	data_t qualityStdDeviation(0);
	for ( auto& value: qualityMetrics[0] ) {
		qualityStdDeviation += (value-meanQualityValue)*(value-meanQualityValue);
	}
	qualityStdDeviation /= qualityMetrics->size()-1;
	qualityStdDeviation = sqrt(qualityStdDeviation);

	*min = *minmaxQuality.first;
       	*max = *minmaxQuality.second;
	*mean = meanQualityValue;
	*stdDev = qualityStdDeviation;
}


void parametricTest()
{

size_t numOfDeg(5), numOfExp(7), numOfSpans(10);
size_t degVariants[numOfDeg] = { 0, 1, 2, 3, 4 };
size_t expVariants[numOfExp] = { 0, 1, 2, 3, 4, 5, 6 };
size_t spanVariants[numOfSpans] = { 6, 8, 10, 12, 15, 18, 20, 25, 30, 40 };

clock_t time = -clock();
// ==================== 2D example =============================

{
app_meshAdapter<data_t> app1;
//Mesh2DDisplacedBoundaryAdapter<double, int, int, 3>* name = new Mesh2DDisplacedBoundaryAdapter<double, int, int, 3>;

app1.readProjectSettings( "Leaf/Leaf_30.msh" );
//app1.createProject_console();

app1.set_initialMesh( );
app1.set_adaptedMesh( );
app1.set_meshAdapter( );
app1.write_nodesInfo ( std::cout );
//app1.adaptMesh();


ContainerBase<data_t, std::vector<data_t> > minQualities; minQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > maxQualities; maxQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > meanQualities; meanQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > stdDevQualities; stdDevQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > badElements; badElements.resize( numOfDeg*numOfExp*numOfSpans );


std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfDeg*numOfExp*numOfSpans); 
size_t caseIndex(0);
for ( size_t deg = 0; deg < numOfDeg; deg++ ) {
	app1.get_meshAdapter()->set_polynomialDegree( degVariants[ deg ] );
	for ( size_t exp = 0; exp < numOfExp; exp++ ) {
		app1.get_meshAdapter()->set_interpolationFactor( expVariants[ exp ] );
		for ( size_t span = 0; span < numOfSpans; span++ )
		{
			caseIndex = deg*numOfExp*numOfSpans +exp*numOfSpans +span;

			app1.get_meshAdapter()->set_weightingFunctionSpan( spanVariants[ span ] );
			app1.set_adaptedMesh();
			app1.adaptMesh();

			find_qualityStats( app1.get_adaptedMesh()->get_qualityMetrics(), minQualities.data()+caseIndex, maxQualities.data()+caseIndex, meanQualities.data()+caseIndex, stdDevQualities.data()+caseIndex );

			badElements[ caseIndex ] = app1.get_meshAdapter()->get_badElements().size();
			std::cout << "\r" << int(caseIndex*nodesFraction) << "% "; std::cout.flush();
		}
	}
}
time += clock();

std::cout << "\r" << "100.00%\n\n";

minQualities.get_content( "minQualities_30.dat" );
maxQualities.get_content( "maxQualities_30.dat"  );
meanQualities.get_content( "meanQualities_30.dat" );
stdDevQualities.get_content( "stdDevQualities_30.dat" );
badElements.get_content( "badElements_30.dat" );
}



time -= clock();
{
app_meshAdapter<data_t> app1;
app1.readProjectSettings( "Leaf/Leaf_45.msh" );

app1.set_initialMesh( );
app1.set_adaptedMesh( );
app1.set_meshAdapter( );
app1.write_nodesInfo ( std::cout );


ContainerBase<data_t, std::vector<data_t> > minQualities; minQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > maxQualities; maxQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > meanQualities; meanQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > stdDevQualities; stdDevQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > badElements; badElements.resize( numOfDeg*numOfExp*numOfSpans );


std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfDeg*numOfExp*numOfSpans); 
size_t caseIndex(0);
for ( size_t deg = 0; deg < numOfDeg; deg++ ) {
	app1.get_meshAdapter()->set_polynomialDegree( degVariants[ deg ] );
	for ( size_t exp = 0; exp < numOfExp; exp++ ) {
		app1.get_meshAdapter()->set_interpolationFactor( expVariants[ exp ] );
		for ( size_t span = 0; span < numOfSpans; span++ )
		{
			caseIndex = deg*numOfExp*numOfSpans +exp*numOfSpans +span;

			app1.get_meshAdapter()->set_weightingFunctionSpan( spanVariants[ span ] );
			app1.set_adaptedMesh();
			app1.adaptMesh();

			find_qualityStats( app1.get_adaptedMesh()->get_qualityMetrics(), minQualities.data()+caseIndex, maxQualities.data()+caseIndex, meanQualities.data()+caseIndex, stdDevQualities.data()+caseIndex );
			badElements[ caseIndex ] = app1.get_meshAdapter()->get_badElements().size();

			std::cout << "\r" << int(caseIndex*nodesFraction) << "% "; std::cout.flush();
		}
	}
}

time += clock();
std::cout << "\r" << "100.00%\n\n";

minQualities.get_content( "minQualities_45.dat" );
maxQualities.get_content( "maxQualities_45.dat"  );
meanQualities.get_content( "meanQualities_45.dat" );
stdDevQualities.get_content( "stdDevQualities_45.dat" );
badElements.get_content( "badElements_45.dat" );
}



time -= clock();

{
app_meshAdapter<data_t> app1;
app1.readProjectSettings( "Leaf/Leaf_60.msh" );

app1.set_initialMesh( );
app1.set_adaptedMesh( );
app1.set_meshAdapter( );
app1.write_nodesInfo ( std::cout );


ContainerBase<data_t, std::vector<data_t> > minQualities; minQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > maxQualities; maxQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > meanQualities; meanQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > stdDevQualities; stdDevQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > badElements; badElements.resize( numOfDeg*numOfExp*numOfSpans );


std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfDeg*numOfExp*numOfSpans); 
size_t caseIndex(0);
for ( size_t deg = 0; deg < numOfDeg; deg++ ) {
	app1.get_meshAdapter()->set_polynomialDegree( degVariants[ deg ] );
	for ( size_t exp = 0; exp < numOfExp; exp++ ) {
		app1.get_meshAdapter()->set_interpolationFactor( expVariants[ exp ] );
		for ( size_t span = 0; span < numOfSpans; span++ )
		{
			caseIndex = deg*numOfExp*numOfSpans +exp*numOfSpans +span;

			app1.get_meshAdapter()->set_weightingFunctionSpan( spanVariants[ span ] );
			app1.set_adaptedMesh();
			app1.adaptMesh();

			find_qualityStats( app1.get_adaptedMesh()->get_qualityMetrics(), minQualities.data()+caseIndex, maxQualities.data()+caseIndex, meanQualities.data()+caseIndex, stdDevQualities.data()+caseIndex );
			badElements[ caseIndex ] = app1.get_meshAdapter()->get_badElements().size();

			std::cout << "\r" << int(caseIndex*nodesFraction) << "% "; std::cout.flush();
		}
	}
}

time += clock();
std::cout << "\r" << "100.00%\n\n";

minQualities.get_content( "minQualities_60.dat" );
maxQualities.get_content( "maxQualities_60.dat"  );
meanQualities.get_content( "meanQualities_60.dat" );
stdDevQualities.get_content( "stdDevQualities_60.dat" );
badElements.get_content( "badElements_60.dat" );
}




time -= clock();


{
app_meshAdapter<data_t> app1;
app1.readProjectSettings( "Leaf/Leaf_80.msh" );

app1.set_initialMesh( );
app1.set_adaptedMesh( );
app1.set_meshAdapter( );
app1.write_nodesInfo ( std::cout );



ContainerBase<data_t, std::vector<data_t> > minQualities; minQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > maxQualities; maxQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > meanQualities; meanQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > stdDevQualities; stdDevQualities.resize( numOfDeg*numOfExp*numOfSpans );
ContainerBase<data_t, std::vector<data_t> > badElements; badElements.resize( numOfDeg*numOfExp*numOfSpans );


std::cout << "\nProgress: " << std::endl << "0%"; float nodesFraction = float(100)/((float) numOfDeg*numOfExp*numOfSpans); 
size_t caseIndex(0);
for ( size_t deg = 0; deg < numOfDeg; deg++ ) {
	app1.get_meshAdapter()->set_polynomialDegree( degVariants[ deg ] );
	for ( size_t exp = 0; exp < numOfExp; exp++ ) {
		app1.get_meshAdapter()->set_interpolationFactor( expVariants[ exp ] );
		for ( size_t span = 0; span < numOfSpans; span++ )
		{
			caseIndex = deg*numOfExp*numOfSpans +exp*numOfSpans +span;

			app1.get_meshAdapter()->set_weightingFunctionSpan( spanVariants[ span ] );
			app1.set_adaptedMesh();
			app1.adaptMesh();

			find_qualityStats( app1.get_adaptedMesh()->get_qualityMetrics(), minQualities.data()+caseIndex, maxQualities.data()+caseIndex, meanQualities.data()+caseIndex, stdDevQualities.data()+caseIndex );
			badElements[ caseIndex ] = app1.get_meshAdapter()->get_badElements().size();

			std::cout << "\r" << int(caseIndex*nodesFraction) << "% "; std::cout.flush();
		}
	}
}

time += clock();
std::cout << "\r" << "100.00%\n\n";

minQualities.get_content( "minQualities_80.dat" );
maxQualities.get_content( "maxQualities_80.dat"  );
meanQualities.get_content( "meanQualities_80.dat" );
stdDevQualities.get_content( "stdDevQualities_80.dat" );
badElements.get_content( "badElements_80.dat" );
}



float time_s = (float)time/CLOCKS_PER_SEC;


std::cout << "The parametric analysis has been successfully completed, with execution time: " << time_s << std::endl;

}

