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

#include"fileFunctions.h"
#include<iostream>
#include<algorithm>
#include<sstream>
#include<iterator>
#include<bitset>
#include<vector>


void getFromTextFile_numOfColumnsAndLines(std::string filePath, size_t* numOfColumns, size_t* numOfLines)
{
	std::ifstream newFileStream;
	newFileStream.open(filePath);
	size_t tempCols[1], tempLines[1];

	if (numOfColumns == 0) numOfColumns = tempCols;
	if (!numOfLines)   numOfLines = tempLines;

	std::string tempLine; std::getline(newFileStream, tempLine);
	*numOfColumns = std::distance(std::istream_iterator<std::string>(std::istringstream(tempLine)>>std::ws), std::istream_iterator<std::string>());

	if ( numOfLines ) *numOfLines = std::count(std::istreambuf_iterator<char>(newFileStream), std::istreambuf_iterator<char>(), '\n') +1;


	newFileStream.close();
}



std::ifstream openFile(std::string filePath, size_t* expectedNumOfColumns, size_t* expectedNumOfLines)
{
	std::ifstream newFileStream; 
	newFileStream.exceptions(std::ios_base::badbit); 
	newFileStream.open(filePath);

	size_t numOfLines, numOfCols;
	if ( expectedNumOfLines ) getFromTextFile_numOfColumnsAndLines(filePath, &numOfCols, &numOfLines);
	else getFromTextFile_numOfColumnsAndLines(filePath, &numOfCols);

	if ( expectedNumOfLines && !(*expectedNumOfLines) ) *expectedNumOfLines = numOfLines;
	if ( expectedNumOfColumns && !(*expectedNumOfColumns) ) *expectedNumOfColumns = numOfCols;

	if ( expectedNumOfColumns && numOfCols != *expectedNumOfColumns ) throw exception_badNumOfColumns();
	if ( expectedNumOfLines && numOfLines != *expectedNumOfLines) throw exception_badNumOfLines();

	return newFileStream;
}


int get_dataLine( std::istream& streamIn, std::string& line )
{

size_t charPosition;
 
while( std::getline(streamIn, line) )
{
	charPosition = line.find_first_not_of(" \t");
	if ( charPosition != std::string::npos ) // check for empty lines
	{
		if ( line[charPosition] == '!' ||  line[charPosition] == '#' || line[charPosition] == '$' || (line[charPosition] == '/' && line[charPosition+1] == '/' ) ) continue; // the line is a comment
		else { return 1; } // data line found and stored to the string variable 'line' 
 	}
	else continue; // the line is empty
}

return 0; // end of file

}


int get_consecutiveData( std::istream& streamIn, std::string& line )
{

size_t charPosition;

//streamIn.ignore();
if ( !std::getline(streamIn, line) )	return 0; // check end of file

charPosition = line.find_first_not_of(" \t");
if ( line.empty() ) 			return 0; // check for empty lines
if ( line[charPosition] == '!' ||  line[charPosition] == '#' || line[charPosition] == '$' || (line[charPosition] == '/' && line[charPosition+1] == '/' ) ) 	return 0; // check for comment line

return 1; // data line found and stored to the string variable 'line' 

}



template<class T>
int convert_objFile( std::string fileName )
{

std::ifstream newFileStream; newFileStream.exceptions(std::ios_base::badbit); 
newFileStream.open( fileName );
std::string line; std::istringstream iss;

std::string inputType;
size_t lineOfInputIndex(0);

T x, y, z; std::vector<T> xv, yv, zv;
std::vector<T> ut, vt, wt;
std::vector<T> xvn, yvn, zvn;

int p1, p2, p3; std::vector<size_t> conn;
int p1t, p2t, p3t; std::vector<size_t> conn_t;
int p1n, p2n, p3n; std::vector<size_t> conn_n;
char tempChar;
size_t tempInt;


while( get_dataLine( newFileStream, line) ) {

	iss.clear(); iss.str(line); lineOfInputIndex++; 
	inputType.clear(); iss >> inputType;

	// read point coords
	if ( inputType.compare("v") == 0 ) {
		iss >> x >> y >> z;
		xv.push_back( x );
		yv.push_back( y );
		zv.push_back( z );
	}

	// read textures
	if ( inputType.compare("vt") == 0 ) {
		iss >> x >> y;	ut.push_back( x ); vt.push_back( y );
		if ( iss.peek() != EOF ) { iss >> z; wt.push_back( z ); }
	}

	// read normal vectors
	if ( inputType.compare("vn") == 0 ) {
		iss >> x >> y >> z;
		xvn.push_back( x ); yvn.push_back( y ); zvn.push_back( z ); }

	// read face connectivities
	if ( inputType.compare("f") == 0 ) {
		iss >> p1; conn.push_back( p1 );
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			if (iss.peek() != '/') { iss >> p1t; conn_t.push_back( p1t ); } 
		}
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			iss >> p1n; conn_n.push_back( p1n ); }

		iss >> p2; conn.push_back( p2 );
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			if (iss.peek() != '/') { iss >> p2t; conn_t.push_back( p2t ); }
		}
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			iss >> p2n; conn_n.push_back( p2n ); }

		iss >> p3; conn.push_back( p3 );
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			if (iss.peek() != '/') { iss >> p3t; conn_t.push_back( p3t ); }
		}
		if ( iss.peek() == '/' ) {
			iss >> tempChar;
			iss >> p3n; conn_n.push_back( p3n ); }
	}
}

std::vector<int> flags( xv.size(), int(0));

std::string nodesFile = fileName.substr(0, fileName.rfind('.')) + ".nod";
std::string connFile = fileName.substr(0, fileName.rfind('.')) + ".ele";
std::ofstream outFileStream; outFileStream.exceptions(std::ios_base::badbit); 


outFileStream.open(nodesFile);
outFileStream << xv.size() << std::endl;
for ( auto el: flags ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: xv ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: yv ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: zv ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: ut ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: vt ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: wt ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: xvn ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: yvn ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: zvn ) outFileStream << el << "\t";
outFileStream << std::endl;
outFileStream.close();


outFileStream.open(connFile);
outFileStream << conn.size()/3 << std::endl;
for ( auto el: conn ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: conn_t ) outFileStream << el << "\t";
outFileStream << std::endl;
for ( auto el: conn_n ) outFileStream << el << "\t";
outFileStream << std::endl;
outFileStream.close();

}

template int convert_objFile<float>( std::string fileName );
template int convert_objFile<double>( std::string fileName );

