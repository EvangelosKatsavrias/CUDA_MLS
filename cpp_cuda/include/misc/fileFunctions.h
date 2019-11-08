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

#ifndef FILEFUNCTIONSHEADERFILE
#define FILEFUNCTIONSHEADERFILE

#include<stdexcept>
#include<iostream>
#include<fstream>
#include<sstream>
#include<errno.h>
#include<string.h>


void getFromTextFile_numOfColumnsAndLines(std::string filePath, size_t* numOfColumns, size_t* numOfLines = 0);

std::ifstream openFile(std::string filePath, size_t* expectedNumOfColumns = 0, size_t* expectedNumOfLines = 0);

class exception_badNumOfColumns: std::exception
{ public: const char* what() const noexcept {return "Bad number of columns in the provided file. ";} };

class exception_badNumOfLines: std::exception
{ public: const char* what() const noexcept {return "Bad number of lines in the provided file. ";} };


int get_dataLine( std::istream& streamIn, std::string& line );
int get_consecutiveData( std::istream& streamIn, std::string& line );

template<class T>
int convert_objFile( std::string fileName );

#endif
