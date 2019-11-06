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
