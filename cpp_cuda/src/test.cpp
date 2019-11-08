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

void clock_test();
void myBLAStest();
void topologyData_test();
void classicLeastSquaresTest();
void weightedLSTest();
void gridDataTest();
int  meshViewerTest( int argc, char *argv[] );
int app3Dtest();
int app3Dtest2( int argc, char *argv[] );
void parametricTest();

int main(  int argc, char *argv[] )
{

//clock_test();
//myBLAStest();
//topologyData_test();
//classicLeastSquaresTest();
//weightedLSTest();
gridDataTest();
//parametricTest();
//meshViewerTest(  argc, argv );
//app3Dtest();
//app3Dtest2( argc, argv );

}
