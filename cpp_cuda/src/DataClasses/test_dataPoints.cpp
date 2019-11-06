#include"Nodes.h"
#include<iostream>
#include<stdexcept>
#include"vectorPlotOutStreamer.h"
#include"StationaryPoints.h"


using namespace std;

void Nodes_test()
{

using seltype = float;

std::vector<seltype> v1(10, 5);
std::vector<seltype> v2(10, 2);
std::vector<seltype> v3(10, 1);

std::cout << v1 << std::endl;
std::cout << v2 << std::endl;
std::cout << v3 << std::endl;

seltype* data[3];
data[0] = v1.data(); data[1]=v2.data(); data[2] = v3.data();


Nodes<seltype> mydata(3, 10, data);

seltype *c1 = mydata[0], *c2 = mydata[1], *c3 = mydata[2];

v1.assign(10, 10);

cout << "number of components " << mydata.get_numOfComponents() << endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c1[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c2[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c3[index] << std::endl;


size_t compind[1]; compind[0] = 0; seltype* tt[1]; tt[0] = v2.data();
mydata.set_dataValues(tt, compind, 0, 10, 1);
for (int index=0; index<10; index++) std::cout << index << ". " << c1[index] << std::endl;
std::cout << v1 << std::endl;


std::vector<seltype>* dd[3]; dd[0] = &v3; dd[1] = &v1; dd[2] = &v2;
size_t compind2[3] = {0, 1, 2}; 
mydata.set_dataValues(dd, compind2, 0, 10, 3);
for (int index=0; index<10; index++) std::cout << index << ". " << c1[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c2[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c3[index] << std::endl;

v1.assign(10, 5);
mydata.set_dataValues(v1, 1, 0, 10);
for (int index=0; index<10; index++) std::cout << index << ". " << c1[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c2[index] << std::endl;

std::string dataSetPath = "../../exe/dataset_1/data2";
size_t compind3[2] = {0, 2}; 
mydata.set_dataValues(dataSetPath, compind3, 0, 10, 2);
for (int index=0; index<10; index++) std::cout << index << ". " << c1[index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << c3[index] << std::endl;

mydata.insert_data(data, 10, int(10));
for (int index=0; index<20; index++) std::cout << index << ". " << mydata[0][index] << std::endl;
for (int index=0; index<20; index++) std::cout << index << ". " << mydata[1][index] << std::endl;
for (int index=0; index<20; index++) std::cout << index << ". " << mydata[2][index] << std::endl;


mydata.insert_data(dd, 10, int(20));
for (int index=0; index<30; index++) std::cout << index << ". " << mydata[0][index] << std::endl;
for (int index=0; index<30; index++) std::cout << index << ". " << mydata[1][index] << std::endl;
for (int index=0; index<30; index++) std::cout << index << ". " << mydata[2][index] << std::endl;


dataSetPath = "../../exe/dataset_1/data3";
mydata.insert_data(dataSetPath, 10, int(10));
for (int index=0; index<40; index++) std::cout << index << ". " << mydata[0][index] << std::endl;
for (int index=0; index<40; index++) std::cout << index << ". " << mydata[1][index] << std::endl;
for (int index=0; index<40; index++) std::cout << index << ". " << mydata[2][index] << std::endl;


mydata.remove_data(10, 30);
for (int index=0; index<10; index++) std::cout << index << ". " << mydata[0][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata[1][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata[2][index] << std::endl;


Nodes<double> mydata2(dataSetPath);
cout << "num of components "<< mydata2.get_numOfComponents() << endl;
cout << "num of points "<< mydata2.get_numOfNodes() << endl;

for (int index=0; index<10; index++) std::cout << index << ". " << mydata2[0][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata2[1][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata2[2][index] << std::endl;


Nodes<seltype> mydata3(3, dd);
cout << "num of components "<< mydata3.get_numOfComponents() << endl;
cout << "num of points "<< mydata3.get_numOfNodes() << endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata3[0][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata3[1][index] << std::endl;
for (int index=0; index<10; index++) std::cout << index << ". " << mydata3[2][index] << std::endl;


Nodes<float> dp1;
StationaryPoints<float> sp1(&dp1);

}
