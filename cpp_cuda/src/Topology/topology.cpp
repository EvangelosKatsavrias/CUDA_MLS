#include"Topology.h"
#include"TopologyAdapters.h"

int main(){

Connectivities<size_t, 3> conn;
Topology<double, 3> topo1(2, 3, 0, &conn), topo2;

TopologyAdapter<double> topoAdapt1( &topo1, &topo2);

}
