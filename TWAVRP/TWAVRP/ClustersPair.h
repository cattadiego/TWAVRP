#pragma once

#include <bitset>
#include "Cluster.h"

using namespace std;

class ClustersPair {

public:
	Cluster *c1;
	Cluster *c2;

	int pos1;	// memorizes the index in the vector of clusters of cluster c1
	int pos2;	// memorizes the index in the vector of clusters of cluster c2

	ClustersPair(Cluster *c1, Cluster *c2, int pos1, int pos2);

	bool intersection(ClustersPair &cp);
};
