#pragma once

#include <algorithm>

#include "Cluster.h"
#include "PbData.h"

class Solver
{
	PbData *pbdata;

	double permute(vector<int> a, int l, int r, double cost);
	double permuteCost(vector<int> a, int l, int r);
	double verifyCostSequence(vector<int> sequence);

public:
	Solver(PbData *pbdata);
	~Solver();

	void generateAllClusters();
};

