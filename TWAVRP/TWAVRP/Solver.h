#pragma once

#include <algorithm>
#include <ilcplex/ilocplex.h>
#include "Cluster.h"
#include "PbData.h"
#include "Util.h"

class Solver
{
	PbData *pbdata;

	double permute(vector<int> a, int l, int r, double cost);
	double permuteCost(vector<int> a, int l, int r);
	double verifyCostSequence(vector<int> sequence);

	vector<Cluster> clusters;
	vector<unordered_map<int, vector<int>>> clusterFeasThatContain;
	vector<vector<int>> clusterFeasForScenario;

	void evaluateAllCluster();
	void generateAllClusters();

public:
	Solver(PbData *pbdata);
	~Solver();

	void enumeration();
};

