#pragma once

#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include "Cluster.h"
#include "PbData.h"
#include "Util.h"

using namespace std;

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

	template<typename T>
	void initializeVariablesUCVRP(IloEnv &env, vector< unordered_map<int, T > > &u);
	template<typename T>
	void initializeObjFunctionUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u);
	template<typename T>
	void initializeCoveringCnstUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector< unordered_map<int, IloRange > > &coveringCnst);
	template<typename T>
	void determineClustersInSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u, vector<vector<int>> &sol);
	template<typename T>
	void printSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u);
	void printSolUCVRP(vector<vector<int>> &sol);
	double solveSeparationUCVRP(vector<vector<int>> &sol);
public:
	Solver(PbData *pbdata);
	~Solver();

	void enumeration();
};

