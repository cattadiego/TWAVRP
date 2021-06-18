#pragma once

#include <algorithm>
#include <ilcplex/ilocplex.h>
#include <iostream>
#include <map>
#include <windows.h>  
#include <set>  
#include "Cluster.h"
#include "ClustersPair.h"
#include "MyParam.h"
#include "MyStructures.h"
#include "PbData.h"
#include "Util.h"

using namespace std;

class Solver
{
	PbData *pbdata;

	double permute(vector<int> a, vector<int> &aBest, int l, int r, double cost, double &bestCost);
	double permuteCost(vector<int> a, int l, int r);
	double verifyCostSequence(vector<int> sequence);

	vector<Cluster> clusters;
	vector<unordered_map<int, vector<int>>> clusterFeasThatContain;
	vector<vector<int>> clusterFeasForScenario;

	int calculateNbClusters();

	void evaluateAllCluster();
	void generateAllCluster();
	void fillClusterStructures();
	void generateEvaluateAllClusters();

	template<typename T>
	void initializeVariablesUCVRP(IloEnv &env, vector< unordered_map<int, T > > &u);
	template<typename T>
	void initializeVariablesUCVRPCutSymmetries(IloEnv &env, vector<T> &y, vector< unordered_map<int, T > > &z);
	void initializeVariablesUCVRPValidCuts(IloEnv &env, unordered_map<int, unordered_map<int, IloNumVar > > &deltaSS);
	template<typename T>
	void initializeObjFunctionUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u);
	template<typename T>
	void initializeObjFunctionUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, IloNumVar &delta);
	template<typename T>
	void initializeCoveringCnstUCVRP(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector< unordered_map<int, IloRange > > &coveringCnst);
	template<typename T>
	void initializeLinkUandYCnstUCVRPCutSymmetries(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector<IloIntVar> &y, vector< unordered_map<int, IloRange > > &linkUandYCnst);
	template<typename T>
	void determineClustersInSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u, vector<vector<int>> &sol);
	template<typename T>
	void determineClustersInSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u, unordered_map<int, vector<int>> &sol);
	template<typename T>
	void printSolUCVRP(IloCplex &cplex, vector< unordered_map<int, T > > &u);
	string printSolUCVRP(vector<vector<int>> &sol);
	string printSolUCVRP(unordered_map<int, vector<int>> &sol);
	double solveSeparationUCVRP(vector<vector<int>> &sol);
	double solveSeparationUCVRP(unordered_map<int, vector<int>> &sol);
	template<typename T>
	int addFeasibilityCuts(IloEnv &env, IloModel &model, vector< unordered_map<int, T > > &u, vector< unordered_map<int, T > > &z, unordered_map<int, vector<int>> &invSol);
	int addInfeasibilityCuts(IloEnv &env, IloModel &model, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol);
	int addValidCutsOnY(IloEnv &env, IloModel &model, IloNumVar &delta, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter);
	int addValidCutsOnY(IloEnv &env, IloModel &model, IloNumVar &delta, unordered_map<int, unordered_map<int, IloNumVar>> &deltaCC, unordered_map<int, unordered_map<int, IloNumVar>> &gammaCC, vector< unordered_map<int, IloIntVar > > &u, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter, UnorderedMapKeyPairOrderedInt& memoryCuts);
	int addValidCutsOnU(IloEnv &env, IloModel &model, IloNumVar &delta, unordered_map<int, unordered_map<int, IloNumVar > > &deltaSS, vector< unordered_map<int, IloIntVar > > &u, unordered_map<int, vector<int>> &sol);
	int addValidCuts(IloEnv &env, IloModel &model, IloNumVar &delta, vector< unordered_map<int, IloIntVar > > &u, vector<IloIntVar> &y, unordered_map<int, vector<int>> &sol, int iter, UnorderedMapKeyPairOrderedInt& memoryCuts);
	/* from scenario as a key and vector of clusters associated with each scenario, 
	we obtain a solution where the key is the cluster id, and the value is the vector of scenarios in which it is present */
	void invertSolution(unordered_map<int, vector<int>> &sol, unordered_map<int, vector<int>> &invSol);


	void determineDisjointPairsOfClusters(IloEnv &env, IloModel &model, double tspCost, double twaCost, unordered_map<int, unordered_map<int, IloNumVar>> &deltaCC, vector<IloIntVar> &y, vector<vector<ClustersPair>> &setpairs, vector<int> &s, vector<int> s1, int c, int c1);
	void determineMultipleOccurencesOfPairs(IloEnv &env, IloModel &model, IloNumVar &delta, double tspCost, double twaCost, unordered_map<int, unordered_map<int, IloNumVar>> &gammaCC, vector< unordered_map<int, IloIntVar > > &u, vector<int> &s, vector<int> s1, int c, int c1);

	void insertClusterPairInSets(ClustersPair &cp, vector<vector<ClustersPair>> &setpairs);

public:
	Solver(PbData *pbdata);
	~Solver();

	void enumeration();
	void enumerationCutSymmetries();
	void test();
};

