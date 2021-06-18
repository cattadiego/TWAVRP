#pragma once

#include <bitset>
#include <unordered_map>
#include <vector>

#include "PbData.h"

using namespace std;

class Cluster
{
	PbData *pbData;
	int isTspFlex = -1;
	bool determineIsTspFlexible();

public:
	Cluster(PbData *pbData);
	~Cluster();

	bool addCustomer(int id);
	bool isFeasibleForScenario(int s);
	string print();
	string print(int scenario);

	vector<int> cluster;
	vector<int> tsp;
	vector<int> tspInSol;
	vector<int> demandScenario;
	vector<int> bitWords;

	bool intersection(const Cluster &cluster);
	
	double cost;
	double costInSol;

	int coeffObjFunction;		// to be used in separation

	bool isTspFlexible();

};

