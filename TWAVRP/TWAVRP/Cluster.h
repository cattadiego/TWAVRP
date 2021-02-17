#pragma once

#include <unordered_map>
#include <vector>

#include "PbData.h"

using namespace std;

class Cluster
{
	PbData *pbData;
	
public:
	Cluster(PbData *pbData);
	~Cluster();

	bool addCustomer(int id);
	bool isFeasibleForScenario(int s);
	string print();
	string print(int scenario);

	vector<int> cluster; 
	vector<int> demandScenario;
	//vector<bool> feasScenario;
	double cost;

};

