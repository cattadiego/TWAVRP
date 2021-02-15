#pragma once

#include <unordered_map>
#include <vector>

#include "PbData.h"

using namespace std;

class Cluster
{
	PbData *pbData;
	vector<int> demandScenario;

public:
	Cluster(PbData *pbData);
	~Cluster();

	bool addCustomer(int id);

	vector<int> cluster;
	double cost;

};

