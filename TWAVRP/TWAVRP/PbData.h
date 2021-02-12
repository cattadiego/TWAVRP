#pragma once

#include <string>
#include <vector>
using namespace std;



class PbData {

public:

	int vehicleCapacity;
	vector<double> scenarioProbability;

	string instanceName;

	void readData();

	PbData(string instanceName);

};