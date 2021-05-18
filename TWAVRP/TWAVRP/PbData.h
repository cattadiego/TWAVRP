#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "Config.h"

using namespace std;

class PbData {

	void readData();
	void readContinuousTW(ifstream &stream);
	void readDiscreteTW(ifstream &stream);
	void readDistancesCosts(ifstream &stream);
	void calculateBigM();
	Config config;
public:

	int vehicleCapacity;
	int nbPoints;
	int nbCustomers;
	int nbScenarios;
	vector<int> id;
	vector<int> idCustomers;
	vector<int> keybit;
	vector<int> wordbit;
	int nbWords;

	vector<double> scenarioProbability;
	unordered_map<int, vector<int>> demandScenarios;
	unordered_map<int, double> x;
	unordered_map<int, double> y;

	/* continuous customer */ 
	unordered_map<int, int> lbExogenousTW;
	unordered_map<int, int> ubExogenousTW;
	unordered_map<int, int> widthEndogenousTW;

	/* discrete customer */
	unordered_map<int, vector<int>> lbTWs;
	unordered_map<int, vector<int>> ubTWs;

	unordered_map <int, unordered_map<int, double> > travelTimes;
	unordered_map <int, unordered_map<int, double> > travelCosts;
	string instanceName;

	double bigM;

	PbData(string instanceName, Config config);
	PbData(string instanceName, Config config, int nbScenarios);
	string printTravelCost(int i, int j);
	string printTravelTime(int i, int j);
};