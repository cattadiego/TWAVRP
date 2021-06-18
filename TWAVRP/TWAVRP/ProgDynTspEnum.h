#pragma once

#include <vector>
#include <iostream>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "MyParam.h"
#include "PbData.h"
#include "Util.h"

using namespace std;

class ProgDynTspEnum {

private:

	int size;
	int Q;
	PbData *pbData;
	int exactcounter;

	//int* lastv;
	//int*  stateDemand;
	//float*  stateCost;
	//float*  stateArr;
	//int*  set;
	//int*  pred;

	//int*  card;
	//int*  next;

	vector<int> lastv;
	vector<int>  stateDemand;
	vector<float>  stateCost;
	vector<float>  stateArr;
	vector<int>  set;
	vector<int>  pred;
	vector<int>  sumId;

	vector<int>  card;
	vector<int>  next;

	vector<vector<vector<vector<int>>>> head;


	vector<int> demand, key;
	vector<vector<float>> cost;
	vector<vector<float>> travelTime;

	vector<unordered_map<int, int>*> scenarioDemands;

	string getSetFromKey(int key);

	bool checkHeuristicDominance(int iState, int currj);
	bool checkExactDominance(int iState, int currj, int q);			//0 not dominated nor dominate, 1 dominated, 2 dominate another
	bool checkArrTime(int iState, int currj);
	bool checkCapacity(int iState, int currj);
	bool checkCapacityScenarioDemand(int iState, int currj, unordered_map<int, int> *demandNState);
	void initializeData();
	void initializeKeys();
	void initializeVectors();
	string printState(int iState);
	string printAllStates();
	void resize();
	void resizeExactDominance();
	void updateState(int iState, int nState, int currj, int quantity);
	void updateStateExactDominance(int iState, int nState);

public:
	ProgDynTspEnum();
	ProgDynTspEnum(PbData *pbData);
	ProgDynTspEnum(vector<int> const &dem, vector<vector<float>> const &cost, int Q);
	void solve();
	void solveWithScenarioDeamnds();
	void test();

};