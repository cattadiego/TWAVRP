#pragma once

#include <iostream>
#include <string>
#include <vector>

#define myEpsilon 0.00000001
#define myCutOff 0.1

using namespace std;

enum instanceType {
	twa = 1,
	dtwa = 2
};

enum user {
	diego = 1,
	lab = 2
};

void subsetsUtil(vector<int>& A, vector<vector<int> >& res, vector<int>& subset, int index);
vector<vector<int>> subsets(vector<int>& A);