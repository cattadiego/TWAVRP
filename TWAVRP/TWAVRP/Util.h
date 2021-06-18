#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <windows.h>  

#define myEpsilon 0.0001
#define myCutOff 0.1

#define INSTANCE_SIZE 10

using namespace std;

enum instanceType {
	twa = 1,
	dtwa = 2
};

enum user {
	diego = 1,
	lab = 2,
	home = 3
};

template<typename T>
string toString(T nt) {
	ostringstream strs;
	strs << nt;
	return strs.str();
}

void subsetsUtil(vector<int>& A, vector<vector<int> >& res, vector<int>& subset, int index);
vector<vector<int>> subsets(vector<int>& A);

void printGreenMessage(string message);
void printGreenMessage(string message, string header);
void printRedMessage(string message);
void printRedMessage(string message, string header);

void writeInStream(string streamName, string message);