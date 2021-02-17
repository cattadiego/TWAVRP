#pragma once

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

template<typename T>
string toString(T nt) {
	ostringstream strs;
	strs << nt;
	return strs.str();
}
