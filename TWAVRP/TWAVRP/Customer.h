#pragma once

#include <vector>

#include "Point.h"

using namespace std;

class Customer : Point{
public:

	vector<int> demandsScenario;

	Customer(int id, double x, double y, vector<int> demandsScenario);
};