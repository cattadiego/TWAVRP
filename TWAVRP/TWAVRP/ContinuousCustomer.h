#pragma once

#include "Customer.h"

class ContinuousCustomer : Customer{
public:

	int lbExogTW; 
	int ubExogTW;
	int twWidth;

	ContinuousCustomer(int id, double x, double y, vector<int> demandsScenario, int lbExogTW, int ubExogTW, int twWidth);
};