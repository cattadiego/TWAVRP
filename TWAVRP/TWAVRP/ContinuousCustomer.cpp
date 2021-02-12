#include "ContinuousCustomer.h"

ContinuousCustomer::ContinuousCustomer(int id, double x, double y, vector<int> demandsScenario, int lbExogTW, int ubExogTW, int twWidth) 
	: Customer(id, x, y, demandsScenario)
{
	this->lbExogTW = lbExogTW;
	this->ubExogTW = ubExogTW;
	this->twWidth = twWidth;
}