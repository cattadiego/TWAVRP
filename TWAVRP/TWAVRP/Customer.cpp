#include "Customer.h"

Customer::Customer(int id, double x, double y, vector<int> demandsScenario) 
	: Point(id, x, y)
{
	this->demandsScenario = demandsScenario;
}