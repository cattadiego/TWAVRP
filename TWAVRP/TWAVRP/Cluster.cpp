#include "Cluster.h"



Cluster::Cluster(PbData *pbData)
{
	this->pbData = pbData;
	this->demandScenario = vector<int>(this->pbData->nbScenarios, 0);
//	this->feasScenario = vector<bool>(this->pbData->nbScenarios, true);
}
Cluster::~Cluster()
{
}
bool Cluster::addCustomer(int id) {
	bool atLeastFeasForOneScenario = false;
	for (int s = 0; s < this->pbData->nbScenarios; ++s) {
		if (this->demandScenario.at(s) + this->pbData->demandScenarios.at(id).at(s) <= this->pbData->vehicleCapacity) 
			atLeastFeasForOneScenario = true;
		this->demandScenario.at(s) += this->pbData->demandScenarios.at(id).at(s);
	}
	this->cluster.push_back(id);
	return atLeastFeasForOneScenario;
}
bool Cluster::isFeasibleForScenario(int s) {
	return this->demandScenario.at(s) <= this->pbData->vehicleCapacity;
}
string Cluster::print() {
	string str = "";
	for (int i = 0; i < this->cluster.size(); ++i) {
		str += toString(this->cluster.at(i)) + "\t";
	}
	str += "\n";
	str += "cost: " + toString(this->cost) + "\n";
	return str;
}
string Cluster::print(int scenario) {
	string str = "";
	for (int i = 0; i < this->cluster.size(); ++i) {
		str += toString(this->cluster.at(i)) + "\t";
	}
	str += "\n";
	str += "quantity: " + toString(this->demandScenario.at(scenario)) + "\t";
	str += "cost: " + toString(this->cost) + "\n";
	return str;
}