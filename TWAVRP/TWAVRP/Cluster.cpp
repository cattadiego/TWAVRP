#include "Cluster.h"



Cluster::Cluster(PbData *pbData)
{
	this->pbData = pbData;
	this->demandScenario = vector<int>(this->pbData->nbScenarios, 0);
}
Cluster::~Cluster()
{
}
bool Cluster::addCustomer(int id) {
	for (int s = 0; s < this->pbData->nbScenarios; ++s) {
		if (this->demandScenario.at(s) + this->pbData->demandScenarios.at(id).at(s) > this->pbData->vehicleCapacity)
			return false;
		this->demandScenario.at(s) += this->pbData->demandScenarios.at(id).at(s);
	}
	this->cluster.push_back(id);
	return true;
}