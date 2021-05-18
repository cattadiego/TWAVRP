#include "Cluster.h"



Cluster::Cluster(PbData *pbData)
{
	this->pbData = pbData;
	this->demandScenario = vector<int>(this->pbData->nbScenarios, 0);
	this->coeffObjFunction = 1;
	this->bitWords = vector<int>(this->pbData->nbWords);
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
	this->tsp.push_back(id);
	this->bitWords.at(this->pbData->wordbit.at(id)) += this->pbData->keybit.at(id);
	return atLeastFeasForOneScenario;
}
bool Cluster::isFeasibleForScenario(int s) {
	return this->demandScenario.at(s) <= this->pbData->vehicleCapacity;
}
string Cluster::print() {
	string str = "clust(";
	for (int i = 0; i < this->cluster.size(); ++i) {
		str += toString(this->cluster.at(i)) + "\t";
	}
	str += ")\ntsp(";
	for (int i = 0; i < this->tsp.size(); ++i) {
		str += toString(this->tsp.at(i)) + "\t";
	}
	str += ")\tcost: " + toString(this->cost) + "\n";
	str += "\tseq sol(";
	for (int i = 0; i < this->tspInSol.size(); ++i) {
		str += toString(this->tspInSol.at(i)) + "\t";
	}
	str += ")\tcost sol: " + toString(this->costInSol) + "\n";
	str += "\n";
	
	return str;
}
string Cluster::print(int scenario) {
	string str = "clust(";
	for (int i = 0; i < this->cluster.size(); ++i) {
		str += toString(this->cluster.at(i)) + "\t";
	}
	str += ")\ttsp(";
	for (int i = 0; i < this->cluster.size(); ++i) {
		str += toString(this->tsp.at(i)) + "\t";
	}
	str += ")\t";
	str += "quantity: " + toString(this->demandScenario.at(scenario)) + "\t";
	str += "cost: " + toString(this->cost) + "\n";

	str += "\tseq in sol(";
	for (int i = 0; i < this->tspInSol.size(); ++i) {
		str += toString(this->tspInSol.at(i)) + "\t";
	}
	str += ")\tcost sol:" + toString(this->costInSol) + "\n";
	
	return str;
}
bool Cluster::isTspFlexible() {
	switch (isTspFlex)
	{
	case -1:
		return determineIsTspFlexible();
	case 1:
		return true;
	case 0:
		return false;
	default:
		return false;
	}
}
bool Cluster::determineIsTspFlexible() {
	double t = this->pbData->lbExogenousTW.at(0);
	int prec = 0;
	for (int i = 0; i < this->tsp.size(); ++i) {
		t += this->pbData->travelTimes.at(prec).at(this->tsp.at(i));
		if (t >= this->pbData->ubExogenousTW.at(this->tsp.at(i))) {
			this->isTspFlex = 0;
			return false;
		}
		prec = this->tsp.at(i);
	}
	t += this->pbData->travelTimes.at(this->tsp.at(this->tsp.size() - 1)).at(0);
	if (t >= this->pbData->ubExogenousTW.at(0)) {
		this->isTspFlex = 0;
		return false;
	}
	this->isTspFlex = 1;
	return true;
}
bool Cluster::intersection(const Cluster &cluster) {
	int res = 0;
	for (int i = 0; i < this->pbData->nbWords; ++i)
		res += this->bitWords.at(i) & cluster.bitWords.at(i);
	return (res > 0);
}