#include "ClusterGeneration.h"

ClusterGeneration::ClusterGeneration(PbData *pbData) {
	this->pbData = pbData;
	generateClusters();
}
void ClusterGeneration::generateClusters() {
	ProgDynTspEnum pDyn(this->pbData);
	pDyn.solve();
}
void ClusterGeneration::generateClustersDemandScenario() {
	ProgDynTspEnum pDyn(this->pbData);
	pDyn.solveWithScenarioDeamnds();
}
void ClusterGeneration::initializeData(vector<int> &demands, vector<vector<float>> &costs) {
	for (int i = 1; i < this->pbData->nbPoints; ++i) {
		int d = this->pbData->demandScenarios.at(i).at(0);
		for (int s = 1; s < this->pbData->nbScenarios; ++s) {
			if (this->pbData->demandScenarios.at(i).at(s) < d)
				d = this->pbData->demandScenarios.at(i).at(s);
		}
		demands[i] = d;
	}

	for (int i = 0; i < this->pbData->nbPoints; ++i) {
		for (int j = 0; j < this->pbData->nbPoints; ++j) {
			costs[i][j] = (float)this->pbData->travelCosts.at(i).at(j);
		}
	}
}