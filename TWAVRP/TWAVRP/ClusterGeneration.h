#pragma once

#include "PbData.h"
#include "ProgDynTspEnum.h"

class ClusterGeneration {
public: 

	ClusterGeneration(PbData *pbData);

private:

	PbData *pbData;

	void generateClusters();
	void generateClustersDemandScenario();
	void initializeData(vector<int> &demands, vector<vector<float>> &costs);
};
