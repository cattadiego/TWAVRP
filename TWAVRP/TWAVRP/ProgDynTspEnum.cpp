#include "ProgDynTspEnum.h"

ProgDynTspEnum::ProgDynTspEnum() {
	this->test();	
}
ProgDynTspEnum::ProgDynTspEnum(PbData *pbData) {
	this->exactcounter = 0;
	this->pbData = pbData;
	this->size = this->pbData->nbPoints;
	this->Q = this->pbData->vehicleCapacity;
	this->initializeData();
	this->initializeVectors();
	this->initializeKeys();
}
ProgDynTspEnum::ProgDynTspEnum(vector<int> const &dem, vector<vector<float>> const &cost, int Q) {
	this->demand = dem;
	this->cost = cost;
	this->Q = Q;
	this->size = this->demand.size();
	this->initializeVectors();
	this->initializeKeys();	
}
bool ProgDynTspEnum::checkHeuristicDominance(int iState, int currj) {
	int i4 = lastv[iState];
	if (i4 == 0) return true;
	int i3 = lastv[pred[iState]];
	if (i3 == 0) return true;
	int i2 = lastv[pred[pred[iState]]];
	if (i2 == 0) return true;

	float delta = cost[i2][i4] + cost[i4][i3] + cost[i3][currj];		//new sequence
	float deltaMinus = cost[i2][i3] + cost[i3][i4] + cost[i4][currj];	//old sequence

	if (delta - deltaMinus < -myEpsilon) return false;

	int i1 = lastv[pred[pred[pred[iState]]]];
	if (i1 == 0) return true;

	delta = cost[i1][i3] + cost[i3][i4] + cost[i4][i2] + cost[i2][currj];
	deltaMinus += cost[i1][i2];
	if (delta - deltaMinus < -myEpsilon) return false;

	delta = cost[i1][i4] + cost[i4][i2] + cost[i2][i3] + cost[i3][currj];
	if (delta - deltaMinus < -myEpsilon) return false;

	delta = cost[i1][i4] + cost[i4][i3] + cost[i3][i2] + cost[i2][currj];
	if (delta - deltaMinus < -myEpsilon) return false;

	return true;
}
bool ProgDynTspEnum::checkExactDominance(int iState, int currj, int q) {
	//int nextState = head[card[iState] + 1][stateDemand[iState] + this->demand[currj]][currj];
	int nextState = head[card[iState] + 1][q][currj][sumId[iState] + currj];
	int currk = set[iState] + this->key[currj];
	float currcost = stateCost[iState] + this->cost[lastv[iState]][currj];
	while (nextState != -1) {
		if (currk == set[nextState]) {
			//here the current state either dominates or is dominated
			//in the first case the old state is replaced, in the second case the new is not kept
			//in both cases the number of states remains the same
			if (stateCost[nextState] > currcost) {
				stateCost[nextState] = currcost;
				pred[nextState] = iState;
			}
			return false;
		}
		nextState = next[nextState];
	}
	//in this case the new state must be memorized
	return true;
}
string ProgDynTspEnum::getSetFromKey(int key) {
	string str = "{";
	for (int i = 0; i < this->key.size(); ++i) {
		if ((key & this->key[i]) != 0) {
			str += toString(i) + "\t";
		}
	}
	str += "}";
	return str;
}
void ProgDynTspEnum::initializeKeys() {
	this->key = vector<int>(this->size);
	for (int j = 0; j < this->size; ++j) {
		this->key[j] = pow(2, j);
	}
}
void ProgDynTspEnum::initializeVectors() {

	int nSize = (pow(2, (size - 1)) + 1);
	//int nSize = 2000000;
	

	this->lastv = vector<int>(nSize);
	this->stateDemand = vector<int>(nSize);
	this->stateCost = vector<float>(nSize);
	this->stateArr = vector<float>(nSize);
	this->set = vector<int>(nSize);
	this->pred = vector<int>(nSize);
	this->sumId = vector<int>(nSize);

	this->card = vector<int>(nSize);
	this->next = vector<int>(nSize);
	
	/*this->lastv = new int(nSize);
	this->stateDemand = new int(nSize);
	this->stateCost = new float(nSize);
	this->stateArr = new float(nSize);
	this->set = new int(nSize);
	this->pred = new int(nSize);;

	this->card = new int(nSize);
	this->next = new int(nSize);*/

	int sumMaxId = this->pbData->nbPoints * (this->pbData->nbPoints + 1) / 2;
	this->head = vector<vector<vector<vector<int>>>>(size, vector<vector<vector<int>>>(Q + 1, vector<vector<int>>(size, vector<int>(sumMaxId, -1))));

	this->lastv[0] = 0;
	this->stateDemand[0] = 0;
	this->stateCost[0] = 0;
	this->stateArr[0] = 0;
	this->set[0] = 0;
	this->pred[0] = -1;
	this->sumId[0] = 0;

	this->card[0] = 0;
	this->next[0] = 0;
	//this->head[0][0][0] = -1;

}
void ProgDynTspEnum::initializeData() {
	this->demand = vector<int>(this->pbData->nbPoints);
	for (int i = 1; i < this->pbData->nbPoints; ++i) {
		int d = this->pbData->demandScenarios.at(i).at(0);
		for (int s = 1; s < this->pbData->nbScenarios; ++s) {
			if (this->pbData->demandScenarios.at(i).at(s) < d)
				d = this->pbData->demandScenarios.at(i).at(s);
		}
		this->demand[i] = d;
	}

	this->cost = vector<vector<float>>(this->pbData->nbPoints, vector<float>(this->pbData->nbPoints));
	this->travelTime = vector<vector<float>>(this->pbData->nbPoints, vector<float>(this->pbData->nbPoints));
	for (int i = 0; i < this->pbData->nbPoints; ++i) {
		for (int j = 0; j < this->pbData->nbPoints; ++j) {
			this->cost[i][j] = (float)this->pbData->travelCosts.at(i).at(j);
			this->travelTime[i][j] = (float)this->pbData->travelTimes.at(i).at(j);
		}
	}
}
string ProgDynTspEnum::printAllStates() {
	string str = "";
	for (int iState = 0; iState < lastv.size(); ++iState) {
		str += printState(iState);
	}
	return str;
}
string ProgDynTspEnum::printState(int iState) {
	string str = "";
	str += "iState: " + toString(iState) + "\n";
	str += getSetFromKey(set[iState]) + "\n";
	str += "cost: " + toString(stateCost[iState]) + "\n";
	str += "demand: " + toString(stateDemand[iState]) + "\n";
	str += "pred: " + toString(pred[iState]) + "\n";
	str += "lastv: " + toString(lastv[iState]) + "\n";
	str += "sumId: " + toString(sumId[iState]) + "\n";
	str += "-------\n";
	return str;
}
void ProgDynTspEnum::resize() {
	int nSize = 2 * next.size();
	lastv.resize(nSize);
	stateDemand.resize(nSize);
	stateCost.resize(nSize);
	stateArr.resize(nSize);
	set.resize(nSize);
	pred.resize(nSize);
	sumId.resize(nSize);
}
void ProgDynTspEnum::resizeExactDominance() {
	resize();
	int nSize = 2 * next.size();
	next.resize(nSize);
	card.resize(nSize);
}
void ProgDynTspEnum::solve() {

	int nState = 0;

	//set[0] = 0;
	//pred[0] = -1;

	for (int iState = 0; iState <= nState; iState++) {
		//cout << iState << endl;
#ifdef CONSOLE_VERBOSE
		printGreenMessage(printState(iState, lastv, demand, cost, set, pred));
#endif // CONSOLE_VERBOSE

		if ((set[iState] & this->key[0]) == 0) {	// the set of the state should not contain the depot
			for (int j = 0; j < size; ++j) {
#ifdef CONSOLE_VERBOSE
				cout << "considering adding " + toString(j) << endl;
#endif // CONSOLE_VERBOSE
				int q = this->stateDemand[iState] + this->demand[j];
				//if (checkCapacity(iState, j)) {
				if (q <= this->Q) {
					if ((set[iState] & this->key[j]) == 0) {
						// checkExactDominance(iState, j, q)&&
						if (checkHeuristicDominance(iState, j)) {

							/*************************************/
							bool check = true;
							int nextState = head[card[iState] + 1][q][j][sumId[iState] + j];
							int currk = set[iState] + this->key[j];
							float currcost = stateCost[iState] + this->cost[lastv[iState]][j];
							while (nextState != -1) {
								if (currk == set[nextState]) {
									//here the current state either dominates or is dominated
									//in the first case the old state is replaced, in the second case the new is not kept
									//in both cases the number of states remains the same
									if (stateCost[nextState] > currcost) {
										stateCost[nextState] = currcost;
										pred[nextState] = iState;
									}
									check = false;
									break;
								}
								nextState = next[nextState];
							}
							//in this case the new state must be memorized
							/*************************************/
							if (check) {
								++nState;
								if (nState > this->set.size() - 1) resizeExactDominance();
								updateState(iState, nState, j, q);
								updateStateExactDominance(iState, nState);
							}
													
//#ifdef CONSOLE_VERBOSE
//							else printRedMessage("exactly dominated: eliminated");
//#endif // CONSOLE_VERBOSE
						}
//#ifdef CONSOLE_VERBOSE
//						else printRedMessage("heuristically dominated: eliminated");
//#endif // CONSOLE_VERBOSE
					}
//#ifdef CONSOLE_VERBOSE
//					else printRedMessage("iState already contains " + toString(j) + " eliminated");
//#endif // CONSOLE_VERBOSE
				}
#ifdef CONSOLE_VERBOSE
				else printRedMessage("capacity violated - eliminated");
#endif // CONSOLE_VERBOSE
			}
		}
#ifdef CONSOLE_VERBOSE
		else printRedMessage("contains zero - not extended");
#endif // CONSOLE_VERBOSE
	}
	ofstream stream("progDyn.dat", ios::app);
	stream << exactcounter << "\t" << nState << "\t";
	stream.close();
}
void ProgDynTspEnum::solveWithScenarioDeamnds() {

	//scenarioDemands = vector<unordered_map<int, int>*>(next.size());

	int nState = 0;

	set[0] = 0;
	pred[0] = -1;
	scenarioDemands[0] = new unordered_map<int, int>();
	for (int s = 0; s < this->pbData->nbScenarios; ++s)
		scenarioDemands[0]->insert(make_pair(s, 0));

	for (int iState = 0; iState <= nState; iState++) {
		//cout << iState << endl;
#ifdef CONSOLE_VERBOSE
		printGreenMessage(printState(iState, lastv, demand, cost, set, pred));
#endif // CONSOLE_VERBOSE

		if ((set[iState] & this->key[0]) == 0) {	// the set of the state should not contain the depot
			for (int j = 0; j < size; ++j) {
#ifdef CONSOLE_VERBOSE
				cout << "considering adding " + toString(j) << endl;
#endif // CONSOLE_VERBOSE
				if ((set[iState] & this->key[j]) == 0) {
					//cout << (set[iState] & this->key[j]) << endl;
					unordered_map<int, int> *demandNState = new unordered_map<int, int>();
					if (checkCapacityScenarioDemand(iState, j, demandNState)) {
						if (!checkHeuristicDominance(iState, j)) {
							//if (!checkExactDominance(iState, j)) {
								//cout << getSetFromKey(set[iState] + this->key[j]) << endl;
							++nState;
							/*if (nState > this->set.size() - 1) {
								resizeExactDominance();
								scenarioDemands.resize(next.size());
							}*/
							updateState(iState, nState, j, 0);
							updateStateExactDominance(iState, nState);
							scenarioDemands[nState] = demandNState;
							//}
#ifdef CONSOLE_VERBOSE
							else printRedMessage("exactly dominated: eliminated");
#endif // CONSOLE_VERBOSE
						}
#ifdef CONSOLE_VERBOSE
						else printRedMessage("heuristically dominated: eliminated");
#endif // CONSOLE_VERBOSE
					}
#ifdef CONSOLE_VERBOSE
					else printRedMessage("iState already contains " + toString(j) + " eliminated");
#endif // CONSOLE_VERBOSE
				}
#ifdef CONSOLE_VERBOSE
				else printRedMessage("capacity violated - eliminated");
#endif // CONSOLE_VERBOSE
			}
		}
#ifdef CONSOLE_VERBOSE
		else printRedMessage("contains zero - not extended");
#endif // CONSOLE_VERBOSE
	}
	ofstream stream("progDyn.dat", ios::app);
	stream << nState << "\t";
	stream.close();
}
void ProgDynTspEnum::test() {
	demand.push_back(0);
	demand.push_back(3);
	demand.push_back(4);
	demand.push_back(3);
	demand.push_back(2);

	cost.push_back({ 0, 3, 4, 2, 1 });
	cost.push_back({ 3, 0, 5, 2, 1 });
	cost.push_back({ 4, 5, 0, 3, 4 });
	cost.push_back({ 2, 2, 3, 0, 5 });
	cost.push_back({ 1, 1, 4, 5, 0 });
	
	size = 5;
	Q = 20;

	this->initializeKeys();

	this->solve();
	
}
void ProgDynTspEnum::updateState(int iState, int nState, int currj, int quantity) {
	lastv[nState] = currj;
	//stateDemand[nState] = stateDemand[iState] + demand[currj];
	stateDemand[nState] = quantity;
	stateCost[nState] = stateCost[iState] + this->cost[lastv[iState]][currj];
	set[nState] = (set[iState] + this->key[currj]);
	pred[nState] = iState;
	sumId[nState] = sumId[iState] + currj;
}
void ProgDynTspEnum::updateStateExactDominance(int iState, int nState) {
	card[nState] = card[iState] + 1;
	sumId[nState] = sumId[iState] + lastv[nState];
	next[nState] = head[card[nState]][stateDemand[nState]][lastv[nState]][sumId[nState]];
	head[card[nState]][stateDemand[nState]][lastv[nState]][sumId[nState]] = nState;
}
bool ProgDynTspEnum::checkArrTime(int iState, int currj) {
	return (this->stateArr[iState] + this->travelTime[lastv[iState]][currj] <= this->pbData->ubExogenousTW.at(currj));
}
bool ProgDynTspEnum::checkCapacity(int iState, int currj) {
	return (stateDemand[iState] + demand[(currj)] <= Q);
}
bool ProgDynTspEnum::checkCapacityScenarioDemand(int iState, int currj, unordered_map<int, int> *demandNState) {
	if (currj == 0) {
		for (auto it = this->scenarioDemands[iState]->begin(); it != this->scenarioDemands[iState]->end(); ++it) {
			demandNState->insert(make_pair(it->first, it->second));			
		}
		return true;
	}
	for (auto it = this->scenarioDemands[iState]->begin(); it != this->scenarioDemands[iState]->end(); ++it) {
		if (it->second + this->pbData->demandScenarios[currj][it->first] <= Q) {
			demandNState->insert(make_pair(it->first, it->second + this->pbData->demandScenarios[currj][it->first]));
		}
	}
	return (demandNState->size() > 0);
}