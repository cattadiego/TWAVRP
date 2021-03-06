#include "PbData.h"

PbData::PbData(string instanceName, Config config) {
	this->instanceName = instanceName;
	this->config = config;
	readData();
	calculateBigM();
}

PbData::PbData(string instanceName, Config config, int nbScenarios) {
	this->instanceName = instanceName;
	this->config = config;
	this->nbScenarios = nbScenarios;
	readData();
	calculateBigM();
}

void PbData::readData() {
	
	string comment, file;

	if(this->instanceName.find("ExtraScenarios") == string::npos)
		file = this->config.filePath + this->config.instanceFolder + this->instanceName + ".txt";
	else {
		int pos = this->instanceName.find("ExtraScenarios");
		string instName = this->instanceName;
		file = this->config.filePath + this->config.instanceFolder + instName.erase(pos, 14) + ".txt";
	}
	ifstream stream; stream.open(file.c_str());

	if (!stream)
	{
		cout << "Error opening the file" << endl;
		cout << "Path: " << this->config.filePath << endl;
		cout << "Name: " << this->instanceName << endl;
		system("pause");
		exit(1);
	}

	stream >> comment;
	stream >> this->nbCustomers;
	this->nbPoints = this->nbCustomers + 1;

	stream >> comment;
	stream >> this->vehicleCapacity;

	stream >> comment;
	if (this->instanceName.find("ExtraScenarios") == string::npos)
		stream >> this->nbScenarios;

	
	do
	{
		stream >> comment;
	} while (comment != "Location_Coordinates:");
	

	int nbBytes = sizeof(int);
	int nbBits = pow(2, nbBytes);
	this->nbWords = ceil((double)nbPoints / (double)(nbBits - 2));

	int i = 0;
	do {
		++i;
		int locId;
		double locX, locY;
		stream >> locId;
		stream >> locX;
		stream >> locY;
		this->id.push_back(locId);

		int word = floor(locId / (double)(nbBits - 2));
		this->wordbit.push_back(word);
		this->keybit.push_back(pow(2, locId - (nbBits - 2) * word));
		if(locId != 0) this->idCustomers.push_back(locId);
		this->x.insert(make_pair(locId, locX));
		this->x.insert(make_pair(locId, locY));
	} while (i < nbPoints);

	int nbOriginalScenarios = this->nbScenarios;
	if (this->instanceName.find("ExtraScenarios") != string::npos) {
		stream.close();
		string instFolder = this->config.instanceFolder;
		instFolder.erase(instFolder.length() - 2, 2);
		instFolder += "ExtraScenarios//";
		file = this->config.filePath + instFolder + this->instanceName + ".txt";
		stream.open(file.c_str());
		do
		{
			stream >> comment;
		} while (comment != "Number_Of_Scenarios:");
		stream >> nbOriginalScenarios;
		this->nbScenarios = min(this->nbScenarios, nbOriginalScenarios);
	}

	

	do
	{
		stream >> comment;
	} while (comment != "Demand_per_scenario:");

	for (int i = 0; i < nbCustomers; ++i) {
		int locId;
		stream >> locId;
		vector< int > d;
		for (int s = 0; s < this->nbScenarios; ++s) {			
			int locD; stream >> locD;
			d.push_back(locD);
		}
		for (int s = this->nbScenarios; s < nbOriginalScenarios; ++s) {
			int locD; stream >> locD;
		}
		this->demandScenarios.insert(make_pair(locId, d));
	}

	do
	{
		stream >> comment;
	} while (comment != "Probability_of_each_scenario_occuring:");

	if (this->nbScenarios == nbOriginalScenarios) {
		for (int s = 0; s < this->nbScenarios; ++s) {
			double ps;
			stream >> ps;
			this->scenarioProbability.push_back(ps);
		}
	}
	else {
		for (int s = 0; s < this->nbScenarios; ++s) {
			double ps = 1.0 / this->nbScenarios;
			this->scenarioProbability.push_back(ps);
		}
	}

	if (this->instanceName.find("ExtraScenarios") != string::npos) {
		stream.close();
		int pos = this->instanceName.find("ExtraScenarios");
		string instName = this->instanceName;
		file = this->config.filePath + this->config.instanceFolder + instName.erase(pos, 14) + ".txt";
		stream.open(file.c_str());
	}

	switch (this->config.instance) {
	case instanceType::twa:
		readContinuousTW(stream);
		break;
	case instanceType::dtwa:
		readDiscreteTW(stream);
		break;
	}
	readDistancesCosts(stream);
}
void PbData::readContinuousTW(ifstream &stream) {
	string comment;
	do
	{
		stream >> comment;
	} while (comment != "Exogenous_time_windows:");
	for (int i = 0; i < this->nbPoints; ++i) {
		stream >> comment;
		int valTW;
		stream >> valTW;
		this->lbExogenousTW.insert(make_pair(i, valTW));
		stream >> comment;
		stream >> valTW;
		this->ubExogenousTW.insert(make_pair(i, valTW));
		stream >> comment;
	}

	do
	{
		stream >> comment;
	} while (comment != "Endogenous_time_window_width:");
	for (int i = 1; i < this->nbPoints; ++i) {
		int valTW;
		stream >> valTW;
		this->widthEndogenousTW.insert(make_pair(i, valTW));
	}
}
void PbData::readDiscreteTW(ifstream &stream) {
	string comment;
	do
	{
		stream >> comment;
	} while (comment != "NumberOfWindows-TimeWindows:");
	for (int i = 0; i < this->nbPoints; ++i) {
		int nbTW;
		stream >> nbTW;
		vector<int> lb, ub;
		for (int j = 0; j < nbTW; ++j) {
			stream >> comment;
			int valTW;
			stream >> valTW;
			lb.push_back(valTW);
			stream >> comment;
			stream >> valTW;
			ub.push_back(valTW);
			stream >> comment;
		}
		this->lbTWs.insert(make_pair(i, lb));
		this->ubTWs.insert(make_pair(i, ub));
	}
}
void PbData::readDistancesCosts(ifstream &stream) {
	string comment;
	do
	{
		stream >> comment;
	} while (comment != "Travel_Costs:");

	for (int i = 0; i < this->nbPoints; ++i) {		
		for (int j = 0; j < this->nbPoints; ++j) {
			double val; 
			stream >> val;
			try {
				this->travelCosts.at(i).insert(make_pair(j, val));
			}
			catch (...) {
				this->travelCosts.insert(make_pair(i, unordered_map<int, double>()));
				this->travelCosts.at(i).insert(make_pair(j, val));
			}
		}
	}

	do
	{
		stream >> comment;
	} while (comment != "Travel_Time:");

	for (int i = 0; i < this->nbPoints; ++i) {
		for (int j = 0; j < this->nbPoints; ++j) {
			double val;
			stream >> val;
			try {
				this->travelTimes.at(i).insert(make_pair(j, val));
			}
			catch (...) {
				this->travelTimes.insert(make_pair(i, unordered_map<int, double>()));
				this->travelTimes.at(i).insert(make_pair(j, val));
			}
		}
	}

}
void PbData::calculateBigM() {
	this->bigM = 0;
	for (int i = 0; i < this->nbPoints; ++i) {
		for (int j = 0; j < this->nbPoints; ++j) {
			this->bigM += this->travelCosts.at(i).at(j);
		}
	}
}
string PbData::printTravelCost(int i, int j) {
	if (i < 0 || i > this->nbPoints) return "index out of bounds";
	if (j < 0 || j > this->nbPoints) return "index out of bounds";
	return "t(" + toString(i) + ", " + toString(j) + "):" + toString(this->travelCosts.at(i).at(j));
}
string PbData::printTravelTime(int i, int j) {
	if (i < 0 || i > this->nbPoints) return "index out of bounds";
	if (j < 0 || j > this->nbPoints) return "index out of bounds";
	return "t(" + toString(i) + ", " + toString(j) + "):" + toString(this->travelTimes.at(i).at(j));
}