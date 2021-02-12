#pragma once

#include <string>

using namespace std;

#define myEpsilon 0.00000001
#define myCutOff 0.1

enum instance {
	twa = 1,
	dtwa = 2
};

enum user {
	diego = 1,
	lab = 2
};

class Config {

public:

	Config(user user, instance instance);

	string filePath;
	string instanceFolder;
};