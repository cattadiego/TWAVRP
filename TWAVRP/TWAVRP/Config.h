#pragma once

#include <string>
#include "Util.h"

using namespace std;




class Config {

public:

	Config(user user, instanceType instanceType);
	Config(Config &config);
	Config();

	instanceType instance;
	string filePath;
	string instanceFolder;
};

