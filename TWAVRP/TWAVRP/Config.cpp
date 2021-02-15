#include "Config.h"

Config::Config() {

}
Config::Config(Config &config) {
	this->filePath = config.filePath;
	this->instance = config.instance;
	this->instanceFolder = config.instanceFolder;
}

Config::Config(user user, instanceType instanceType) {

	this->instance = instanceType;
	switch (user)
	{
	case diego:
		this->filePath = "C://Users//cattaruz//Desktop//Diego//Instances//TWAVRP//";
		break;
	case lab:
		break;
	default:
		break;
	}

	switch (instanceType)
	{
	case twa:
		this->instanceFolder = "TWAVRPInstances//";
		break;
	case dtwa:
		this->instanceFolder = "DTWAVRPInstances//";
		break;
	default:
		break;
	}
}