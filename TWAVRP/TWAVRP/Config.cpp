#include "Config.h"

Config::Config(user user, instance instance) {
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

	switch (instance)
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