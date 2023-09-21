#pragma once

#include <stdio.h>
#include "extern/json.hpp"

enum ACTYPE { Wall = 1, Ceiling = 2, Stand = 3, Tower = 4 };
enum WINDTYPE { MANUAL, POWER, FOREST, AUTO, AIRGUIDE, HIGHCEILING, FOCUS, WIDECARE = 1, SPACEDIV, XXXX, LEFT, RIGHT, DIV_L, DIV_R
};


struct AirConditioner {
	ACTYPE acType;
	int pos[3];
	float rot[3];

	WINDTYPE windType;
	int ventLevel; //up-mid-down-swing
	int ventSpeed; //low-mid-high

	int direction;	//0:+x 1:+x-z 2:-z 3:-x-z 4:-x 5:-x+z 6:+z 7:+x+z (wind direction)

	//float curTemp;
	float tarTemp;
};

struct MachineSetting
{
	int type = 2;
	int posX = 32;
	int posY = 32;
	int posZ = 32;

	int ventDir = 2;
	int ventSpeed = 1;
	int ventDist = 4; // distance from center;

	float tarTem = 21.f;	// raw data
	float target = 0.7f;	// normal data

	WINDTYPE windType;
	void setWindType(WINDTYPE wind);
	int direction;

	float initAngle = 20.f;
	float interval = 10.f;
	float swingSpeed = 0.001f;

	void postProcess(float min, float max)
	{
		target = 1.f - (tarTem - min) / (max - min);
	}
};

struct EulerianSetting
{
	size_t gridX = 64, gridY = 32, gridZ = 64;
	float timeStep = 0.005f;
	float diffuseCoef = 1.0e-6f;   // for material
	float viscosityCoef = 1.0e-7f; // for velocity
	float buoyancyCoef = 0.007f;
	float gravity = -9.8f;
	float force = 0.83f;

	// raw data
	float srcTem = 15.f;
	float minTem = 12.f;
	float maxTem = 31.f;
	float curTem = 30.f;

	// normal data
	float source = 15.f;
	float rstTem = 30.f; // rest temperature

	void postProcess()
	{
		//source = 1.f - (srcTem - minTem) / (maxTem - minTem);
		//rstTem = 1.f - (curTem - minTem) / (maxTem - minTem);
		source = srcTem;
		rstTem = curTem;
		//printf("min max src cur: %f %f %f %f\n", minTem, maxTem, srcTem, curTem);
		//printf("src: %f, rst: %f\n", source, rstTem);
	}
};


struct EnvSetting
{
	bool saveData = false;
	bool saveAvg = false;
	bool temHistory = false;
	int saveInterval = 100;
	float scale = 10.f;
	float fastScale = 4.f;
	float endTime = 0.0f;

	float sAttenuationR = 0.8f;  // rate
	float sAttenuationP = 100.f;  // power

	bool useUpSampling = true;
	float drawScale = 0.1f;
	bool drawSecond = true;
};


class Setting
{
public:
	Setting();
	~Setting();
	static Setting* instance;

	//static Setting* Instance()
	//{
	//	if (instance == NULL)
	//		instance = new Setting();

	//	return instance;
	//}

	bool LoadSetting(const char* fn);

	EulerianSetting eulerianSetting;
	EnvSetting envSetting;
	int machineCount;
	MachineSetting* machineSetting;
	std::string settingFileLocation;

private:

	bool LoadEulerain(nlohmann::json& jsonData);
	bool LoadEnvironment(nlohmann::json& jsonData);
	bool LoadMachine(nlohmann::json& jsonData);

	template <typename T>
	static bool ReadValue(const nlohmann::json& j, const std::string& key, T& v) {
		if (j.find(key) == j.end())
			return false;
		v = j[key].get<T>();
		return true;
	}
};

