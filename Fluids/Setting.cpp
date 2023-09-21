#include "Setting.h"

#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>

Setting* Setting::instance = NULL;

void MachineSetting::setWindType(WINDTYPE wind)
{
	windType = wind;
}

Setting::Setting()
{
	machineCount = 1;
	machineSetting = new MachineSetting[machineCount];
}

Setting::~Setting()
{

}


bool Setting::LoadSetting(const char* fn)
{
	nlohmann::json m_json;
	std::ifstream input_file(fn);

	if (!input_file.is_open())
	{
		printf("setting file loading error\n");
		Setting s;
		*this = s;
		return false;
	}

	try {
		m_json << input_file;
	}
	catch (const std::exception& e) {
		printf("%s\n", e.what());
		exit(1);
	}

	settingFileLocation = fn;

	if (m_json.find("Eulerian") != m_json.end()) {
		nlohmann::json eul = m_json["Eulerian"];
		LoadEulerain(eul);
	}

	if (m_json.find("Environment") != m_json.end()) {
		nlohmann::json env = m_json["Environment"];
		LoadEnvironment(env);
	}

	if (m_json.find("Machine") != m_json.end()) {
		nlohmann::json mac = m_json["Machine"];
		LoadMachine(mac);
	}

	eulerianSetting.postProcess();
	for (int i = 0; i < machineCount; i++)
		machineSetting[i].postProcess(eulerianSetting.minTem, eulerianSetting.maxTem);

	return true;
}

bool Setting::LoadEulerain(nlohmann::json& jsonData)
{
	ReadValue(jsonData, "gridX", eulerianSetting.gridX);
	ReadValue(jsonData, "gridY", eulerianSetting.gridY);
	ReadValue(jsonData, "gridZ", eulerianSetting.gridZ);
	ReadValue(jsonData, "timeStep", eulerianSetting.timeStep);
	ReadValue(jsonData, "diffuse", eulerianSetting.diffuseCoef);
	ReadValue(jsonData, "viscosity", eulerianSetting.viscosityCoef);
	ReadValue(jsonData, "buoyancy", eulerianSetting.buoyancyCoef);
	ReadValue(jsonData, "gravity", eulerianSetting.gravity);
	ReadValue(jsonData, "force", eulerianSetting.force);
	ReadValue(jsonData, "srcTem", eulerianSetting.srcTem);
	ReadValue(jsonData, "minTem", eulerianSetting.minTem);
	ReadValue(jsonData, "maxTem", eulerianSetting.maxTem);
	ReadValue(jsonData, "curTem", eulerianSetting.curTem);

	//printf("Eul Setting---------------\n");
	//printf("grid (%d,%d,%d)\n", eulerianSetting.gridX, eulerianSetting.gridY, eulerianSetting.gridZ);
	//printf("timestep %f, diffuse %f, viscosity %f,\n", eulerianSetting.timeStep, eulerianSetting.diffuseCoef, eulerianSetting.viscosityCoef);
	//printf("buoyancy %f, gravity %f, force %f,\n", eulerianSetting.buoyancyCoef, eulerianSetting.gravity, eulerianSetting.force);
	//printf("srcTem %f, minTem %f, maxTem %f, curTem %f,\n", eulerianSetting.srcTem, eulerianSetting.minTem, eulerianSetting.maxTem, eulerianSetting.curTem);

	return true;
}

bool Setting::LoadEnvironment(nlohmann::json& jsonData)
{
	ReadValue(jsonData, "saveData", envSetting.saveData);
	ReadValue(jsonData, "saveAvg", envSetting.saveAvg);
	ReadValue(jsonData, "temHistory", envSetting.temHistory);
	ReadValue(jsonData, "saveInterval", envSetting.saveInterval);
	ReadValue(jsonData, "scale", envSetting.scale);
	ReadValue(jsonData, "fastScale", envSetting.fastScale);
	ReadValue(jsonData, "endtime", envSetting.endTime);
	ReadValue(jsonData, "sAttenuationR", envSetting.sAttenuationR);
	ReadValue(jsonData, "sAttenuationP", envSetting.sAttenuationP);
	ReadValue(jsonData, "useUpSampling", envSetting.useUpSampling);
	ReadValue(jsonData, "drawSecond", envSetting.drawSecond);
	ReadValue(jsonData, "drawScale", envSetting.drawScale);

	//printf("Env Setting---------------\n");
	//printf("saveData %s, saveAvg %s, temHistory %s,\n", envSetting.saveData ? "true" : "false", envSetting.saveAvg ? "true" : "false", envSetting.temHistory ? "true" : "false");
	//printf("saveInterval %d, scale %f, fastScale %f, endtime %f,\n", envSetting.saveInterval, envSetting.scale, envSetting.fastScale, envSetting.endTime);
	//printf("sAttenuationR %f, sAttenuationP %f, useUpSampling %s,\n", envSetting.sAttenuationR, envSetting.sAttenuationP, envSetting.useUpSampling ? "true" : "false");
	//printf("drawSecond %s, drawScale %f,\n", envSetting.drawSecond ? "true" : "false", envSetting.drawScale);

	return true;
}

bool Setting::LoadMachine(nlohmann::json& jsonData)
{
	ReadValue(jsonData, "machineCount", machineCount);
	machineSetting = new MachineSetting[machineCount];

	if (jsonData.find("MachineInfo") == jsonData.end())
		return false;

	std::vector<nlohmann::json> jsonData_mac;
	if (!ReadValue(jsonData, "MachineInfo", jsonData_mac))
		return false;

	
	for (int i = 0; i < machineCount; i++) {
		nlohmann::json& mac_i = jsonData_mac[i];
		MachineSetting& macSetting_i = machineSetting[i];
		
		char type[12];
		if (mac_i.find("type") == mac_i.end()) {
			printf("Unknown air conditioner type\n");
			continue;
		}
		else {
			strcpy(type, mac_i["type"].get<std::string>().c_str());
			if (strcmp(type, "wall") == 0)
				macSetting_i.type = ACTYPE::Wall;
			else if (strcmp(type, "ceiling") == 0)
				macSetting_i.type = ACTYPE::Ceiling;
			else if (strcmp(type, "stand") == 0)
				macSetting_i.type = ACTYPE::Stand;
			else if (strcmp(type, "tower") == 0)
				macSetting_i.type = ACTYPE::Tower;
			else
				printf("Unknown air conditioner type\n");
		}

		if (mac_i.find("windType") == mac_i.end())
		{
			macSetting_i.windType = WINDTYPE::MANUAL;
		}
		else
		{
			strcpy(type, mac_i["windType"].get<std::string>().c_str());

			if (strcmp(type, "power") == 0)
				macSetting_i.windType = WINDTYPE::POWER;
			else if (strcmp(type, "highceiling") == 0)
				macSetting_i.windType = WINDTYPE::HIGHCEILING;
			else if (strcmp(type, "airguide") == 0)
				macSetting_i.windType = WINDTYPE::AIRGUIDE;
			else if (strcmp(type, "forest") == 0)
				macSetting_i.windType = WINDTYPE::FOREST;
			else if (strcmp(type, "auto") == 0)
				macSetting_i.windType = WINDTYPE::AUTO;
			else if (strcmp(type, "manual") == 0)
				macSetting_i.windType = WINDTYPE::MANUAL;

			else if (strcmp(type, "widecare") == 0)
				macSetting_i.windType = WINDTYPE::WIDECARE;
			else if (strcmp(type, "spacediv") == 0)
				macSetting_i.windType = WINDTYPE::SPACEDIV;
			else if (strcmp(type, "xxxx") == 0)
				macSetting_i.windType = WINDTYPE::XXXX;
			else if (strcmp(type, "left") == 0)
				macSetting_i.windType = WINDTYPE::LEFT;
			else if (strcmp(type, "right") == 0)
				macSetting_i.windType = WINDTYPE::RIGHT;
			else if (strcmp(type, "div_l") == 0)
				macSetting_i.windType = WINDTYPE::DIV_L;
			else if (strcmp(type, "div_r") == 0)
				macSetting_i.windType = WINDTYPE::DIV_R;
			else{
				printf("Unknown air conditioner windType\n");
				macSetting_i.windType = WINDTYPE::MANUAL;
			}
		}

		ReadValue(mac_i, "posX", macSetting_i.posX);
		ReadValue(mac_i, "posY", macSetting_i.posY);
		ReadValue(mac_i, "posZ", macSetting_i.posZ);
		ReadValue(mac_i, "ventDir", macSetting_i.ventDir);
		ReadValue(mac_i, "ventSpeed", macSetting_i.ventSpeed);
		ReadValue(mac_i, "ventDist", macSetting_i.ventDist);
		ReadValue(mac_i, "swingSpeed", macSetting_i.swingSpeed);
		ReadValue(mac_i, "tarTem", macSetting_i.tarTem);
		ReadValue(mac_i, "initAngle", macSetting_i.initAngle);
		ReadValue(mac_i, "interval", macSetting_i.interval);
		ReadValue(mac_i, "direction", macSetting_i.direction);

		if(macSetting_i.type==ACTYPE::Ceiling){
			if (macSetting_i.windType == WINDTYPE::POWER) {
				macSetting_i.ventSpeed = 3;
			}
			else if (macSetting_i.windType == WINDTYPE::HIGHCEILING) {
				macSetting_i.ventSpeed = 3;
			}
			else if (macSetting_i.windType == WINDTYPE::AIRGUIDE) {
				macSetting_i.ventSpeed = 3;
			}
			else if (macSetting_i.windType == WINDTYPE::AUTO) {
				macSetting_i.ventSpeed = 3;
			}
			else if (macSetting_i.windType == WINDTYPE::FOREST) {
				macSetting_i.ventSpeed = 2;
			}
		}

		//printf("%d\n", macSetting_i.windType);
		//printf("Mac Setting---------------\n");
		//printf("saveData %d, pos (%d,%d,%d),\n", macSetting_i.type, macSetting_i.posX, macSetting_i.posY, macSetting_i.posZ);
		//printf("ventDir %d, ventSpeed %d, ventDist %d,\n", macSetting_i.ventDir, macSetting_i.ventSpeed, macSetting_i.ventDist);
		//printf("swingSpeed %f, tarTem %f, initAngle %f, interval %f,\n", macSetting_i.swingSpeed, macSetting_i.tarTem, macSetting_i.initAngle, macSetting_i.interval);
	}

	

	return true;
}

