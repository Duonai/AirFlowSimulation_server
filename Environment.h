#pragma once

#include <iostream>
#include <fstream>


using namespace std;


class EnvironmentBase
{
public : 
	string name;
};


class FluidsEnv : public EnvironmentBase 
{
public:



};


class Environment
{
public : 
	static Environment& instance() {
		
		if (instance_ = nullptr)
			instance_ = new Environment;
			
		return *instance_;
	}

private :
	Environment(){};

	void Load();
	

private : 
	static Environment* instance_;

};


class LoadParam
{
public :
	LoadParam(string filename)
	{
	
	};

private:

};

void readSection();
void readVariable();

