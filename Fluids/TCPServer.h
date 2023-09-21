#pragma once

#include "extern/Packet.h"

#include <vector>
#include <thread>
#include <mutex>

#include <WinSock2.h>
#include "solver.cuh"

using namespace tcp_boost;
using namespace std;


class TCPServer
{
public:
	TCPServer();
	~TCPServer();

	void memoryFree();
	void open();
	void start();
	void update(int solvernum, const float* velocity, const float* temperature);

	int getGridX(int solvernum);
	int getGridY(int solvernum);
	int getGridZ(int solvernum);
	//const bool* getOccupancyGrid();
	//const int* getACPos();
	//const float* getACVec();
	//const float getCurTem();
	//const float getTarTem();
	const vector<thread*> getClientList();

	void setSolver(Solver** sol, int n);
	//bool isAllocated();
	bool isSimReset();
	void setSimReset(bool val);

	void setTimeFlies();
	void reAllocate(int solvernum, int x, int y, int z);

private: 
	Packet printChar(Packet& str);
	Packet processSpace(int solvernum, Packet& str);
	Packet processOccupancy(int solvernum, Packet& str);
	Packet processACInfo(int solvernum, Packet& str);
	Packet processMakeData(int solvernum);
	Packet processMakeDataNew(int solvernum);
	Packet processStoreBIM(Packet& str);
	Packet processBIM(int solvernum, Packet &str);
	Packet processGraphData(int solvernum, Packet& str);
	Packet process(int solvernum, Packet& str);

	int getnewsolver();

	void work(int solvernum, SOCKET clientSock, SOCKADDR_IN clientAddr);

	void reset(int solvernum);

	void index3D(int solvernum, int idx, int& x, int& y, int& z, float scale = 1.f);
	int index1D(int solvernum, int x, int y, int z);
private:
	bool initialized;
	//bool allocated;
	bool simReset;
	//bool occInitialized;

	vector<thread*> clientList;
	SOCKET serverSock;
	mutex mut;
	//mutex procmut;

	//float* vel;
	//float* tem;
	//bool* occ;
	//int ACPosPre[3];
	//int ACPos[3];
	//float ACVec[3];
	Solver* solver[SOLVER_NUM];
	bool solverOccupied[SOLVER_NUM];

	//float curTem;	// Current Temperature
	//float tarTem;	// Target Temperature

	//int gridSize;
	//int gridx, gridy, gridz;
	//int rawx, rawy, rawz;
	//int upsampleScale;

	// timelater 
	int timeLater;
	float startTime;
};

string getCurrentTime();
