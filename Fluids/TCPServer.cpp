#pragma comment(lib, "ws2_32")
#pragma warning(disable:4996)

#include <stdio.h>
#include "TCPServer.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime> 
#include <sstream>
#include <math.h>
#include <exception>
#include <direct.h>

#include<fstream>
#include<string>


using std::this_thread::sleep_for;
using namespace std;

using namespace tcp_boost;
using namespace std::chrono;

#define BUFFERSIZE 1024

enum packetType { occupancy = 0x11, requestData = 0x12, printChar = 0x13, ACInfo = 0x14, saveBIM = 0x16, loadBIM = 0x17, drawGraph = 0x18, None = 0xff };
mutex sendmut;

TCPServer::TCPServer()
{
	initialized = false;
	//allocated = false;
	simReset = false;
	//occInitialized = false;
	
	startTime = 0.0f;

	//ACPosPre[0] = 0;
	//ACPosPre[1] = 0;
	//ACPosPre[2] = 0;

}

TCPServer::~TCPServer()
{
	memoryFree();


	for (int i = 0; i < clientList.size(); i++)
	{
		delete clientList[i];
	}
}

void TCPServer::memoryFree()
{
	//if (vel)
	//{
	//	delete vel;
	//	vel = nullptr;
	//}
	//if (tem)
	//{
	//	delete tem;
	//	tem = nullptr;
	//}
	//if (occ)
	//{
	//	delete occ;
	//	occ = nullptr;
	//}
}

void TCPServer::open()
{
	WSADATA wsaData; //소켓 데이터
	if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0) //소켓 실행
	{
		std::cout << "Socket - Initialization fail" << endl;
		return;
	}

	serverSock = socket(PF_INET, SOCK_STREAM, 0); //소켓 생성
	SOCKADDR_IN addr; //소켓 주소

	memset(&addr, 0, sizeof(addr));
	addr.sin_family = AF_INET;
	addr.sin_addr.s_addr = htonl(INADDR_ANY);
	addr.sin_port = htons(1235);

	if (::bind(serverSock, (SOCKADDR*)&addr, sizeof(SOCKADDR_IN)) == SOCKET_ERROR)
	{
		std::cout << "Socket - Bind error" << endl;
		return;
	}
	if (listen(serverSock, SOMAXCONN) == SOCKET_ERROR)
	{
		std::cout << "Socket - Listen error" << endl;
		return;
	}

	//curTem = 0.f;
	//tarTem = -1.f;

	initialized = true;
}

void TCPServer::reAllocate(int solvernum, int x, int y, int z)
{
	//allocated = false;
	
	//rawx = x;
	//rawy = y;
	//rawz = z;
	const int max = 128 * 128 * 128;
	int upsampleScale, gridSize;
	int gridx, gridy, gridz;

	//if (solver[solvernum]->setting->envSetting.useUpSampling)
	//{
	//	upsampleScale = cbrt(max / (x * y * z));
	//	if (upsampleScale < 1)
	//		upsampleScale = 1;
	//}
	//else
	//{
	//	upsampleScale = 1;
	//}
	upsampleScale = 1;
	solver[solvernum]->setUpSamplingScale(upsampleScale);
		
	gridx = x * upsampleScale;
	gridy = y * upsampleScale;
	gridz = z * upsampleScale;

	printf("----Client grid : %d x %d x %d\n", x, y, z);
	printf("----Server grid : %d x %d x %d\n", gridx, gridy, gridz);

	//gridSize = gridx * gridy * gridz;

	//try
	//{
	//	mut.lock();
	//	memoryFree();
	//	vel = new float[gridSize * 3];
	//	tem = new float[gridSize];
	//	//occ = new bool[gridSize];
	//	mut.unlock();
	//}
	//catch (exception& ex)
	//{
	//	printf("-- TCPServer.reAllocate (vel, temp, occ) Error\n%s\n", ex.what());
	//}

	//try
	//{
	//	memset(vel, 0, sizeof(bool) * gridSize * 3);
	//	memset(tem, 0, sizeof(bool) * gridSize);
	//	//memset(occ, 0, sizeof(bool) * gridSize);
	//}
	//catch (exception& ex)
	//{
	//	printf("-- TCPServer.reAllocate memset Error\n%s\n", ex.what());
	//}
	solver[solvernum]->init(gridx, gridy, gridz);
	//solver[solvernum]->setUpSamplingScale(upsampleScale);

	simReset = true;

	printf("-- Server is allocated..\n");
}

void TCPServer::index3D(int solvernum, int idx, int& x, int& y, int& z, float scale)
{
	int rawx = solver[solvernum]->nx / solver[solvernum]->upsampleScale;
	int rawy = solver[solvernum]->ny / solver[solvernum]->upsampleScale;
	int rawz = solver[solvernum]->nz / solver[solvernum]->upsampleScale;

	z = (idx % rawz) * scale;
	y = ((idx / rawz) % rawy) * scale;
	x = (idx / (rawz * rawy)) * scale;
}

int TCPServer::index1D(int solvernum, int x, int y, int z)
{
	return z + y * solver[solvernum]->nz + x * solver[solvernum]->nz * solver[solvernum]->ny;
}

void TCPServer::start()
{
	if (!initialized)
	{
		printf("Socket is not initialized\n");
		return;
	}

	printf("Socket - Start\n");

	while (1) //다중 접속 & 서버 유지
	{
		int len = sizeof(SOCKADDR_IN);
		SOCKADDR_IN clientAddr;
		SOCKET clientSock = accept(serverSock, (SOCKADDR*)&clientAddr, &len); //접속한 클라 소켓 받음
		clientList.push_back(new thread(&TCPServer::work, this, getnewsolver(), clientSock, clientAddr)); //쓰레드 리스트에 삽입
	}
}
int TCPServer::getnewsolver() {
	for (int i = 0; i < SOLVER_NUM; i++) {
		if (!solverOccupied[i]) {
			solverOccupied[i] = true;
			return i;
		}
	}
	return -1;
}
void TCPServer::update(int solvernum, const float* velocity, const float* temperature)
{
	//try
	//{
	//	mut.lock();
	//	memcpy(vel, velocity, solver[solvernum]->n3 * 3 * sizeof(float));
	//	memcpy(tem, temperature, solver[solvernum]->n3 * sizeof(float));
	//	mut.unlock();
	//}
	//catch (exception& ex)
	//{
	//	printf("-- TCPServer.update Error\n%s\n", ex.what());
	//}
}

int TCPServer::getGridX(int solvernum)
{
	return solver[solvernum]->nx;
}

int TCPServer::getGridY(int solvernum)
{
	return solver[solvernum]->ny;
}

int TCPServer::getGridZ(int solvernum)
{
	return solver[solvernum]->nz;
}

//const bool* TCPServer::getOccupancyGrid()
//{
//	return occ;
//}
//
//const int* TCPServer::getACPos()
//{
//	return ACPos;
//}
//
//const float* TCPServer::getACVec()
//{
//	return ACVec;
//}
//
//const float TCPServer::getCurTem()
//{
//	return curTem;
//}
//
//const float TCPServer::getTarTem()
//{
//	return tarTem;
//}

const vector<thread*> TCPServer::getClientList()
{
	return clientList;
}

void TCPServer::setSolver(Solver** sol, int n)
{
	if (n > SOLVER_NUM)
		n = SOLVER_NUM;
	for (int i = 0; i < n; i++) {
		solver[i] = sol[i];
		solverOccupied[i] = false;
	}
}

//bool TCPServer::isAllocated()
//{
//	return allocated;
//}

bool TCPServer::isSimReset()
{
	return simReset;
}

void TCPServer::setSimReset(bool val)
{
	simReset = val;
}

void TCPServer::setTimeFlies()
{
	timeLater = 0;
}

Packet TCPServer::printChar(Packet& str)
{
	int size = str.get_body_length();
	char* ret = new char[size - 1 + 10];

	cout << size << endl;

	memcpy(ret, "echo - ", 7);

	memcpy(ret + 7, str.get_body() + 1, size - 1); //type 용 1byte 제외
	memcpy(ret + 7 + size - 1, "\0", 1);

	std::cout << "From Client message : ";
	for (int i = 0; i < size; i++)
	{
		cout << ret[7 + i];
	}

	std::cout << endl;
	memcpy(ret + 7 + size - 1, "\n>\0", 3);

	Packet tmp = Packet(size + 10 + sizeof(int) + 1);
	tmp.push_byte(packetType::printChar);
	//cout << tmp.get_buffer_length() << endl;

	for (int i = 0; i < size + 10; i++)
	{
		tmp.push_byte(ret[i]);
	}

	delete[] ret;
	return tmp;
}

Packet TCPServer::processSpace(int solvernum, Packet& str) //수신 데이터를 콘솔 출력 및 echo 데이터 생성용 함수
{
	int size = str.get_body_length() - 1;
	char* ret = new char[20];

	size_t x = str.pop_int32();
	size_t y = str.pop_int32();
	size_t z = str.pop_int32();

	reAllocate(solvernum, x, y, z);

	memcpy(ret, "Space success\n", 20);

	Packet tmp = Packet(20 + sizeof(int) + 1);
	tmp.push_byte(packetType::occupancy);

	for (int i = 0; i < 20; i++)
	{
		tmp.push_byte(ret[i]);
	}

	delete[] ret;
	return tmp;
}

Packet TCPServer::processOccupancy(int solvernum, Packet& str) //수신 데이터를 콘솔 출력 및 echo 데이터 생성용 함수
{
	int size = str.get_body_length() - 1;
	char* ret = new char[20];
	//mut.lock();
	// 수신 데이터가 버퍼 사이즈를 넘었을 경우.
	int upsampleScale = solver[solvernum]->upsampleScale;


	if (solver[solvernum]->occInitialized)
	{
		//bool* occ = solver[solvernum]->h_obstacle;

		// update occupancy using index;
		//const std::lock_guard<std::mutex> lock(mut);
		for (int i = 0; i < size / 4; i++)
		{
			auto idx = str.pop_uint32();

			int ux, uy, uz;
			index3D(solvernum, idx, ux, uy, uz, upsampleScale);

			for (int i = ux; i < ux + upsampleScale; i++) {
				for (int j = uy; j < uy + upsampleScale; j++) {
					for (int k = uz; k < uz + upsampleScale; k++) {

						solver[solvernum]->flipObstacle(i + 1, j + 1, k + 1);
					}
				}
			}
		}
	}
	else
	{
		bool* occ = new bool[solver[solvernum]->n3];
		memset(occ, 0, solver[solvernum]->n3 * sizeof(bool));
		// init occupancy using index;
		for (int idx = 0; idx < size; idx++)
		{
			int ux, uy, uz;
			index3D(solvernum, idx, ux, uy, uz, upsampleScale);
			bool val = str.pop_bool();

			for (int i = ux; i < ux + upsampleScale; i++) {
				for (int j = uy; j < uy + upsampleScale; j++) {
					for (int k = uz; k < uz + upsampleScale; k++) {
						
						occ[index1D(solvernum, i, j, k)] = val;
					}
				}
			}
		}
		solver[solvernum]->setObstacle(occ);
		solver[solvernum]->occInitialized = true;
		//solver[solvernum]->start();
	}
	//solver[solvernum]->setObstacle(occ);
	//mut.unlock();


#ifdef DEBUG
	for (int i = 0; i < size; i++)
		cout << out[i] << "\t";;
	cout << endl;
#endif

	memcpy(ret, "Occupancy success\n", 20);

	Packet tmp = Packet(20 + sizeof(int) + 1);
	tmp.push_byte(packetType::occupancy);

	for (int i = 0; i < 20; i++)
	{
		tmp.push_byte(ret[i]);
	}

	delete[] ret;
	return tmp;
}

Packet TCPServer::processACInfo(int solvernum, Packet& str)
{
	char* ret = new char[20];

	int tp = str.pop_int32();
	solver[solvernum]->setACType(0, tp);
	if (tp == ACTYPE::Ceiling) {
		//solver[solvernum]->
	}

	AirConditioner& ac0 = solver[solvernum]->airconditioners[0];

	int ACPos[3] = { 0, };
	int ACPosPre[3] = { 0, };	
	int ACDirPre = ac0.direction;
	ACPosPre[0] = ac0.pos[0];
	ACPosPre[1] = ac0.pos[1];
	ACPosPre[2] = ac0.pos[2];

	int dir = str.pop_int32();

	for (int j = 0; j < 3; j++)
	{
		ACPos[j] = (str.pop_int32() + 1/*(tp == ACTYPE::Tower || (tp == ACTYPE::Ceiling && j != 1) || tp == ACTYPE::Wall ? 1 : 0)*/) * solver[solvernum]->upsampleScale;
	}
	if (tp == ACTYPE::Ceiling && ACPos[1] >= solver[solvernum]->ny)
		ACPos[1] = solver[solvernum]->ny - 1;
	solver[solvernum]->setACPosition(0, ACPos[0], ACPos[1], ACPos[2]);
	bool acPosModified = (ACPosPre[0] != ACPos[0]) || (ACPosPre[1] != ACPos[1]) || (ACPosPre[2] != ACPos[2]) || ACDirPre != dir;

	int ventLevel = str.pop_int32();
	int ventSpeed = str.pop_int32();
	bool reset = str.pop_bool();
	//dir = 7;
	solver[solvernum]->setACDirection(0, dir);
	solver[solvernum]->setACVentDir(0, ventLevel);
	solver[solvernum]->setACVentMag(0, ventSpeed);

	/*if(ac0.acType==ACTYPE::Ceiling){
		if (ventLevel == 1) {
			solver[solvernum]->setACWindType(0, WINDTYPE::POWER);
			ac0.ventSpeed = 3;
		}
		else if (ventLevel == 2) {
			solver[solvernum]->setACWindType(0, WINDTYPE::FOREST);
			ac0.ventSpeed = 2;
		}
		else if (ventLevel == 3) {
			solver[solvernum]->setACWindType(0, WINDTYPE::AUTO);
			ac0.ventSpeed = 3;
		}
		else if (ventLevel == 4) {
			solver[solvernum]->setACWindType(0, WINDTYPE::AIRGUIDE);
			ac0.ventSpeed = 3;
		}
		else if (ventLevel == 5) {
			solver[solvernum]->setACWindType(0, WINDTYPE::HIGHCEILING);
			ac0.ventSpeed = 3;
		}
		else if (ventLevel == 6) {
			solver[solvernum]->setACWindType(0, WINDTYPE::FOCUS);
			ac0.ventSpeed = 3;
		}
		else {
			ac0.ventSpeed = ventSpeed;
		}
	}

	if (ac0.acType == ACTYPE::Tower) {
		if (ventLevel == 1)
			solver[solvernum]->setACWindType(0, WINDTYPE::WIDECARE);
		else if (ventLevel == 2)
			solver[solvernum]->setACWindType(0, WINDTYPE::SPACEDIV);
		else if (ventLevel == 3)
			solver[solvernum]->setACWindType(0, WINDTYPE::XXXX);
		else if (ventLevel == 4)
			solver[solvernum]->setACWindType(0, WINDTYPE::LEFT);
		else if (ventLevel == 5)
			solver[solvernum]->setACWindType(0, WINDTYPE::RIGHT);
		else if (ventLevel == 6)
			solver[solvernum]->setACWindType(0, WINDTYPE::DIV_L);
		else if (ventLevel == 7)
			solver[solvernum]->setACWindType(0, WINDTYPE::DIV_R);
		else
			solver[solvernum]->setACWindType(0, WINDTYPE::MANUAL);
	}*/


	if (reset || acPosModified)
	{
		solver[solvernum]->stop();
		solver[solvernum]->reset(false);
		solver[solvernum]->start();
	}
	//solver[solvernum]->start();
	//printf("--");
	solver[solvernum]->curTem = str.pop_single();
	solver[solvernum]->tarTem = str.pop_single();

	printf("\n-- Receive : type: %s // position: (%d, %d, %d) // level: %d // speed: %d // reset: %s\n", (tp == 1 ? "wall" : (tp == 2 ? "ceiling" : (tp == 3 ? "stand" : (tp == 4 ? "Tower" : "unknown")))), ACPos[0], ACPos[1], ACPos[2], ventLevel, ventSpeed, reset ? "true" : "false");
	memcpy(ret, "AC Info success\n", 20);

	Packet tmp = Packet(20 + sizeof(int) + 1);
	tmp.push_byte(packetType::ACInfo);

	for (int i = 0; i < 20; i++)
	{
		tmp.push_byte(ret[i]);
	}

	delete[] ret;
	return tmp;
}

Packet TCPServer::processMakeData(int solvernum)
{
	while (timeLater != 0)
	{
		// simulation is more faster realdata than 15.8 times. 
		auto curSec = solver[solvernum]->getSimulationTime();
		if (timeLater * 60 < (curSec - startTime) * 15.f)
		{
			auto scale = solver[solvernum]->setting->envSetting.fastScale;
			solver[solvernum]->setTimeStep(solver[solvernum]->getTimeStep() / scale);
			solver[solvernum]->setForce(solver[solvernum]->getForce() / scale);
			timeLater = 0;
			printf("\n\n End fast forward...........\n");
			break;
		}
		else
		{
			printf("\n-- %3.2f => %3.2f (%3.2f  x 15) \n", startTime, curSec, timeLater * 4.f);
		}
		Sleep(500);  // _sleep
	}
	//int acx = solver[solvernum]->airconditioners[0].pos[0] - 7;
	//int acy = solver[solvernum]->airconditioners[0].pos[1] - 3;
	//int acz = solver[solvernum]->airconditioners[0].pos[2] - 5;
	//printf("\nvelocity: (%f,%f,%f) ", solver[solvernum]->getVelocityGPUU(acx, acy, acz), solver[solvernum]->getVelocityGPUV(acx, acy, acz), solver[solvernum]->getVelocityGPUW(acx, acy, acz));
	//printf("temperature: %f", solver[solvernum]->getDensityGPU(acx, acy, acz));

	float* vel = new float[solver[solvernum]->n3 * 3];
	float* tem = new float[solver[solvernum]->n3];

	solver[solvernum]->getVelocityGPUN(vel);
	solver[solvernum]->getDensityGPUN(tem);


	const int gridsizeDown = solver[solvernum]->n3 / (solver[solvernum]->upsampleScale * solver[solvernum]->upsampleScale * solver[solvernum]->upsampleScale);
	//const int gridsizeDown = gridSize / (upsampleScale * upsampleScale * upsampleScale);
	size_t dirSize = (gridsizeDown * 1.5f + 0.6f); // odd number will be made even number adding 0.6;
	//unique_ptr<unsigned char[]> velDir(new unsigned char[dirSize]);
	//unique_ptr<unsigned char[]> velMag(new unsigned char[gridsize3]);
	//unique_ptr<unsigned char[]> temp(new unsigned char[gridsize3]);


	unsigned char *velDir = new unsigned char[dirSize];
	unsigned short *velMag2 = new unsigned short[gridsizeDown];
	//unsigned char* velMag1 = new unsigned char[gridsize3];
	unsigned char *temp   = new unsigned char[gridsizeDown];

	auto norDir = [](float v, float length)  -> float {
		return (v / length + 1.f) / 2.f * 15.f + 0.5f;
	};

#pragma omp parallel for
	for (int p = 0; p < gridsizeDown; p++)
	{
		float x0, y0, z0;
		float x1, y1, z1;
		float length;

		unsigned char ux0, uy0, uz0;
		unsigned char ux1, uy1, uz1;

		// 4bit =  unit a axis value; 0~15
		// 12 bit = a set 
		// 24 bit = 3byte

		int px, py, pz, up;
		index3D(solvernum, p, px, py, pz, solver[solvernum]->upsampleScale);
		up = index1D(solvernum, px, py, pz);

		x0 = vel[3 * up];
		y0 = vel[3 * up + 1];
		z0 = vel[3 * up + 2];
		length = sqrt(x0 * x0 + y0 * y0 + z0 * z0);

		ux0 = unsigned int(norDir(x0, length));
		uy0 = unsigned int(norDir(y0, length));
		uz0 = unsigned int(norDir(z0, length));

		if (p % 2 == 0)
		{
			velDir[3 * p / 2] = (ux0 << 4) | uy0;
			velDir[3 * p / 2 + 1] = (uz0 << 4);
		}
		else
		{
			velDir[3 * p / 2] |= ux0;
			velDir[3 * p / 2 + 1] = (uy0 << 4) | uz0;
		}

		velMag2[p] = unsigned int(length * 30000) < 65001 ? unsigned short(length * 30000) : 65000;
	}

	EulerianSetting& eu = solver[solvernum]->setting->eulerianSetting;

#pragma omp parallel for
	for (int i = 0; i < gridsizeDown; i++)
	{
		int px, py, pz, up;
		index3D(solvernum, i, px, py, pz);
		up = index1D(solvernum, px * solver[solvernum]->upsampleScale, py * solver[solvernum]->upsampleScale, pz * solver[solvernum]->upsampleScale);
		

		tem[up] = 1 - (tem[up] - eu.minTem) / (eu.maxTem - eu.minTem);
		temp[i] = int(tem[up] * 256) < 256 ? int(tem[up] * 256) : 255;
	}

	Packet tmp = Packet( int(4.5f * gridsizeDown + 0.6f) + sizeof(int) + 1);  // odd number must will be even number adding 0.6;
	tmp.push_byte(packetType::requestData);

	tmp.push_byte_array(reinterpret_cast<char*>(velDir), dirSize);
	//tmp.push_byte_array(reinterpret_cast<char*>(velMag), gridsize3);
	tmp.push_byte_array(reinterpret_cast<char*>(velMag2), gridsizeDown * 2);
	tmp.push_byte_array(reinterpret_cast<char*>(temp), gridsizeDown);

	try {
		delete[] velDir;
		//delete velMag1;
		delete[] velMag2;
		delete[] temp;
		delete[] vel;
		delete[] tem;
	}
	catch (exception &ex)
	{
		printf("\n-- TCPServer.ProcessMakedata \nex.what()\n");
	}
	return tmp;
}

Packet TCPServer::processMakeDataNew(int solvernum)
{
	int dataSize = sizeof(float);
	Packet tmp = Packet(solver[solvernum]->n3 * dataSize * 4 + sizeof(int) + 1);
	tmp.push_byte(packetType::requestData);

	//for (int p = 0; p < solver[solvernum]->n3 * 3; p++)
	//{
	//	tmp.push_single(vel[p]);
	//}

	//for (int p = 0; p < solver[solvernum]->n3; p++)
	//{
	//	tmp.push_single(tem[p]);
	//}

	return tmp;
}

Packet TCPServer::processStoreBIM(Packet& str)
{
	uint16_t gridSizeX = str.pop_uint16();
	uint16_t gridSizeY = str.pop_uint16();
	uint16_t gridSizeZ = str.pop_uint16();

	float distXL = str.pop_single();
	float distXR = str.pop_single();
	float distYU = str.pop_single();
	float distYD = str.pop_single();

	float roomSizeX = str.pop_single();
	float roomSizeY = str.pop_single();
	float roomSizeZ = str.pop_single();

	uint16_t roomNumber = str.pop_uint16();

	int size = gridSizeX * gridSizeY * gridSizeZ;

	char* occForBIM = new char[size];

	for (int i = 0; i < size; i++)
	{
		occForBIM[i] = str.pop_bool() ? '1' : '0';
	}

	string fileName = to_string(roomNumber);
	fileName = "BIM/" + fileName + ".bin";
	_mkdir("BIM");
	/*
	ofstream fp;

	fp.open(fileName);  //= fopen(pchar, "w");

	fp.write(to_string(gridSizeX).c_str(), to_string(gridSizeX).size());
	fp.write(to_string(gridSizeY).c_str(), to_string(gridSizeY).size());
	fp.write(to_string(gridSizeZ).c_str(), to_string(gridSizeZ).size());

	fp.write(to_string(distX).c_str(), to_string(distX).size());
	fp.write(to_string(distY).c_str(), to_string(distY).size());

	fp.write(to_string(roomSizeX).c_str(), to_string(roomSizeX).size());
	fp.write(to_string(roomSizeY).c_str(), to_string(roomSizeY).size());
	fp.write(to_string(roomSizeZ).c_str(), to_string(roomSizeZ).size());

	fp.write(occForBIM, sizeof(occForBIM));

	fp.close();
	*/
	FILE* fp;

	fp = fopen(fileName.c_str(), "wb");

	fwrite(&gridSizeX, sizeof(uint16_t), 1, fp);
	fwrite(&gridSizeY, sizeof(uint16_t), 1, fp);
	fwrite(&gridSizeZ, sizeof(uint16_t), 1, fp);

	fwrite(&distXL, sizeof(float), 1, fp);
	fwrite(&distXR, sizeof(float), 1, fp);
	fwrite(&distYU, sizeof(float), 1, fp);
	fwrite(&distYD, sizeof(float), 1, fp);

	fwrite(&roomSizeX, sizeof(float), 1, fp);
	fwrite(&roomSizeY, sizeof(float), 1, fp);
	fwrite(&roomSizeZ, sizeof(float), 1, fp);

	fwrite(occForBIM, 1, size, fp);

	fclose(fp);

	char* ret = new char[17];
	memcpy(ret, "BIMsave success\n", 17);

	Packet tmp = Packet(18 + sizeof(int));
	tmp.push_byte(packetType::saveBIM);

	for (int i = 0; i < 17; i++)
	{
		tmp.push_byte(ret[i]);
	}

	delete fp;
	delete[] occForBIM;
	delete[] ret;

	return tmp;
}

Packet TCPServer::processBIM(int solvernum, Packet& str)
{
	uint16_t roomNumber = str.pop_uint16();
	string fileName = to_string(roomNumber);
	fileName = "BIM/" + fileName + ".bin";

	uint16_t gridSizeX = 0;
	uint16_t gridSizeY = 0;
	uint16_t gridSizeZ = 0;

	float distXL;
	float distXR;
	float distYU;
	float distYD;

	float roomSizeX;
	float roomSizeY;
	float roomSizeZ;

	FILE* fp;

	fp = fopen(fileName.c_str(), "rb");

	fread(&gridSizeX, sizeof(uint16_t), 1, fp);
	fread(&gridSizeY, sizeof(uint16_t), 1, fp);
	fread(&gridSizeZ, sizeof(uint16_t), 1, fp);

	fread(&distXL, sizeof(float), 1, fp);
	fread(&distXR, sizeof(float), 1, fp);
	fread(&distYU, sizeof(float), 1, fp);
	fread(&distYD, sizeof(float), 1, fp);

	fread(&roomSizeX, sizeof(float), 1, fp);
	fread(&roomSizeY, sizeof(float), 1, fp);
	fread(&roomSizeZ, sizeof(float), 1, fp);

	int size = gridSizeX * gridSizeY * gridSizeZ;

	char* occForBIM = new char[size];

	fread(occForBIM, sizeof(char), size, fp);

	fclose(fp);

	reAllocate(solvernum, gridSizeX, gridSizeY, gridSizeZ);

	Packet tmp = Packet(sizeof(uint16_t) * 3 + sizeof(float) * 7 + size + sizeof(int) + 1);

	tmp.push_byte(packetType::loadBIM);

	tmp.push_uint16(gridSizeX);
	tmp.push_uint16(gridSizeY);
	tmp.push_uint16(gridSizeZ);

	tmp.push_single(distXL);
	tmp.push_single(distXR);
	tmp.push_single(distYU);
	tmp.push_single(distYD);

	tmp.push_single(roomSizeX);
	tmp.push_single(roomSizeY);
	tmp.push_single(roomSizeZ);

	for (int i = 0; i < size; i++)
	{
		if (occForBIM[i] == '1')
		{
			tmp.push_byte((char)true);
		}
		else
		{
			tmp.push_byte((char)false);
		}
	}

	delete[] occForBIM;

	return tmp;
}

Packet TCPServer::processGraphData(int solvernum, Packet& str)
{
	uint16_t x = str.pop_uint16();
	uint16_t y = str.pop_uint16();
	uint16_t z = str.pop_uint16();

	auto datas = solver[solvernum]->getTempDatas();

	//auto conv1D = [&](int x, int y, int z) -> int {
	//	return z + y * solver[solvernum]->nz + x * solver[solvernum]->nz * solver[solvernum]->ny;
	//};

	auto tempConvert = [&](float t) -> float {
		auto set = solver[solvernum]->setting->eulerianSetting;
		return set.minTem + (set.maxTem - set.minTem) * (1.f-t);
	};

	int numData = 10;
	Packet tmp = Packet(sizeof(float) * numData + sizeof(int) + 1);
	tmp.push_byte(packetType::drawGraph);

	if (solver[solvernum]->setting->envSetting.temHistory)
	{
		printf("\n Send tmeperature history : ");

		// if the datas is not enough to request size (numData);
		float dx = ((datas.size() < numData) ? 1 : (datas.size() / float(numData)));
		int cell = index1D(solvernum, x, y, z);
		for (int i = 0; i < numData; i++)
		{
			int idx = i * dx;
			if (idx >= datas.size())
				idx = datas.size() - 1;

			auto t1 = datas[idx][cell];
			auto sendT = tempConvert(t1);
			tmp.push_single(sendT);

			printf("%2.2f /", sendT);
		}
		printf("\n");
	}
	else
	{
		printf("\n the temperature history option is \"false\". Turn on that option\n");
		for (int i = 0; i < numData; i++)
		{
			tmp.push_single(0.f);
		}
	}

	return tmp;
}

Packet TCPServer::process(int solvernum, Packet& str)
{
	//lock_guard<mutex> lock(procmut);
	try
	{
		char type = str.pop_byte();
		if(type!=(char)18)
			printf("\n%d process %d ", solvernum, type);
		switch (type)
		{
		case 0x10: // 16
			return processSpace(solvernum, str);
		case 0x11: // 17
			//printf("post-Occupancy\n");
			return processOccupancy(solvernum, str);
		case 0x12: // 18
			//printf("post-makeData");
			printf(".");
			timeLater = (int)str.pop_byte(); //need unsigned check
			if (timeLater > 0)
			{
				auto scale = solver[solvernum]->setting->envSetting.fastScale;

				startTime = solver[solvernum]->getSimulationTime();
				solver[solvernum]->setTimeStep(solver[solvernum]->getTimeStep() * scale);
				solver[solvernum]->setForce(solver[solvernum]->getForce() * scale);
				printf("\n\n Start fast forward...........\n");
			}
			return processMakeData(solvernum);
		case 0x13: // 19
			//printf("post-char\n");
			return printChar(str);
		case 0x14: // 20
			//printf("post-acInfo\n");
			return processACInfo(solvernum, str);
		case 0x16:
			return processStoreBIM(str);
		case 0x17:
			return processBIM(solvernum, str);
		case 0x18:
			return processGraphData(solvernum, str);
		default:
			std::cout << "type error! type: " << type << endl;
			return nullptr;
		}
	}
	catch (exception& ex)
	{
		printf("-- TCPServer.Work process Error\n%s\n", ex.what());
	}
	return nullptr;
}

string getCurrentTime()
{
	auto now = std::chrono::system_clock::now();
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
	auto in_time_t = std::chrono::system_clock::to_time_t(now);
	tm* tstruct = new tm();
	localtime_s(tstruct, &in_time_t);

	std::stringstream ss;
	ss << std::put_time(tstruct, "%Y-%m-%d %X");
	ss << ":" << std::setw(3) << std::setfill('0') << ms.count();

	delete tstruct;
	return ss.str();
}


void TCPServer::work(int solvernum, SOCKET clientSock, SOCKADDR_IN clientAddr)
{
	solver[solvernum]->occInitialized = false;
	//bool flag = true;
	//setsockopt(clientSock, IPPROTO_TCP, TCP_NODELAY, (char*)&flag, sizeof(flag));

	//u_int mode = 1;  // 1 to enable non-blocking socket
	//ioctlsocket(clientSock, FIONBIO, &mode);

	auto curTime = time(NULL);
	struct tm* pLocal = localtime(&curTime);

	printf("[%s] Client connected IP address = %s : %d\n", getCurrentTime().c_str(), inet_ntoa(clientAddr.sin_addr), ntohs(clientAddr.sin_port));

	int curPos = 0;
	unsigned int len = 0;
	bool lenDone = false;

	unsigned char* rcvbuf = new unsigned char[1048600];
	char* tmpbuf = new char[1048600];

	//printf("Start Loop\n");

	while (1)
	{
		Sleep(10);
		//printf("l");
		if ((!lenDone && curPos < 4) || (lenDone && curPos < len + sizeof(int))) {
			int recvLen;
			try
			{
				//printf("r");
				recvLen = recv(clientSock, tmpbuf, 262150, 0); // maximmum for occupancy  4 + 1 + 64 * 64* 64 (+ 1)
				//printf(" %d ",recvLen);
			}
			catch (exception& ex)
			{
				printf("-- TCPServer.Work Receive Error\n%s\n", ex.what());
				break;
			}

			if (recvLen == 0)
			{
				cout << "\nGraceful Socket Close" << endl;
				break;
			}
			else if (recvLen < 0)
			{
				cout << "\nProblem in Ln: " << __LINE__ << ", Error Code: " << WSAGetLastError() << endl; //10014
				break;
			}

			memcpy(rcvbuf + curPos, tmpbuf, recvLen);
			curPos += recvLen;
		}

		if (!lenDone && curPos >= 4)
		{
			//memcpy(&len, rcvbuf, sizeof(short));
			len = (rcvbuf[0] << 24) | (rcvbuf[1] << 16) | (rcvbuf[2] << 8) | rcvbuf[3];
			lenDone = true;
		}

		// post processing
		if (lenDone && curPos >= len + sizeof(int))
		{
			Packet recvData = Packet(rcvbuf, sizeof(int) + len);
			try
			{
				Packet echo = process(solvernum, recvData);
				int real_len = (echo.get_buffer()[0] << 24) | (echo.get_buffer()[1] << 16) | (echo.get_buffer()[2] << 8) | echo.get_buffer()[3];
				//if (send(clientSock, echo.get_buffer(), real_len + sizeof(int), 0) == SOCKET_ERROR)
				//printf("sendmut ");
				sendmut.lock();
				//printf("lock ");
				if (send(clientSock, echo.get_buffer(), echo.get_total_length(), 0) == SOCKET_ERROR)
				{
					cout << "\nProblem in Ln: " << __LINE__ << ", Error Code: " << WSAGetLastError() << endl;
					sendmut.unlock();
					break;
				}
				sendmut.unlock();
			}
			catch (exception& ex)
			{
				printf("-- TCPServer.Work Send Error\n%s\n", ex.what());
				break;
			}

			//힙 메모리 해제
			//curPos = 0;
			memmove(rcvbuf, rcvbuf + len + sizeof(int), curPos - len - sizeof(int));
			curPos -= len + sizeof(int);

			len = 0;
			lenDone = false;
			//printf("cont\n");
			continue;
		}
		//printf("_cont\n");
	}
	delete[] rcvbuf;
	delete[] tmpbuf;

	solver[solvernum]->stop();
	solver[solvernum]->reset(true);
	solverOccupied[solvernum] = false;
	printf("freed solver #%d\n", solvernum);
	Sleep(500);
	try
	{
		closesocket(clientSock);
		std::cout << "Client disconnected IP address = " << inet_ntoa(clientAddr.sin_addr) << ":" << ntohs(clientAddr.sin_port) << endl;
	}
	catch (exception& ex)
	{
		printf("-- TCPServer.Work closesocket Error\n%s\n", ex.what());
	}
	
	//reset();
	//allocated = false;
	for (auto ptr = clientList.begin(); ptr < clientList.end(); ptr++)
	{
		if ((*ptr)->get_id() == this_thread::get_id())
		{
			clientList.erase(ptr); //해당 클라이언트(쓰레드) 제거
			break;
		}
	}
}

void TCPServer::reset(int solvernum)
{
	mut.lock();
	//memset(vel, 0, solver[solvernum]->n3 * 3 * sizeof(float));
	//memset(tem, 0, solver[solvernum]->n3 * sizeof(float));
	//memset(solver[solvernum]->h_obstacle, 0, solver[solvernum]->n3b * sizeof(bool));
	mut.unlock();
}
