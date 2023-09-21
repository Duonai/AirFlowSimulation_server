#pragma once

#include <iostream>
#include <chrono>
#include <iomanip>
#include <ctime> 
#include <thread>
#include <sstream>
#include <cstdlib>
#include <ctime>
//#include <string>

#include "freeglut.h"

//#include "solver.cuh"
#include "TCPServer.h"
//#include "Setting.h"

using std::this_thread::sleep_for;

using namespace std::chrono;
using namespace std;

//
//  main.c
//  Stable Fluids 3D
//


class Visualizer {


public:
    int windowWidth;
    int windowHeight;
    int windowOffsetW;
    int windowOffsetH;

    GLfloat translationX;
    GLfloat translationY;
    GLfloat translationZ;

    GLfloat rotationX;
    GLfloat rotationY;

    Setting* setting;

    size_t nx = 60;
    size_t ny = 27;
    size_t nz = 60;
    size_t maxN;

    float minTem = 12;
    float maxTem = 35;
    float srcTem = 15;
    float curTem = 30;

    //bool drawSecond = false;
    //int saveInterval = 10;
    //bool temHistory = false;
    //float endTime = 0;
    //bool saveSimData = false;


    Solver* solver;
    //TCPServer server;
    bool pause;
    bool pauseflag;
    bool debug;


    GLfloat alpha = 0.02f;
    int drawMode;
    bool draw;
    bool drawPerspective;
    bool drawAxis;
    bool drawObstacle;

    int logFrameStep;

    float total;
    float oldTotal;
    char filename[100] = { 0, };

    float drawScale;
    int drawThick;
    float drawVolStart;

    float preAvgTemp;

    int px, py, pz;
    bool simlog;
    char simlogname[512];
    uint64_t tick = GetTickCount64();

    int drawY;

    int room;

    const unsigned int lineLength = 30;
    const unsigned int lifespan = 2500;
    int visible = 0;
    struct ParticleTrace {
        float particlePos[100][3];
        float originalPos[3];
        int traceLength = 0;    //current length
        int age = 0;

    };
    const unsigned int numParticles = 1000;
    ParticleTrace particles[1000];

    float particleSpeed;

    float temAvg;
    int ventDir;

    bool realTime = true;
    const bool FPP = false;
    int prevframe = -1; //wuIdle

public:
    Visualizer();

    ~Visualizer();

    void initParticles();

    string getTime();

    float norm2tem(float t);

    float getLocalTem();

    float getAvgTem();

    void saveData(char* fn);

    void writeParam(FILE* f);

    bool OpenFile(char* fn);

    bool AppendData(char* fn, float time, float scale);

    void help();

    void wuInitialize();

    void wuDrawGrid();

    void wuTemptoGLColor(float temp, float opacity);

    float wuVorticity(float px, float py, float pz);
    bool wuInVortex(float px, float py, float pz);
    //float Dist(float* v0, float* v1);
    //float Curv(float* u, float* v, float* w);

    float wuGetDensity(float px, float py, float pz);

    void wuGetVelocity(float px, float py, float pz, float& vx, float& vy, float& vz);
    void wuGetVelocityRK4(float px, float py, float pz, float& vx, float& vy, float& vz);
    void wuGetNextStreamPosition(float px, float py, float pz, float& vx, float& vy, float& vz);
    
    bool wu_no_stick(float& x, float& y, float& z);

    void wuDrawDensity();
    void DrawObstacle();
    void wuDrawVelocity();

    void advectParticles();

    void advectParticlesRK4();

    float areaQuad(float* v0, float* v1, float* v2, float* v3);

    void wuDrawStreamline(bool animated);
    void wuDrawStreamsurface();
    void wuDrawStreamsurfaceWithSplit();
    void wuDrawStreamsurfaceWithLines();

    void wuDrawStreakline();
    void wuDrawStreaksurface();

    void wuDisplay();

    void wuReshape(GLint offsetW, GLint offsetH, GLint width, GLint height);

    void TemperatureDistribution(const float* t);

    void wuSpecialKeys(int key, int x, int y);

    void wuKeyBoard(unsigned char key, int x, int y);

    void wuIdle(void* a, int simNo);

    void setSettingVals(Setting* s);

};
