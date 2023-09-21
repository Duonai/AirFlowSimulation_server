#pragma once

#include <vector>
#include <stdio.h>
#include <mutex>

#include "device_launch_parameters.h"

#include "Setting.h"
//#include "Visualizer.h"
#define SOLVER_NUM 10



class Solver
{
public : 
    Solver(Setting *s);
    ~Solver();

    void init(int n);
    void init(int nx, int ny, int nz);
    bool isInitialized();
    void start();
    void stop();
    bool isPlay();
    void fastPlay();
    void getGridSize(size_t &x, size_t &y, size_t &z);

    void update();
    void addSource();
    void reset(bool occ = true);
    void allocate();

    void extendY();
    //void outletSwitch();
    //bool getOutletState();
    //void swingSwitch();
    //bool getSwingState();
    
    double getSimulationTime();
    int getFramestep();
    void setTimeStep(float ts);
    float getTimeStep();

    void setForce(float f);
    float getForce();
 
    void setObstacle(int x, int y, int z, bool val = true);
    void setObstacle(const bool* obstacle);
    void flipObstacle(int x, int y, int z);
    bool getObstacle(int x, int y, int z);
    bool getObstacleU(int x, int y, int z);
    bool getObstacleV(int x, int y, int z);
    bool getObstacleW(int x, int y, int z);

    void getObstacle(bool* o);
    //void getObstacleU(bool* o);
    //void getObstacleV(bool* o);
    //void getObstacleW(bool* o);

    void setACType(int acIdx, int i);
    void setACDirection(int acIdx, int dir);
    void setACWindType(int acIdx, WINDTYPE type);
    void setACPosition(int acIdx, int x, int y, int z);
    void setACFlowDirect(int acIdx, float x, float y, float z);
    void setACVentMag(int acIdx, int mag);
    void setACVentDir(int acIdx, int dir);

    // fluids param
    void setBuoyancyCoef(float coef);
    void setDiffusionCoef(float coef);
    void setViscosity(float coef);

    float getFPS();

    // graph data
    void setUpSamplingScale(int s);
    const std::vector<float*>& getTempDatas();
    void pushData();

    void generateObstacle();

    float getAvgTem();

 private:
    
    void hostAlloc();
    
    // for cuda
    void cudaAlloc();
    void cudaRelease();
    void cudaReset();

    void add_sourceGPU(float *dst, float* src, bool copy);
    void add_buoyancyGPU(float *dst, float* src);
    void advectGPU(float* dst, float* src, float* u, float* v, float* w, int boundaryType, bool useDirichlet = false, float* target = NULL);
    void diffuseGPU(float* dst, float* src, float coef, int boundaryType, bool useDirichlet = false, float* target = NULL);
    void linSolverGPU(float* dst, float* src, float a, bool sum, int boundaryType, bool useDirichlet = false, float* target = NULL);
    void projectGPU(float* u, float* v, float* w, float* p, float* div);

    void wuSourceGPU();
    void updateObstacleGPU();
    void velocityStepGPU();
    void densityStepGPU();

public :
    float getVelocityGPUU(int x, int y, int z);
    float getVelocityGPUV(int x, int y, int z);
    float getVelocityGPUW(int x, int y, int z);
    const float* getVelocityGPUN();
    void getVelocityGPUN(float* v);
    void getVelocityGPUN2(float* v);
    
    float getDensityAvgGPU();
    float getDensityGPU(int x, int y, int z);
    const float* getDensityGPUN();
    void getDensityGPUN(float* d);

public :

    bool initialized;
    bool simStart;

    size_t nx, ny, nz;
    size_t n3; 
    size_t n3b;
    int framestep;
    double totalTime;
    int swingStep;
    uint64_t tick = GetTickCount64();
    bool fpsflag = true;

    Setting* setting;
    bool resetFlag = false;
    bool DidReset();

    float timeStep;
    float diffuseCoef;   // for material
    float viscosityCoef; // for velocity
    float buoyancyCoef;
    float gravity;
    float force;
    float source;
    float sourceInc;

    bool addDensity;
    bool fastForward;

    bool occInitialized = false;
    float curTem;
    float tarTem;


    // ac variable
    //bool singleOutlet;
    //bool swingOutlet;

    int numAC;
    AirConditioner* airconditioners;

    // for graph data
    std::vector<float*> sampledTemps;
    int upsampleScale;

    float fps;

    char outputfile[128];
    int fileCnt;


    /// <summary>
    /// for cuda variable
    /// </summary>

    float* h_velocityU;
    float* h_velocityV;
    float* h_velocityW;
    float* h_velocityN;

    float* h_preVelocityU;
    float* h_preVelocityV;
    float* h_preVelocityW; 

    float* h_density;
    float* h_densityN;
    float* h_preDensity;

    bool* h_obstacle;
    void SetHalfBoundary();
    void ResetHalfBoundary();
    //bool* h_obstacleU;
    //bool* h_obstacleV;
    //bool* h_obstacleW;

    float* d_velocityU;
    float* d_velocityV;
    float* d_velocityW;

    //float* d_sourceVelocityU;
    //float* d_sourceVelocityV;
    //float* d_sourceVelocityW;

    float* d_preVelocityU;
    float* d_preVelocityV;
    float* d_preVelocityW;

    float* d_density;
    float* d_preDensity;

    bool* d_obstacle;
    //bool* d_obstacleU;
    //bool* d_obstacleV;
    //bool* d_obstacleW;

    // GPU
    int3 dim;
    int memSize;
    dim3 block;
    dim3 grid;

    //std::chrono::time_point<std::chrono::system_clock> start;

    /*static*/ std::mutex mut;
    //static std::mutex kernelMut;

};