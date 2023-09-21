#pragma once

//#define IX(i,j,k) ((i)+(nx+2)*(j) + (nx+2)*(ny+2)*(k))
#define IX(i,j,k) ((ny+2)*(nz+2)*(i)+(nz+2)*(j) + (k))
#define SWAP(x0,x) {float * tmp=x0;x0=x;x=tmp;}
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define FOR_IJK for(int i=0; i<=nx; i++) { for(int j=0; j<=ny; j++) { for(int k=1; k<=nz; k++) {
#define FOR_XYZ for(int x=1; x<=nx; x++) { for(int y=1; y<=ny; y++) { for(int z=1; z<=nz; z++) {
#define FOR_XYZ0 for(int x=0; x<=nx+1; x++) { for(int y=0; y<=ny+1; y++) { for(int z=0; z<=nz+1; z++) {
#define FOR_END }}}

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <chrono>
//Cuda
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>

#include "solver.cuh"
#include "fluid_kernels.cuh"



Solver::Solver(Setting* s)
{
    initialized = false;
    simStart = false;
    fastForward = false;
    swingStep = 0;
    totalTime = 0;

    const EulerianSetting& eu = s->eulerianSetting;
    
    setting = s;

    timeStep = eu.timeStep;
    diffuseCoef = eu.diffuseCoef;       // material
    viscosityCoef = eu.viscosityCoef;      // velosity
    buoyancyCoef = eu.buoyancyCoef;
    gravity = eu.gravity;
    source = eu.source;
    force = eu.force;
    sourceInc = 0.f;
    
    addDensity = false;
    //singleOutlet = false;
    //swingOutlet = false;
    resetFlag = false;
    fps = 1.f;

    // AC variable
    numAC = s->machineCount;
    airconditioners = new AirConditioner[numAC];
    //printf("%d AC", numAC);

    occInitialized = false;
    curTem = s->eulerianSetting.curTem;
    


    h_velocityU = h_velocityV = h_velocityW = h_velocityN = NULL;
    h_preVelocityU = h_preVelocityV = h_preVelocityW = NULL;
    h_density = h_densityN = h_preDensity = NULL;
    h_obstacle = NULL;

    d_velocityU = d_velocityV = d_velocityW = NULL;
    d_preVelocityU = d_preVelocityV = d_preVelocityW = NULL;
    d_density = d_preDensity = NULL;
    d_obstacle = NULL;

    //h_obstacleU = h_obstacleV = h_obstacleW = d_obstacleU = d_obstacleV = d_obstacleW = NULL;
    //d_sourceVelocityU = d_sourceVelocityV = d_sourceVelocityW = NULL;

    init(eu.gridX, eu.gridY, eu.gridZ);
    upsampleScale = 1;

    for (int i = 0; i < numAC; i++) {
        const MachineSetting& mac = s->machineSetting[i];
        AirConditioner& ac = airconditioners[i];

        ac.acType = (ACTYPE)(mac.type);

        setACVentDir(i, mac.ventDir);
        ac.ventSpeed = mac.ventSpeed;
        ac.tarTemp = mac.tarTem;
        ac.windType = mac.windType;
        setACPosition(i, mac.posX, mac.posY, mac.posZ);
        setACDirection(i, mac.direction);
    }
}
 
Solver::~Solver()
{
    for (auto e : sampledTemps)
        delete[] e;
    sampledTemps.clear();

    delete  h_obstacle;

    cudaRelease();
}

void Solver::init(int n)
{
    init(n, n, n);
}

void Solver::init(int nx, int ny, int nz)
{
    initialized = false;
    simStart = false;

    printf("\n-- Start solver initialize\n");

    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    n3 = nx * ny * nz;
    n3b = (nx + 2) * (ny + 2) * (nz + 2);

    // for GPU
    dim = make_int3(nx + 2, ny + 2, nz + 2);
    memSize = n3b * sizeof(float);

    // Setup execution parameters
    int threads = 16;
    int dualThreads = 1024 / threads / threads;
    if (nx > nz)
    {
        block = dim3(dualThreads, threads, threads);
        grid = dim3((dim.x + dualThreads - 1) / dualThreads,
            (dim.y + threads - 1) / threads,
            (dim.z + threads - 1) / threads);
    }
    else
    {
        block = dim3(threads, threads, dualThreads);
        grid = dim3((dim.x + threads - 1) / threads,
            (dim.y + threads - 1) / threads,
            (dim.z + dualThreads - 1) / dualThreads);
    }

    allocate();

    printf("\nSimulation is initialized with size : %d x %d x %d\n", nx, ny, nz);

    initialized = true;
}

bool Solver::isInitialized()
{
    return initialized;
}

void Solver::start()
{
    const std::lock_guard<std::mutex> lock(mut);
    simStart = true;
    fpsflag = true;
    framestep = 0;
}

void Solver::stop()
{
    const std::lock_guard<std::mutex> lock(mut);
    simStart = false;
}

bool Solver::isPlay()
{
    return simStart;
}

void Solver::fastPlay()
{
    auto scale = setting->envSetting.fastScale;

    fastForward = !fastForward;
    if (fastForward)
    {
        timeStep *= scale;
        force *= scale;
        printf("\n-- accelerating simulation : %3.2f\n", scale);
    }
    else
    {
        timeStep /= scale;
        force /= scale;
        printf("\n-- restore simulation speed\n");
    }
}

void Solver::getGridSize(size_t& x, size_t& y, size_t& z)
{
    x = nx;
    y = ny;
    z = nz;
}

std::mutex kernelMut;


void Solver::update()
{
    if (!initialized)
        return;
    if (!simStart)
        return;


    // air source rnage is 8~12 degree.
    //float avgTem = getDentiyAvgGPU();
    //sourceInc = setting->envSetting.sAttenuationR * pow(setting->envSetting.sAttenuationP, 1 * (avgTem - 1.f));

    //auto start = std::chrono::system_clock::now();
    framestep++;
    totalTime += timeStep;

    wuSourceGPU();
    const std::lock_guard<std::mutex> lock(kernelMut);
    const std::lock_guard<std::mutex> lock2(mut);
    updateObstacleGPU();
    densityStepGPU();
    velocityStepGPU();


    //checkCudaErrors(cudaMemcpyAsync(d_preVelocityU, h_preVelocityU, memSize, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaMemcpyAsync(d_preVelocityV, h_preVelocityV, memSize, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaMemcpyAsync(d_preVelocityW, h_preVelocityW, memSize, cudaMemcpyHostToDevice));
    //add_sourceGPU(d_velocityU, d_preVelocityU, true);
    //add_sourceGPU(d_velocityV, d_preVelocityV, true);
    //add_sourceGPU(d_velocityW, d_preVelocityW, true);
    //checkCudaErrors(cudaMemcpyAsync(h_velocityU, d_velocityU, memSize, cudaMemcpyDeviceToHost));
    //checkCudaErrors(cudaMemcpyAsync(h_velocityV, d_velocityV, memSize, cudaMemcpyDeviceToHost));
    //checkCudaErrors(cudaMemcpyAsync(h_velocityW, d_velocityW, memSize, cudaMemcpyDeviceToHost));


    if (fpsflag) {
        tick = GetTickCount64();
        fpsflag = false;
    }
    fps = framestep * 1000.0f / (GetTickCount64() - tick + 1);
}

void Solver::addSource()
{
    addDensity = true;
}

void Solver::reset(bool occ)
{
    const std::lock_guard<std::mutex> lock(kernelMut);
    const std::lock_guard<std::mutex> lock2(mut);
    framestep = 0;
    fpsflag = true;
    swingStep = 0;
    totalTime = 0;
    memSize = n3b * sizeof(float);

    if (occ)
        memset(h_obstacle, 0, n3b * sizeof(bool));

    //for (int i = 10; i <= 50; i++)
    //    for (int j = 10; j <= 50; j++)
    //        for (int y = 2; y <= 20; y++) {
    //            if (i >= 28 && i <= 32 && j >= 28 && j <= 32)
    //                continue;
    //            if ((i * 7 + j * 5 + y * 3) % 31 == 1)
    //                setObstacle(i, y, j);
    //        }
    //printf("OCCRANDOM\n");
    //generateObstacle();

    for (auto e : sampledTemps)
        delete[] e;
    sampledTemps.clear();

    cudaReset();

    resetFlag = true;
}

void Solver::allocate()
{
    mut.lock();
    hostAlloc();
    cudaAlloc();
    mut.unlock();
    reset();
}

void Solver::extendY()
{

}

//void Solver::outletSwitch()
//{
//    singleOutlet = !singleOutlet;
//}

//bool Solver::getOutletState()
//{
//    return singleOutlet;
//}

//void Solver::swingSwitch()
//{
//    swingOutlet = !swingOutlet;
//}

//bool Solver::getSwingState()
//{
//    return swingOutlet;
//}

double Solver::getSimulationTime()
{
    return totalTime; //timeStep * framestep;
}

int Solver::getFramestep()
{
    return framestep;
}

void Solver::setTimeStep(float ts)
{
    timeStep = ts;
}

float Solver::getTimeStep()
{
    return timeStep;
}

void Solver::setForce(float f)
{
    force = f;
}

float Solver::getForce()
{
    return force;
}

void Solver::setObstacle(int x, int y, int z, bool val)
{
    //const std::lock_guard<std::mutex> lock(mut);
    if (h_obstacle) h_obstacle[IX(x, y, z)] = val;
    //if (val) {
    //    h_obstacleU[IX(x, y, z)] = true;
    //    h_obstacleV[IX(x, y, z)] = true;
    //    h_obstacleW[IX(x, y, z)] = true;

    //    if (x < dim.x - 1)
    //        h_obstacleU[IX(x + 1, y, z)] = true;
    //    if (y < dim.y - 1)
    //        h_obstacleV[IX(x, y + 1, z)] = true;
    //    if (z < dim.z - 1)
    //        h_obstacleW[IX(x, y, z + 1)] = true;
    //}
    //else {
    //    if (x > 0 && !h_obstacle[IX(x - 1, y, z)])
    //        h_obstacleU[IX(x, y, z)] = false;
    //    if (y > 0 && !h_obstacle[IX(x, y - 1, z)])
    //        h_obstacleV[IX(x, y, z)] = false;
    //    if (z > 0 && !h_obstacle[IX(x, y, z - 1)])
    //        h_obstacleW[IX(x, y, z)] = false;

    //    if (x < dim.x - 1 && !h_obstacle[IX(x + 1, y, z)])
    //        h_obstacleU[IX(x + 1, y, z)] = false;
    //    if (y < dim.y - 1 && !h_obstacle[IX(x, y + 1, z)])
    //        h_obstacleV[IX(x, y + 1, z)] = false;
    //    if (z < dim.z - 1 && !h_obstacle[IX(x, y, z + 1)])
    //        h_obstacleW[IX(x, y, z + 1)] = false;
    //}
}

void Solver::flipObstacle(int x, int y, int z)
{
    setObstacle(x, y, z, !h_obstacle[IX(x, y, z)]);
    //const std::lock_guard<std::mutex> lock(mut);
    //h_obstacle[IX(x, y, z)] = !h_obstacle[IX(x, y, z)];
}

//outdated
void Solver::setObstacle(const bool* occ)
{
    const std::lock_guard<std::mutex> lock(mut);
    int c = 0;
    cudaMemset(h_obstacle, 0, n3b * sizeof(bool));
    //cudaMemset(h_obstacleU, 0, n3b * sizeof(bool));
    //cudaMemset(h_obstacleV, 0, n3b * sizeof(bool));
    //cudaMemset(h_obstacleW, 0, n3b * sizeof(bool));
    try
    {
        FOR_XYZ
            bool val = occ[c++];
            if (val) {
                h_obstacle[IX(x, y, z)] = val;
                //h_obstacleU[IX(x, y, z)] = true;
                //h_obstacleV[IX(x, y, z)] = true;
                //h_obstacleW[IX(x, y, z)] = true;

                //if (x < dim.x - 1)
                //    h_obstacleU[IX(x + 1, y, z)] = true;
                //if (y < dim.y - 1)
                //    h_obstacleV[IX(x, y + 1, z)] = true;
                //if (z < dim.z - 1)
                //    h_obstacleW[IX(x, y, z + 1)] = true;
                
            }
        FOR_END;
        //for (int i = 0; i < numAC; i++) {
        //    AirConditioner& ac = airconditioners[i];
        //    //setObstacle(ac.pos[0], ac.pos[1], ac.pos[2], false);
        //}
    }
    catch (std::exception& ex)
    {
        printf("-- Solver.setObstacle Error\n%s\n", ex.what());
    }

    //SetHalfBoundary();
}

bool Solver::getObstacle(int x, int y, int z)
{
    const std::lock_guard<std::mutex> lock(mut);
    return h_obstacle[IX(x, y, z)];
}

//bool Solver::getObstacleU(int x, int y, int z)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    return h_obstacleU[IX(x, y, z)];
//}
//bool Solver::getObstacleV(int x, int y, int z)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    return h_obstacleV[IX(x, y, z)];
//}
//bool Solver::getObstacleW(int x, int y, int z)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    return h_obstacleW[IX(x, y, z)];
//}
void Solver::getObstacle(bool* o)
{
    const std::lock_guard<std::mutex> lock(mut);
    int c = 0;

#pragma omp parallel for
    FOR_XYZ0
        c = IX(x, y, z);
        o[c] = h_obstacle[c];
    FOR_END
}
//void Solver::getObstacleU(bool* o)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    int c = 0;
//
//#pragma omp parallel for
//    FOR_XYZ0
//        c = IX(x, y, z);
//        o[c] = h_obstacleU[c];
//    FOR_END
//}
//void Solver::getObstacleV(bool* o)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    int c = 0;
//
//#pragma omp parallel for
//    FOR_XYZ0
//        c = IX(x, y, z);
//        o[c] = h_obstacleV[c];
//    FOR_END
//}
//void Solver::getObstacleW(bool* o)
//{
//    const std::lock_guard<std::mutex> lock(mut);
//    int c = 0;
//
//#pragma omp parallel for
//    FOR_XYZ0
//        c = IX(x, y, z);
//        o[c] = h_obstacleW[c];
//    FOR_END
//}

void Solver::setACType(int acIdx, int i)
{
    AirConditioner& ac = airconditioners[acIdx];

    ac.acType = ACTYPE(i);
    //printf("\nset actype :  %d\n", i);
}

void Solver::setACDirection(int acIdx, int dir)
{
    AirConditioner& ac = airconditioners[acIdx];

    ac.direction = dir;

    ResetHalfBoundary();
    SetHalfBoundary();
}

void Solver::setACPosition(int acIdx, int x, int y, int z)
{
    AirConditioner& ac = airconditioners[acIdx];

    ac.pos[0] = x;
    ac.pos[1] = y;
    ac.pos[2] = z;

    //ResetHalfBoundary();
}

bool Solver::DidReset() {
    if (!resetFlag)
        return false;
    resetFlag = false;
    return true;
}

void Solver::setACFlowDirect(int acIdx, float x, float y, float z)
{
    //AirConditioner& ac = airconditioners[acIdx];

    //acFlowDirX = x;
    //acFlowDirY = y;
    //acFlowDirZ = z;
    assert(0);
}

void Solver::setACWindType(int acIdx, WINDTYPE type)
{
    AirConditioner& ac = airconditioners[acIdx];
    ac.windType = type;
}

void Solver::setACVentMag(int acIdx, int mag)
{
    AirConditioner& ac = airconditioners[acIdx];

    if ( mag > 0 && mag <= 10)
        ac.ventSpeed = mag;
}

void Solver::setACVentDir(int acIdx, int dir)
{
    AirConditioner& ac = airconditioners[acIdx];
    if (dir < 0)
        dir = 0;
    if (dir > 90)
        dir = 90;
    ac.ventLevel = dir;

    // 1,2,3, // swing
    //if ( dir > 3 )
    //{
    //    swingOutlet = true;
    //    ac.ventLevel = 4;
    //}
    //else if(dir > 0 &&   dir < 4 )
    //{
    //    ac.ventLevel = dir;
    //    swingOutlet = false;
    //}

    swingStep = 0;
}

void Solver::setBuoyancyCoef(float coef)
{
    buoyancyCoef = coef;
}

void Solver::setDiffusionCoef(float coef)
{
    diffuseCoef = coef;
}

void Solver::setViscosity(float coef)
{
    viscosityCoef = coef;
}

float Solver::getFPS()
{
    return fps;
}

void Solver::setUpSamplingScale(int s)
{
    upsampleScale = s;
}

const std::vector<float*>& Solver::getTempDatas()
{
    return sampledTemps;
}

void Solver::pushData()
{
    float* d = new float[n3 /pow(upsampleScale, 3)];
    //memcpy(d, h_densityN, n3 * sizeof(float));

    int idx = 0;
    for (int i = 0; i<nx; i+= upsampleScale)
        for (int j = 0; j < ny; j += upsampleScale)
            for (int k = 0; k < nz; k += upsampleScale)
                d[idx++] = h_density[IX(i, j, k)];

    sampledTemps.push_back(d);
    return;

}

void Solver::hostAlloc()
{

}


void Solver::generateObstacle()
{
    int s = 5;
    for (int i = -s; i < s; i++)
        for (int j = -s; j < s; j++)
            for (int k = -s; k < s; k++)
            {
                setObstacle(nx/2+i, ny/2+j, nz/2+k);
            }

    for (int xx = 1; xx <= nx; xx++)
        for (int yy = 1; yy <= nz/2; yy++)
        {
            setObstacle(xx, ny/2, yy);
        }
}

void Solver::cudaAlloc()
{
    cudaRelease();
    checkCudaErrors(cudaMallocHost((void**)&h_velocityU, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_velocityV, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_velocityW, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_velocityN, n3 * sizeof(float) * 3));

    checkCudaErrors(cudaMallocHost((void**)&h_preVelocityU, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_preVelocityV, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_preVelocityW, memSize));

    checkCudaErrors(cudaMallocHost((void**)&h_density, memSize));
    checkCudaErrors(cudaMallocHost((void**)&h_densityN, n3 * sizeof(float) ));
    checkCudaErrors(cudaMallocHost((void**)&h_preDensity, memSize));

    int sizeBool = n3b * sizeof(bool);
    checkCudaErrors(cudaMallocHost((void**)&h_obstacle, sizeBool));
    //checkCudaErrors(cudaMallocHost((void**)&h_obstacleU, sizeBool));
    //checkCudaErrors(cudaMallocHost((void**)&h_obstacleV, sizeBool));
    //checkCudaErrors(cudaMallocHost((void**)&h_obstacleW, sizeBool));

    // Allocate device memory
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_velocityU), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_velocityV), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_velocityW), memSize));

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_preVelocityU), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_preVelocityV), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_preVelocityW), memSize));

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_density), memSize));
    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_preDensity), memSize));

    checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_obstacle), sizeBool));
    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_obstacleU), sizeBool));
    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_obstacleV), sizeBool));
    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_obstacleW), sizeBool));

    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_sourceVelocityU), memSize));
    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_sourceVelocityV), memSize));
    //checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_sourceVelocityW), memSize));
}

void Solver::cudaRelease()
{
    checkCudaErrors(cudaFreeHost(h_velocityU));
    checkCudaErrors(cudaFreeHost(h_velocityV));
    checkCudaErrors(cudaFreeHost(h_velocityW));
    checkCudaErrors(cudaFreeHost(h_velocityN));
    checkCudaErrors(cudaFreeHost(h_preVelocityU));
    checkCudaErrors(cudaFreeHost(h_preVelocityV));
    checkCudaErrors(cudaFreeHost(h_preVelocityW));
    checkCudaErrors(cudaFreeHost(h_density));
    checkCudaErrors(cudaFreeHost(h_densityN));
    checkCudaErrors(cudaFreeHost(h_preDensity));

    checkCudaErrors(cudaFreeHost(h_obstacle));
    //checkCudaErrors(cudaFreeHost(h_obstacleU));
    //checkCudaErrors(cudaFreeHost(h_obstacleV));
    //checkCudaErrors(cudaFreeHost(h_obstacleW));

    checkCudaErrors(cudaFree(d_velocityU));
    checkCudaErrors(cudaFree(d_velocityV));
    checkCudaErrors(cudaFree(d_velocityW));
    checkCudaErrors(cudaFree(d_preVelocityU));
    checkCudaErrors(cudaFree(d_preVelocityV));
    checkCudaErrors(cudaFree(d_preVelocityW));
    checkCudaErrors(cudaFree(d_density));
    checkCudaErrors(cudaFree(d_preDensity));
    checkCudaErrors(cudaFree(d_obstacle));

    //checkCudaErrors(cudaFree(d_obstacleU));
    //checkCudaErrors(cudaFree(d_obstacleV));
    //checkCudaErrors(cudaFree(d_obstacleW));
    //checkCudaErrors(cudaFree(d_sourceVelocityU));
    //checkCudaErrors(cudaFree(d_sourceVelocityV));
    //checkCudaErrors(cudaFree(d_sourceVelocityW));
}

void Solver::cudaReset()
{
    cudaMemset(h_velocityU, 0, memSize);
    cudaMemset(h_velocityV, 0, memSize);
    cudaMemset(h_velocityW, 0, memSize);
    cudaMemset(h_velocityN, 0, nx * ny * nz * sizeof(float) * 3);

    checkCudaErrors(cudaMemcpy(d_velocityU, h_velocityU, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_velocityV, h_velocityV, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_velocityW, h_velocityW, memSize, cudaMemcpyHostToDevice));

    cudaMemset(h_preVelocityU, 0, memSize);
    cudaMemset(h_preVelocityV, 0, memSize);
    cudaMemset(h_preVelocityW, 0, memSize);

    checkCudaErrors(cudaMemcpy(d_preVelocityU, h_preVelocityU, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_preVelocityV, h_preVelocityV, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_preVelocityW, h_preVelocityW, memSize, cudaMemcpyHostToDevice));

    //checkCudaErrors(cudaMemcpy(d_sourceVelocityU, d_preVelocityU, memSize, cudaMemcpyDeviceToDevice));
    //checkCudaErrors(cudaMemcpy(d_sourceVelocityV, d_preVelocityV, memSize, cudaMemcpyDeviceToDevice));
    //checkCudaErrors(cudaMemcpy(d_sourceVelocityW, d_preVelocityW, memSize, cudaMemcpyDeviceToDevice));

    //memset(h_density, 0, memSize);
    for (int i = 0; i < n3b; i++) {
        h_density[i] = setting->eulerianSetting.curTem;
        h_preDensity[i] = setting->eulerianSetting.curTem;
    }
    //cudaMemset(h_preDensity, setting->eulerianSetting.rstTem, memSize);
    checkCudaErrors(cudaMemcpy(d_density, h_density, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_preDensity, h_preDensity, memSize, cudaMemcpyHostToDevice));
}

void Solver::add_sourceGPU(float* dst, float* src, bool copy)
{
    addSourceKernel <<< grid, block >>> (dst, src, dim, timeStep, copy);
}

void Solver::add_buoyancyGPU(float* dst, float* src)
{
    addBuoyancyKernel <<< grid, block >>> (dst, src, dim, timeStep, 0.f, buoyancyCoef, gravity, d_obstacle);
}

//TD
void Solver::advectGPU(float* dst, float* src, float* u, float* v, float* w, int boundaryType, bool useDirichlet, float* target)
{
    advectionKernel <<< grid, block >>> (dst, src, u, v, w, dim, timeStep);

    if (useDirichlet)
        setBoundaryKernelDir <<< grid, block >>> (dst, d_obstacle, dim, boundaryType, curTem);
        //setBoundaryKernelDir <<< grid, block >>> (dst, d_obstacle, dim, boundaryType, target);
    else
        setBoundaryKernel <<< grid, block >>> (dst, d_obstacle, dim, boundaryType);
}

void Solver::diffuseGPU(float* dst, float* src, float coef, int boundaryType, bool useDirichlet, float* target)
{
    float a = timeStep * coef * n3;
    linSolverGPU(dst, src, -a, true, boundaryType, useDirichlet, target);
}

float Solver::getAvgTem() {
    double temp = 0.f;

    // for 2nd test room
    float a = -4.78571f;
    float b = 72.78571f;

    float a2 = 6.0909f;
    float b2 = -413.182f;

    //float temp = 0.f;
    int aaa = 0;

    for (int x = 1; x < nx; x++)
    {
        //int minz = (x - b) / a;
        //int maxz = (x - b2) / a2;
        for (int z = 1; z < nz; z++)
        {
            //if (minz < z && z < maxz)
            {
                for (int y = 1; y < ny; y++)
                {
                    temp += getDensityGPU(x, y, z);
                    aaa++;
                }
            }
        }
    }

    return (temp / aaa);
}

//TD
void Solver::linSolverGPU(float* dst, float* src, float a, bool sum, int boundaryType, bool useDirichlet, float* target)
{
    static float* d_temp = NULL;
    if(d_temp == NULL)
        checkCudaErrors(cudaMalloc(reinterpret_cast<void**>(&d_temp), memSize));

    for (int k = 0; k < 20; k++)
    {
        lin_solveKernel <<< grid, block >>> (d_temp, dst, src, dim, a, sum, d_obstacle);

        std::swap(dst, d_temp);

        if (useDirichlet)
            setBoundaryKernelDir <<< grid, block >>> (dst, d_obstacle, dim, boundaryType, curTem);
        else
            setBoundaryKernel <<< grid, block >>> (dst, d_obstacle, dim, boundaryType);
    }
}

//TD
void Solver::projectGPU(float* u, float* v, float* w, float* p, float* div)
{
    cudaMemsetAsync(p, 0, memSize);
    cudaMemsetAsync(div, 0, memSize);

    divergenceKernel <<< grid, block >>> (div, u, v, w, dim, d_obstacle);

    setBoundaryKernel <<< grid, block >>> (div, d_obstacle, dim, 0);
    setBoundaryKernel <<< grid, block >>> (p, d_obstacle, dim, 0);

    linSolverGPU(p, div, 1, false, 0);

    applyPressureKernel <<< grid, block >>> (u, v, w, p, dim, d_obstacle);

    setBoundaryKernel <<<grid, block >>> (u, d_obstacle, dim, 1);
    setBoundaryKernel <<<grid, block >>> (v, d_obstacle, dim, 2);
    setBoundaryKernel <<<grid, block >>> (w, d_obstacle, dim, 3);
}

static int wuSourceGPUCount = 0;
void Solver::wuSourceGPU()
{
    float* density = h_preDensity;
    float* u = h_preVelocityU;
    float* v = h_preVelocityV;
    float* w = h_preVelocityW;
    auto& setmc = setting->machineSetting;

    const int size = n3b * sizeof(float);

    cudaMemset(u, 0, size);
    cudaMemset(v, 0, size);
    cudaMemset(w, 0, size);
    cudaMemset(density, 0, size);


    if (addDensity == true)
    {
        for (int m_i = 0; m_i < numAC; m_i++) {
            AirConditioner& ac = airconditioners[m_i];
            float outV = force * (ac.ventSpeed * 0.05f + 0.4f);
            if (ac.acType == ACTYPE::Wall)
            {
                float swing = ac.ventLevel; //20.f + (10.f * (ac.ventLevel == 4 ? 2 : (ac.ventLevel - 1))) + (ac.ventLevel == 4 ? (10 * sin(framestep * setmc[m_i].swingSpeed)) : 0.f);
                float angle = swing * 3.141592f / 180.f;
                float angleH = outV * cos(angle);
                float angleV = outV * sin(angle);

                int widthSide = 3;
                int heightSide = 1;

                if (wuSourceGPUCount++ == 0)
                {
                    if(ac.pos[0] > 4 && ac.pos[0] < nx - 5 && ac.pos[2] > 4 && ac.pos[2] < nz - 5)
                        printf("The AC position is not correct\n");
                }
                //if (ac.pos[0] < 5)
                //    ac.direction = 0;
                //else if (ac.pos[0] > nx - 5)
                //    ac.direction = 4;
                //else if (ac.pos[2] < 5)
                //    ac.direction = 6;
                //else if (ac.pos[2] > nz - 5)
                //    ac.direction = 2;
                //else {
                //    printf("The AC position is not correct\n");
                //    continue;
                //}

                float* nor = u;     // normal of ac
                float* tan = w;     // tangential of ac
                float* ver = v;
                int acW = 0;
                int acT = widthSide;
                int acH = heightSide;

                bool isACPosX = (ac.direction % 4 == 0);
                int reverse = (ac.direction == 0 || ac.direction == 6) ? 1 : -1;
                int offset = 0;// 2;


                if (ac.direction == 0 || ac.direction == 4) {
                    nor = u;
                    tan = w;
                    acW = 0;
                    acT = widthSide;
                }
                else if (ac.direction == 2 || ac.direction == 6) {
                    nor = w;
                    tan = u;
                    acW = widthSide;
                    acT = 0;
                }


                for (int i = -acW; i <= acW; i++)
                {
                    for (int j = -acH; j <= acH; j++)
                    {
                        for (int k = -acT; k <= acT; k++)
                        {
                            auto tangent = angleH * (isACPosX ? float(k) : float(i)) / float(widthSide * 3);
                            density[IX(ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k)] = source;
                            nor[IX(ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k)] = angleH * reverse;
                            tan[IX(ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k)] = tangent;
                            ver[IX(ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k)] = -angleV;

                            //printf("%d %d %d ==  %d %d %d \n", i, j, k, ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k);
                        }
                    }
                }


            }
            else if (ac.acType == ACTYPE::Stand)
            {
                float swing = ac.ventLevel - 45; //-30.f + (30.f * ((ac.ventLevel == 4 ? 2 : ac.ventLevel) - 1)) + (ac.ventLevel == 4 ? (30 * sin(framestep * setmc[m_i].swingSpeed)) : 0.f);
                float angle = swing * 3.141592f / 180.f;
                float angleH = outV * cos(angle);
                float angleV = outV * sin(angle);

                int widthSide = 2;
                int heightSide = 1;


                int x_wall_dist = std::min(ac.pos[0], (int)nx - ac.pos[0] - 1);
                int z_wall_dist = std::min(ac.pos[2], (int)nz - ac.pos[2] - 1);

                float3 ac_rel = make_float3(0.5f - (float)ac.pos[0] / nx, 0.5f - (float)ac.pos[1] / ny, 0.5f - (float)ac.pos[2] / nz);
                //bool isACPosX = abs(ac_rel.x) > abs(ac_rel.z) ? true : false;
                bool isACPosX = x_wall_dist < z_wall_dist ? true : false;

                int reverse = (isACPosX ? ac_rel.x : ac_rel.z) < 0 ? -1 : 1;

                if (wuSourceGPUCount++ == 0)
                {
                    std::cout << (reverse == -1 ? '-' : '+') << (isACPosX ? 'x' : 'z') << std::endl;
                    if (isACPosX && reverse)
                        ac.direction = 0;
                    else if (isACPosX && !reverse)
                        ac.direction = 4;
                    else if (!isACPosX && !reverse)
                        ac.direction = 2;
                    else if (!isACPosX && reverse)
                        ac.direction = 6;
                }

                float* nor = (isACPosX ? u : w);     // normal of ac
                float* tan = (isACPosX ? w : u);     // tangential of ac
                float* ver = v;
                int acW = (isACPosX ? 1 : widthSide);
                int acT = (isACPosX ? widthSide : 1);
                int acH = heightSide;

                for (int i = -acW; i <= acW; i++)
                {
                    for (int j = -acH; j <= acH; j++)
                    {
                        for (int k = -acT; k <= acT; k++)
                        {
                            int win_x = ac.pos[0] + i;
                            int win_y = ac.pos[1] + j;
                            int win_z = ac.pos[2] + k;

                            auto tangent = angleH * (isACPosX ? float(k) : float(i)) / float(widthSide * 3.f);
                            density[IX(win_x, win_y, win_z)] = source;
                            nor[IX(win_x, win_y, win_z)] = angleH * reverse;
                            tan[IX(win_x, win_y, win_z)] = tangent;
                            ver[IX(win_x, win_y, win_z)] = -angleV;

                            //printf("%d %d %d ==  %d %d %d \n", i, j, k, ac.pos[0] + i, ac.pos[1] + j, ac.pos[2] + k);
                        }
                    }
                }
            }
            else if (ac.acType == ACTYPE::Ceiling)
            {
                // 4way 32m^3/min -> a vein 8m^3/min -> 0.133m^3/sec
                // a vein has 16cell with 0.01m^2 face. -> 0.16m^2 
                // the velocity 0.83 m/s = 0.133 / 0.16

                //if ( framestep % set.envSetting.saveInterval == 0 )
                //    printf("-- source : %f / inc : %f / total : %f \n", source, inc, source-inc);

                //float swing = setmc.initAngle + (setmc.interval * ac.ventLevel)  + (swingOutlet ? ( 15 * sin(framestep * setmc.swingSpeed)) : 0.f) ;  // 1256 cycle  

                float swing = setmc[m_i].initAngle + (setmc[m_i].interval * (ac.ventLevel - 1));
                bool flipDirs = false;

                if (ac.ventLevel == 4)
                {
                    swing = ac.ventLevel; //setmc[m_i].initAngle + setmc[m_i].interval * (-cos(swingStep++ * setmc[m_i].swingSpeed));
                }
                if (ac.windType == WINDTYPE::AIRGUIDE)
                {
                    swing = 15;
                }
                else if (ac.windType == WINDTYPE::HIGHCEILING)
                {
                    //float outV = force * (3 * 0.05f + 0.4f);
                    swing = 85;
                    //printf("85");
                }
                else if (ac.windType == WINDTYPE::POWER)
                {
                    int localMode = ((int)(totalTime / 2.f) / 5) % 4;
                    if (localMode == 0) {
                        swing = 85;
                    }
                    else if (localMode == 2) {
                        swing = 15;
                    }
                    else {
                        float timeAt = remainder(totalTime / 2.f, 20);
                        if (timeAt < 0) timeAt += 20;

                        if (localMode == 3) {
                            timeAt = 20 - timeAt;
                        }
                        else
                            timeAt -= 5;

                        if (timeAt >= 0 && timeAt < 1)
                            swing = 85 - timeAt * 20;
                        else if (timeAt >= 4 && timeAt < 5)
                            swing = 115 - timeAt * 20;
                        else
                            swing = 50 + 15 * cos((timeAt - 1) * 3.141592);
                    }
                    flipDirs = true;
                }
                else if (ac.windType == WINDTYPE::AUTO)
                {
                    //printf("AvgTem: %f, tarTemp : %f\n", getAvgTem(), ac.tarTemp);
                    if (getAvgTem() < ac.tarTemp + 3) {
                        swing = 15;
                    }
                    else {
                        int localMode = ((int)(totalTime / 2.f) / 5) % 4;
                        if (localMode == 0) {
                            swing = 85;
                        }
                        else if (localMode == 2) {
                            swing = 15;
                        }
                        else {
                            float timeAt = remainder(totalTime / 2.f, 20);
                            if (timeAt < 0) timeAt += 20;

                            if (localMode == 3) {
                                timeAt = 20 - timeAt;
                            }
                            else
                                timeAt -= 5;

                            if (timeAt >= 0 && timeAt < 1)
                                swing = 85 - timeAt * 20;
                            else if (timeAt >= 4 && timeAt < 5)
                                swing = 115 - timeAt * 20;
                            else
                                swing = 50 + 15 * cos((timeAt - 1) * 3.141592);
                        }
                        flipDirs = true;
                    }
                }
                else if (ac.windType == WINDTYPE::FOREST)
                {
                    swing = 50 + 15 * sin(totalTime / 2.0 * 3.141592);
                    //swing = 50;
                }
                else if (ac.windType == WINDTYPE::MANUAL) {
                    //swing = 68;
                    outV = force * (ac.ventSpeed * 0.015f + 0.4f);
                    swing = 90 - ac.ventLevel;
                }



                float angle = swing * 3.141592f / 180.f;
                float angleH = outV * cos(angle);
                float angleV = outV * sin(angle);

                int offset = setmc[m_i].ventDist;
                int halfcell = offset / 2;
                //auto wsource = (source - sourceInc) * (ac.ventSpeed * 0.2f + 0.4f);
                auto wsource = source;
                auto acHeight = ac.pos[1] - 1;
                //float3 ac_rel = make_float3(0.5f - (float)ac.pos[0] / nx, 0.5f - (float)ac.pos[1] / ny, 0.5f - (float)ac.pos[2] / nz);
                //bool isACPosX = abs(ac_rel.x) > abs(ac_rel.z) ? true : false;
                //int reverse = (isACPosX ? ac_rel.x : ac_rel.z) < 0 ? 1 : -1;
                for (int j = -1; j < 2; j++)
                {
                    for (int i = -halfcell + abs(j); i <= halfcell - abs(j); i++)
                    {
                        auto tangent = float(i) / float(offset)/10.0f; // *  float(offset-2) / float(offset) ;
                        auto absRate = fabs(float(i) / float(offset));

                        //right
                        density[IX(ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i)] = wsource;
                        u[IX(ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i)] = angleH * (1.f - absRate);
                        v[IX(ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i)] = -angleV;
                        w[IX(ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i)] = tangent;
                        //printf("(%2d,%2d,%2d)\n", ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i);
                        // left
                        density[IX(ac.pos[0] - (offset + j), acHeight, ac.pos[2] + i)] = wsource;
                        u[IX(ac.pos[0] - (offset + j), acHeight, ac.pos[2] + i)] = -angleH * (1.f - absRate);
                        v[IX(ac.pos[0] - (offset + j), acHeight, ac.pos[2] + i)] = -angleV;
                        w[IX(ac.pos[0] - (offset + j), acHeight, ac.pos[2] + i)] = tangent;

                        if (flipDirs) {
                            float temp = angleH;
                            angleH = angleV;
                            angleV = temp;
                        }
                        // top
                        density[IX(ac.pos[0] + i, acHeight, ac.pos[2] + (offset + j))] = wsource;
                        w[IX(ac.pos[0] + i, acHeight, ac.pos[2] + (offset + j))] = angleH * (1.f - absRate);
                        v[IX(ac.pos[0] + i, acHeight, ac.pos[2] + (offset + j))] = -angleV;
                        u[IX(ac.pos[0] + i, acHeight, ac.pos[2] + (offset + j))] = tangent;

                        // bottom
                        density[IX(ac.pos[0] + i, acHeight, ac.pos[2] - (offset + j))] = wsource;
                        w[IX(ac.pos[0] + i, acHeight, ac.pos[2] - (offset + j))] = -angleH * (1.f - absRate);
                        v[IX(ac.pos[0] + i, acHeight, ac.pos[2] - (offset + j))] = -angleV;
                        u[IX(ac.pos[0] + i, acHeight, ac.pos[2] - (offset + j))] = tangent;

                        if (flipDirs) {
                            float temp = angleH;
                            angleH = angleV;
                            angleV = temp;
                        }
                    }
                }
            }
            else if (ac.acType == ACTYPE::Tower)
            {
                const int acHeight = 15 + ac.pos[1];
                float angleGuard_tr = ac.ventLevel;
                float angleGuard_tl = ac.ventLevel;
                float angleGuard_br = ac.ventLevel;
                float angleGuard_bl = ac.ventLevel;
                float angleCircle = 10;

                outV = force * 0.25f; //0.25 for scale
                //Cscale: circle wind speed, Gscale: guard wind speed
                float Cscale = (1.f + 0.1f * (ac.ventSpeed - 1) * 2 / 9.0f);
                float Gscale = (1.f + 0.25f * (ac.ventSpeed - 1) * 2 / 9.0f);

                float angleIn = ac.direction % 2 == 0 ? 22 : 12;

                /*if (ac.windType == WINDTYPE::WIDECARE)
                {
                    angleGuard_tr = 53;
                    angleGuard_tl = 53;
                    angleGuard_br = 53;
                    angleGuard_bl = 53;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::SPACEDIV)
                {
                    angleGuard_tr = angleIn;
                    angleGuard_tl = angleIn;
                    angleGuard_br = 53;
                    angleGuard_bl = 53;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::XXXX)
                {
                    angleGuard_tr = angleIn;
                    angleGuard_tl = angleIn;
                    angleGuard_br = angleIn;
                    angleGuard_bl = angleIn;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::LEFT)
                {
                    angleGuard_tr = angleIn;
                    angleGuard_tl = 53;
                    angleGuard_br = angleIn;
                    angleGuard_bl = 53;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::RIGHT)
                {
                    angleGuard_tr = 53;
                    angleGuard_tl = angleIn;
                    angleGuard_br = 53;
                    angleGuard_bl = angleIn;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::DIV_L)
                {
                    angleGuard_tr = 0;
                    angleGuard_tl = 53;
                    angleGuard_br = 0;
                    angleGuard_bl = 53;
                    angleCircle = 10;
                }
                else if (ac.windType == WINDTYPE::DIV_R)
                {
                    angleGuard_tr = 53;
                    angleGuard_tl = 0;
                    angleGuard_br = 53;
                    angleGuard_bl = 0;
                    angleCircle = 10;
                }*/

                //axis direction
                if (ac.direction % 2 == 0 && ac.direction < 8) {
                    bool isXdir = ac.direction % 4 == 0;
                    int flipped = isXdir ? (ac.direction == 0 ? 1 : -1) : (ac.direction == 6 ? 1 : -1);
                    
                    if (angleCircle != 0)
                    {   //circle
                        static const float COS8 = sqrt(2 + sqrt(2)) / 2;
                        static const float SIN8 = sqrt(2 - sqrt(2)) / 2;
                        float angleA = angleCircle * 3.141592f / 180.f;
                        float outVL = Cscale * outV * COS8 * sin(angleA);
                        float outVS = Cscale * outV * SIN8 * sin(angleA);
                        float outF = Cscale * outV * cos(angleA);
                        
                        if (isXdir) {
                            density[IX(ac.pos[0], acHeight, ac.pos[2] - flipped * 1)] = source;
                            u[IX(ac.pos[0], acHeight, ac.pos[2] - flipped * 1)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight, ac.pos[2] - flipped * 1)] = -outVL;
                            w[IX(ac.pos[0], acHeight, ac.pos[2] - flipped * 1)] = -flipped * outVS;

                            density[IX(ac.pos[0], acHeight, ac.pos[2])] = source;
                            u[IX(ac.pos[0], acHeight, ac.pos[2])] = flipped * outF;
                            v[IX(ac.pos[0], acHeight, ac.pos[2])] = -outVL;
                            w[IX(ac.pos[0], acHeight, ac.pos[2])] = flipped * outVS;


                            density[IX(ac.pos[0], acHeight + 1, ac.pos[2] - flipped * 2)] = source;
                            u[IX(ac.pos[0], acHeight + 1, ac.pos[2] - flipped * 2)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 1, ac.pos[2] - flipped * 2)] = -outVS;
                            w[IX(ac.pos[0], acHeight + 1, ac.pos[2] - flipped * 2)] = -flipped * outVL;

                            density[IX(ac.pos[0], acHeight + 1, ac.pos[2] + flipped * 1)] = source;
                            u[IX(ac.pos[0], acHeight + 1, ac.pos[2] + flipped * 1)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 1, ac.pos[2] + flipped * 1)] = -outVS;
                            w[IX(ac.pos[0], acHeight + 1, ac.pos[2] + flipped * 1)] = flipped * outVL;

                            density[IX(ac.pos[0], acHeight + 2, ac.pos[2] - flipped * 2)] = source;
                            u[IX(ac.pos[0], acHeight + 2, ac.pos[2] - flipped * 2)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 2, ac.pos[2] - flipped * 2)] = outVS;
                            w[IX(ac.pos[0], acHeight + 2, ac.pos[2] - flipped * 2)] = -flipped * outVL;

                            density[IX(ac.pos[0], acHeight + 2, ac.pos[2] + flipped * 1)] = source;
                            u[IX(ac.pos[0], acHeight + 2, ac.pos[2] + flipped * 1)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 2, ac.pos[2] + flipped * 1)] = outVS;
                            w[IX(ac.pos[0], acHeight + 2, ac.pos[2] + flipped * 1)] = flipped * outVL;


                            density[IX(ac.pos[0], acHeight + 3, ac.pos[2] - flipped * 1)] = source;
                            u[IX(ac.pos[0], acHeight + 3, ac.pos[2] - flipped * 1)] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 3, ac.pos[2] - flipped * 1)] = outVL;
                            w[IX(ac.pos[0], acHeight + 3, ac.pos[2] - flipped * 1)] = -flipped * outVS;

                            density[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = source;
                            u[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = flipped * outF;
                            v[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = outVL;
                            w[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = flipped * outVS;

                        }
                        else {
                            density[IX(ac.pos[0] + flipped * 1, acHeight, ac.pos[2])] = source;
                            u[IX(ac.pos[0] + flipped * 1, acHeight, ac.pos[2])] = flipped * outVS;
                            v[IX(ac.pos[0] + flipped * 1, acHeight, ac.pos[2])] = -outVL;
                            w[IX(ac.pos[0] + flipped * 1, acHeight, ac.pos[2])] = flipped * outF;

                            density[IX(ac.pos[0], acHeight, ac.pos[2])] = source;
                            u[IX(ac.pos[0], acHeight, ac.pos[2])] = -flipped * outVS;
                            v[IX(ac.pos[0], acHeight, ac.pos[2])] = -outVL;
                            w[IX(ac.pos[0], acHeight, ac.pos[2])] = flipped * outF;


                            density[IX(ac.pos[0] + flipped * 2, acHeight + 1, ac.pos[2])] = source;
                            u[IX(ac.pos[0] + flipped * 2, acHeight + 1, ac.pos[2])] = flipped * outVL;
                            v[IX(ac.pos[0] + flipped * 2, acHeight + 1, ac.pos[2])] = -outVS;
                            w[IX(ac.pos[0] + flipped * 2, acHeight + 1, ac.pos[2])] = flipped * outF;

                            density[IX(ac.pos[0] - flipped * 1, acHeight + 1, ac.pos[2])] = source;
                            u[IX(ac.pos[0] - flipped * 1, acHeight + 1, ac.pos[2])] = -flipped * outVL;
                            v[IX(ac.pos[0] - flipped * 1, acHeight + 1, ac.pos[2])] = -outVS;
                            w[IX(ac.pos[0] - flipped * 1, acHeight + 1, ac.pos[2])] = flipped * outF;

                            density[IX(ac.pos[0] + flipped * 2, acHeight + 2, ac.pos[2])] = source;
                            u[IX(ac.pos[0] + flipped * 2, acHeight + 2, ac.pos[2])] = flipped * outVL;
                            v[IX(ac.pos[0] + flipped * 2, acHeight + 2, ac.pos[2])] = outVS;
                            w[IX(ac.pos[0] + flipped * 2, acHeight + 2, ac.pos[2])] = flipped * outF;

                            density[IX(ac.pos[0] - flipped * 1, acHeight + 2, ac.pos[2])] = source;
                            u[IX(ac.pos[0] - flipped * 1, acHeight + 2, ac.pos[2])] = -flipped * outVL;
                            v[IX(ac.pos[0] - flipped * 1, acHeight + 2, ac.pos[2])] = outVS;
                            w[IX(ac.pos[0] - flipped * 1, acHeight + 2, ac.pos[2])] = flipped * outF;


                            density[IX(ac.pos[0] + flipped * 1, acHeight + 3, ac.pos[2])] = source;
                            u[IX(ac.pos[0] + flipped * 1, acHeight + 3, ac.pos[2])] = flipped * outVS;
                            v[IX(ac.pos[0] + flipped * 1, acHeight + 3, ac.pos[2])] = outVL;
                            w[IX(ac.pos[0] + flipped * 1, acHeight + 3, ac.pos[2])] = flipped * outF;

                            density[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = source;
                            u[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = -flipped * outVS;
                            v[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = outVL;
                            w[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = flipped * outF;

                        }
                    }   //end circle

                    {   //guard

                        //Gscale += 0.8f;

                        float angle_r = angleGuard_br * 3.141592f / 180.f;
                        float angle_l = angleGuard_bl * 3.141592f / 180.f;
                        float angleF_r = (angleGuard_br < 0) ? 0 : Gscale * outV * cos(angle_r);
                        float angleS_r = (angleGuard_br < 0) ? 0 : Gscale * outV * sin(angle_r);
                        float angleF_l = (angleGuard_bl < 0) ? 0 : Gscale * outV * cos(angle_l);
                        float angleS_l = (angleGuard_bl < 0) ? 0 : Gscale * outV * sin(angle_l);

                        if (isXdir) {
                            for (int i = 5; i <= 9; i++)
                            {
                                density[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = source;
                                u[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = flipped * angleF_r;
                                v[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = 0.000001f;
                                w[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = -flipped * angleS_r;
                                 //printf("- %d %d %d: %8f %8f %8f\n", ac.pos[0], i, ac.pos[2] - flipped * 2, u[IX(ac.pos[0], i, ac.pos[2] - flipped * 2)], v[IX(ac.pos[0], i, ac.pos[2] - flipped * 2)], w[IX(ac.pos[0], i, ac.pos[2] - flipped * 2)]);
                                
                                density[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = source;
                                u[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = flipped * angleF_l;
                                v[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = 0.000001f;
                                w[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = flipped * angleS_l;
                                 //printf("- %d %d %d: %8f %8f %8f\n", ac.pos[0], i, ac.pos[2] + flipped * 1, u[IX(ac.pos[0], i, ac.pos[2] + flipped * 1)], v[IX(ac.pos[0], i, ac.pos[2] + flipped * 1)], w[IX(ac.pos[0], i, ac.pos[2] + flipped * 1)]);

                            }
                        }
                        else {
                            for (int i = 5; i <= 9; i++)
                            {
                                density[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = source;
                                u[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = flipped * angleS_r;
                                v[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = 0.000001f;
                                w[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = flipped * angleF_r;

                                density[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = source;
                                u[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = -flipped * angleS_l;
                                v[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = 0.000001f;
                                w[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = flipped * angleF_l;
                            }
                        }

                        angle_r = angleGuard_tr * 3.141592f / 180.f;
                        angle_l = angleGuard_tl * 3.141592f / 180.f;
                        angleF_r = (angleGuard_tr < 0) ? 0 : Gscale * outV * cos(angle_r);
                        angleS_r = (angleGuard_tr < 0) ? 0 : Gscale * outV * sin(angle_r);
                        angleF_l = (angleGuard_tl < 0) ? 0 : Gscale * outV * cos(angle_l);
                        angleS_l = (angleGuard_tl < 0) ? 0 : Gscale * outV * sin(angle_l);
                        
                        if (isXdir) {
                            for (int i = 10; i <= 14; i++)
                            {
                                density[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = source;
                                u[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = flipped * angleF_r;
                                v[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = 0.000001f;
                                w[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] - flipped * 2)] = -flipped * angleS_r;

                                density[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = source;
                                u[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = flipped * angleF_l;
                                v[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = 0.000001f;
                                w[IX(ac.pos[0], ac.pos[1] + i, ac.pos[2] + flipped * 1)] = flipped * angleS_l;
                            }
                        }
                        else {
                            for (int i = 10; i <= 14; i++)
                            {
                                density[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = source;
                                u[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = flipped * angleS_r;
                                v[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = 0.000001f;
                                w[IX(ac.pos[0] + flipped * 2, ac.pos[1] + i, ac.pos[2])] = flipped * angleF_r;

                                density[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = source;
                                u[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = -flipped * angleS_l;
                                v[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = 0.000001f;
                                w[IX(ac.pos[0] - flipped * 1, ac.pos[1] + i, ac.pos[2])] = flipped * angleF_l;
                            }
                        }
                    }   //end guard
                    
                    
                }
                //diagonal direction
                else if (ac.direction < 8){
                    int xDiff_r = (ac.direction == 5 || ac.direction == 7) ? 1 : -1;
                    int zDiff_r = (ac.direction == 3 || ac.direction == 5) ? 1 : -1;
                    
                    angleCircle = 5;
                    {   //circle
                        float angleA = angleCircle * 3.141592f / 180.f;
                        Cscale += 0.5f;

                        density[IX(ac.pos[0], acHeight, ac.pos[2])] = source;
                        u[IX(ac.pos[0], acHeight, ac.pos[2])] = -2 * Cscale * outV * zDiff_r / sqrt(2) * cos(angleA);
                        v[IX(ac.pos[0], acHeight, ac.pos[2])] = -2 * Cscale * outV * sin(angleA);
                        w[IX(ac.pos[0], acHeight, ac.pos[2])] = 2 * Cscale * outV * xDiff_r / sqrt(2) * cos(angleA);
                        
                        density[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = source;
                        u[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = -2 * Cscale * outV * zDiff_r / sqrt(2) * cos(angleA);
                        v[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = 2 * Cscale * outV * sin(angleA);
                        w[IX(ac.pos[0], acHeight + 3, ac.pos[2])] = 2 * Cscale * outV * xDiff_r / sqrt(2) * cos(angleA);

                        density[IX(ac.pos[0], acHeight + 1, ac.pos[2])] = source;
                        u[IX(ac.pos[0], acHeight + 1, ac.pos[2])] = -2 * Cscale * outV * zDiff_r / sqrt(2) * cos(angleA);
                        v[IX(ac.pos[0], acHeight + 1, ac.pos[2])] = -2 * Cscale * outV * sin(angleA);
                        w[IX(ac.pos[0], acHeight + 1, ac.pos[2])] = 2 * Cscale * outV * xDiff_r / sqrt(2) * cos(angleA);

                        density[IX(ac.pos[0], acHeight + 2, ac.pos[2])] = source;
                        u[IX(ac.pos[0], acHeight + 2, ac.pos[2])] = -2 * Cscale * outV * zDiff_r / sqrt(2) * cos(angleA);
                        v[IX(ac.pos[0], acHeight + 2, ac.pos[2])] = 2 * Cscale * outV * sin(angleA);
                        w[IX(ac.pos[0], acHeight + 2, ac.pos[2])] = 2 * Cscale * outV * xDiff_r / sqrt(2) * cos(angleA);


                        float outVX = Cscale * outV / sqrt(2) * (sqrt(3) / 2 * xDiff_r * sin(angleA) - zDiff_r * cos(angleA));
                        float outVZ = Cscale * outV / sqrt(2) * (sqrt(3) / 2 * zDiff_r * sin(angleA) + xDiff_r * cos(angleA));
                        float outVY = Cscale * outV / 2 * sin(angleA);

                        density[IX(ac.pos[0] + xDiff_r, acHeight + 1, ac.pos[2] + zDiff_r)] = source;
                        u[IX(ac.pos[0] + xDiff_r, acHeight + 1, ac.pos[2] + zDiff_r)] = outVX;
                        v[IX(ac.pos[0] + xDiff_r, acHeight + 1, ac.pos[2] + zDiff_r)] = -outVY;
                        w[IX(ac.pos[0] + xDiff_r, acHeight + 1, ac.pos[2] + zDiff_r)] = outVZ;

                        density[IX(ac.pos[0] + xDiff_r, acHeight + 2, ac.pos[2] + zDiff_r)] = source;
                        u[IX(ac.pos[0] + xDiff_r, acHeight + 2, ac.pos[2] + zDiff_r)] = outVX;
                        v[IX(ac.pos[0] + xDiff_r, acHeight + 2, ac.pos[2] + zDiff_r)] = outVY;
                        w[IX(ac.pos[0] + xDiff_r, acHeight + 2, ac.pos[2] + zDiff_r)] = outVZ;

                        int flipped = xDiff_r == zDiff_r ? -1 : 1;
                        density[IX(ac.pos[0] - xDiff_r, acHeight + 1, ac.pos[2] - zDiff_r)] = source;
                        u[IX(ac.pos[0] - xDiff_r, acHeight + 1, ac.pos[2] - zDiff_r)] = outVZ * flipped;
                        v[IX(ac.pos[0] - xDiff_r, acHeight + 1, ac.pos[2] - zDiff_r)] = -outVY;
                        w[IX(ac.pos[0] - xDiff_r, acHeight + 1, ac.pos[2] - zDiff_r)] = outVX * flipped;

                        density[IX(ac.pos[0] - xDiff_r, acHeight + 2, ac.pos[2] - zDiff_r)] = source;
                        u[IX(ac.pos[0] - xDiff_r, acHeight + 2, ac.pos[2] - zDiff_r)] = outVZ * flipped;
                        v[IX(ac.pos[0] - xDiff_r, acHeight + 2, ac.pos[2] - zDiff_r)] = outVY;
                        w[IX(ac.pos[0] - xDiff_r, acHeight + 2, ac.pos[2] - zDiff_r)] = outVX * flipped;

                    }   //end circle

                    {   //guard
                        float angle_r = angleGuard_br * 3.141592f / 180.f;
                        float angle_l = angleGuard_bl * 3.141592f / 180.f;
                        float angleX_r = angleGuard_br < 0 ? 0 : (Gscale * outV / sqrt(2) * (xDiff_r * sin(angle_r) - zDiff_r * cos(angle_r)));
                        float angleZ_r = angleGuard_br < 0 ? 0 : (Gscale * outV / sqrt(2) * (zDiff_r * sin(angle_r) + xDiff_r * cos(angle_r)));
                        float angleX_l = angleGuard_bl < 0 ? 0 : (Gscale * outV / sqrt(2) * (-xDiff_r * sin(angle_l) - zDiff_r * cos(angle_l)));
                        float angleZ_l = angleGuard_bl < 0 ? 0 : (Gscale * outV / sqrt(2) * (-zDiff_r * sin(angle_l) + xDiff_r * cos(angle_l)));

                        int xx = ac.pos[0] + xDiff_r;
                        int zz = ac.pos[2] + zDiff_r;

                        for (int i = 5; i <= 9; i++)
                        {
                            density[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = source;
                            u[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = angleX_r;
                            v[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = 0.000001f;
                            w[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = angleZ_r;

                            density[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = source;
                            u[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = angleX_l;
                            v[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = 0.000001f;
                            w[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = angleZ_l;
                        }

                        angle_r = angleGuard_tr * 3.141592f / 180.f;
                        angle_l = angleGuard_tl * 3.141592f / 180.f;
                        angleX_r = angleGuard_tr < 0 ? 0 : (Gscale * outV / sqrt(2) * (xDiff_r * sin(angle_r) - zDiff_r * cos(angle_r)));
                        angleZ_r = angleGuard_tr < 0 ? 0 : (Gscale * outV / sqrt(2) * (zDiff_r * sin(angle_r) + xDiff_r * cos(angle_r)));
                        angleX_l = angleGuard_tl < 0 ? 0 : (Gscale * outV / sqrt(2) * (-xDiff_r * sin(angle_l) - zDiff_r * cos(angle_l)));
                        angleZ_l = angleGuard_tl < 0 ? 0 : (Gscale * outV / sqrt(2) * (-zDiff_r * sin(angle_l) + xDiff_r * cos(angle_l)));
                         
                        for (int i = 10; i <= 14; i++)
                        {
                            density[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = source;
                            u[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = angleX_r;
                            v[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = 0.000001f;
                            w[IX(ac.pos[0] + xDiff_r, ac.pos[1] + i, ac.pos[2] + zDiff_r)] = angleZ_r;

                            density[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = source;
                            u[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = angleX_l;
                            v[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = 0.000001f;
                            w[IX(ac.pos[0] - xDiff_r, ac.pos[1] + i, ac.pos[2] - zDiff_r)] = angleZ_l;
                        }

                    }   //end guard

                }
            }
        }

        addDensity = false;
        return;

    }

}

void Solver::updateObstacleGPU()
{
    //h_obstacle = obstacle;
    int size = n3b * sizeof(bool);

    checkCudaErrors(cudaMemcpyAsync(d_obstacle, h_obstacle, size, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaMemcpyAsync(d_obstacleU, h_obstacleU, size, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaMemcpyAsync(d_obstacleV, h_obstacleV, size, cudaMemcpyHostToDevice));
    //checkCudaErrors(cudaMemcpyAsync(d_obstacleW, h_obstacleW, size, cudaMemcpyHostToDevice));
}

void Solver::velocityStepGPU()
{
    // Allocate CUDA events that we'll use for timing
    checkCudaErrors(cudaMemcpyAsync(d_preVelocityU, h_preVelocityU, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyAsync(d_preVelocityV, h_preVelocityV, memSize, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpyAsync(d_preVelocityW, h_preVelocityW, memSize, cudaMemcpyHostToDevice));

    // add source.
    add_sourceGPU(d_velocityU, d_preVelocityU, true);
    add_sourceGPU(d_velocityV, d_preVelocityV, true);
    add_sourceGPU(d_velocityW, d_preVelocityW, true);
    
    // Gravity and Bouyancy;
    add_buoyancyGPU(d_velocityV, d_density);

    std::swap(d_velocityU, d_preVelocityU);
    std::swap(d_velocityV, d_preVelocityV);
    std::swap(d_velocityW, d_preVelocityW);

    // diffuse 
    diffuseGPU(d_velocityU, d_preVelocityU, viscosityCoef, 1);
    diffuseGPU(d_velocityV, d_preVelocityV, viscosityCoef, 2);
    diffuseGPU(d_velocityW, d_preVelocityW, viscosityCoef, 3);
    projectGPU(d_velocityU, d_velocityV, d_velocityW, d_preVelocityU, d_preVelocityV);

    std::swap(d_velocityU, d_preVelocityU);
    std::swap(d_velocityV, d_preVelocityV);
    std::swap(d_velocityW, d_preVelocityW);
    

    advectGPU(d_velocityU, d_preVelocityU, d_preVelocityU, d_preVelocityV, d_preVelocityW, 1);
    advectGPU(d_velocityV, d_preVelocityV, d_preVelocityU, d_preVelocityV, d_preVelocityW, 2);
    advectGPU(d_velocityW, d_preVelocityW, d_preVelocityU, d_preVelocityV, d_preVelocityW, 3);

    projectGPU(d_velocityU, d_velocityV, d_velocityW, d_preVelocityU, d_preVelocityV);
    checkCudaErrors(cudaMemcpyAsync(h_velocityU, d_velocityU, memSize, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpyAsync(h_velocityV, d_velocityV, memSize, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpyAsync(h_velocityW, d_velocityW, memSize, cudaMemcpyDeviceToHost));
}

void Solver::densityStepGPU()
{
    // Allocate CUDA events that we'll use for timing

    checkCudaErrors(cudaMemcpyAsync(d_preDensity, h_preDensity, memSize, cudaMemcpyHostToDevice));

    add_sourceGPU(d_density, d_preDensity, true);

    //std::swap(d_density, d_preDensity);
    checkCudaErrors(cudaMemcpyAsync(d_preDensity, d_density, memSize, cudaMemcpyDeviceToDevice));
    
    diffuseGPU(d_density, d_preDensity, diffuseCoef,  0, true);
    std::swap(d_density, d_preDensity);

    advectGPU(d_density, d_preDensity, d_velocityU, d_velocityV, d_velocityW, 0, true);

    checkCudaErrors(cudaMemcpyAsync(h_density, d_density, memSize, cudaMemcpyDeviceToHost));
}

float Solver::getVelocityGPUU(int x, int y, int z)
{
    const std::lock_guard<std::mutex> lock(mut);
    return h_velocityU[IX(x, y, z)];
}

float Solver::getVelocityGPUV(int x, int y, int z)
{
    const std::lock_guard<std::mutex> lock(mut);
    return h_velocityV[IX(x, y, z)];
}

float Solver::getVelocityGPUW(int x, int y, int z)
{
    const std::lock_guard<std::mutex> lock(mut);
    return h_velocityW[IX(x, y, z)];
}

const float* Solver::getVelocityGPUN()
{
    //const std::lock_guard<std::mutex> lock(mut);
    int c = 0;

#pragma omp parallel for
    FOR_XYZ
        h_velocityN[c] = getVelocityGPUU(x, y, z); c++;
        h_velocityN[c] = getVelocityGPUV(x, y, z); c++;
        h_velocityN[c] = getVelocityGPUW(x, y, z); c++;
    FOR_END

    return h_velocityN;
}

void Solver::getVelocityGPUN(float* v)
{
    const std::lock_guard<std::mutex> lock(mut);
    int c = 0;

#pragma omp parallel for
    FOR_XYZ
        v[c] = h_velocityU[IX(x, y, z)]; c++;
        v[c] = h_velocityV[IX(x, y, z)]; c++;
        v[c] = h_velocityW[IX(x, y, z)]; c++;
    FOR_END
    //auto vel = getVelocityGPUN();

    //memcpy(v, vel, n3 * sizeof(float) * 3);
}

void Solver::getVelocityGPUN2(float* v)
{
    const std::lock_guard<std::mutex> lock(mut);
    int c = 0;

#pragma omp parallel for
    FOR_XYZ0
        c = IX(x, y, z);
        v[3 * c] = h_velocityU[c];
        v[3 * c + 1] = h_velocityV[c];
        v[3 * c + 2] = h_velocityW[c];
    FOR_END
}


float Solver::getDensityAvgGPU()
{
    const std::lock_guard<std::mutex> lock(mut);

    float v = 0.0f;
#pragma omp parallel for
    FOR_XYZ
        v += h_density[IX(x, y, z)];
    FOR_END

    return v / (nx*ny*nz);
}

float Solver::getDensityGPU(int x, int y, int z)
{
    const std::lock_guard<std::mutex> lock(mut);
    return h_density[IX(x, y, z)];
}

const float* Solver::getDensityGPUN()
{
    //const std::lock_guard<std::mutex> lock(mut);
    int c = 0;

#pragma omp parallel for
    FOR_XYZ
        h_densityN[c] = getDensityGPU(x, y, z);
        c++;
    FOR_END

        return h_densityN;
}

void Solver::getDensityGPUN(float* d)
{
    const std::lock_guard<std::mutex> lock(mut);
    //auto den = getDensityGPUN();
    int c = 0;

#pragma omp parallel for
    FOR_XYZ
        d[c] = h_density[IX(x, y, z)];
        c++;
    FOR_END

    //memcpy(d, den, n3 * sizeof(float));
}

void Solver::ResetHalfBoundary()
{
    AirConditioner& ac = airconditioners[0];
    if (ac.acType == ACTYPE::Tower){
        for (int dx = -1; dx <= 1; dx++) {
            for (int dz = -1; dz <= 1; dz++) {
                setObstacle(ac.pos[0] + dx, ac.pos[1] + 16, ac.pos[2] + dz, false);
                setObstacle(ac.pos[0] + dx, ac.pos[1] + 17, ac.pos[2] + dz, false);
            }
        }
    }
    else if (ac.acType == ACTYPE::Stand){

        int widthSide = 2;
        int heightSide = 1;

        bool isACPosX = (ac.direction % 4 == 0);
        int reverse = ac.direction == 0 || ac.direction == 6 ? 1 : -1;

        int acW = (isACPosX ? 1 : widthSide);
        int acT = (isACPosX ? widthSide : 1);
        int acH = heightSide;

        for (int i = -acW; i <= acW; i++)
        {
            for (int j = -acH; j <= acH; j++)
            {
                for (int k = -acT; k <= acT; k++)
                {
                    int win_x = ac.pos[0] + i;
                    int win_y = ac.pos[1] + j;
                    int win_z = ac.pos[2] + k;

                    setObstacle(win_x, win_y, win_z, false);
                }
            }
        }
    }
    else if (ac.acType == ACTYPE::Wall) {

        int widthSide = 3;
        int heightSide = 1;

        bool isACPosX = (ac.direction % 4 == 0);
        int reverse = ac.direction == 0 || ac.direction == 6 ? 1 : -1;

        int acW = (isACPosX ? 0 : widthSide);
        int acT = (isACPosX ? widthSide : 0);
        int acH = heightSide;

        for (int i = -acW; i <= acW; i++)
        {
            for (int j = -acH; j <= acH; j++)
            {
                for (int k = -acT; k <= acT; k++)
                {
                    int win_x = ac.pos[0] + i;
                    int win_y = ac.pos[1] + j;
                    int win_z = ac.pos[2] + k;

                    setObstacle(win_x, win_y, win_z, false);
                }
            }
        }

    }
    else if (ac.acType == ACTYPE::Ceiling) {

        int offset = 4;
        int halfcell = offset / 2;

        auto acHeight = ac.pos[1] - 1;

        for (int j = -1; j < 2; j++)
        {
            for (int i = -halfcell + abs(j); i <= halfcell - abs(j); i++)
            {
                setObstacle(ac.pos[0] + (offset + j), acHeight, ac.pos[2] + i, false);
                setObstacle(ac.pos[0] - (offset + j), acHeight, ac.pos[2] + i, false);
                setObstacle(ac.pos[0] + i, acHeight, ac.pos[2] + (offset + j), false);
                setObstacle(ac.pos[0] + i, acHeight, ac.pos[2] - (offset + j), false);
            }
        }
    }
}

void Solver::SetHalfBoundary()
{
    AirConditioner& ac = airconditioners[0];
    if (ac.acType != ACTYPE::Tower) return;
    
    //printf("(%d,%d,%d), %d\n", ac.pos[0], ac.pos[1], ac.pos[2], ac.direction);
    if (ac.direction < 8 && ac.direction % 2 == 0) {
        //bool isXdir = ac.direction % 4 == 0;
        //int flipped = isXdir ? (ac.direction == 0 ? 1 : -1) : (ac.direction == 6 ? 1 : -1);

        //if (isXdir) {
        //    int min = ac.pos[2] - (flipped ? 1 : 2);
        //    int max = ac.pos[2] + (flipped ? 2 : 1);
        //    for (int y = 6; y <= 19; y++) {
        //        for (int z = min; z <= max; z++) {
        //            if (y <= 15) {
        //                h_obstacleU[IX(ac.pos[0], y, z)] = true;
        //                h_obstacleU[IX(ac.pos[0] + 1, y, z)] = true;
        //            }
        //        }
        //        h_obstacleW[IX(ac.pos[0], y, min)] = true;
        //        h_obstacleW[IX(ac.pos[0], y, max + 1)] = true;
        //    }
        //    for (int z = min; z <= max; z++) {
        //        h_obstacleV[IX(ac.pos[0], 6, z)] = true;
        //        h_obstacleV[IX(ac.pos[0], 20, z)] = true;
        //    }

        //}
        //else {
        //    int min = ac.pos[0] - (flipped ? 1 : 2);
        //    int max = ac.pos[0] + (flipped ? 2 : 1);
        //    for (int y = 6; y <= 19; y++) {
        //        for (int x = min; x <= max; x++) {
        //            if (y <= 15) {
        //                h_obstacleW[IX(x, y, ac.pos[2])] = true;
        //                h_obstacleW[IX(x, y, ac.pos[2] + 1)] = true;
        //            }
        //        }
        //        h_obstacleU[IX(min, y, ac.pos[2])] = true;
        //        h_obstacleU[IX(max + 1, y, ac.pos[2])] = true;
        //    }
        //    for (int x = min; x <= max; x++) {
        //        h_obstacleV[IX(x, 6, ac.pos[2])] = true;
        //        h_obstacleV[IX(x, 20, ac.pos[2])] = true;
        //    }

        //}

        int dx = (ac.direction == 0 || ac.direction == 6) ? +1 : -1;
        int dz = (ac.direction == 0 || ac.direction == 2) ? -1 : +1;

        setObstacle(ac.pos[0], ac.pos[1] + 16, ac.pos[2] + dz);
        setObstacle(ac.pos[0], ac.pos[1] + 17, ac.pos[2] + dz);
        setObstacle(ac.pos[0], ac.pos[1] + 16, ac.pos[2]);
        setObstacle(ac.pos[0], ac.pos[1] + 17, ac.pos[2]);

        setObstacle(ac.pos[0] + dx, ac.pos[1] + 16, ac.pos[2] + dz);
        setObstacle(ac.pos[0] + dx, ac.pos[1] + 17, ac.pos[2] + dz);
        setObstacle(ac.pos[0] + dx, ac.pos[1] + 16, ac.pos[2]);
        setObstacle(ac.pos[0] + dx, ac.pos[1] + 17, ac.pos[2]);


    }
    else if (ac.direction < 8 && ac.direction % 2 == 1) {
        //diagonal direction


    }
    //if (ac.direction % 2 == 0 && ac.direction < 8) {
    //    bool isXdir = ac.direction % 4 == 0;
    //    int flipped = isXdir ? (ac.direction == 0 ? 1 : -1) : (ac.direction == 6 ? 1 : -1);
    //    for (int i = 6; i <= 15; i++) {



    //        //ac.pos[0], acHeight, ac.pos[2]
    //    }
    //}
}
