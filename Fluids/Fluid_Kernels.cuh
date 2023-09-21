#pragma once

#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>


__global__ void advectionKernel(float* __restrict__ dst, float* __restrict__ src,
                          float* __restrict__ xu, float* __restrict__ xv, float* __restrict__ xw, 
                          int3 dim, float time_step);

__global__ void addSourceKernel(float* __restrict__ dst, float* __restrict__ src,  int3 dim, float time_step, bool copy);

__global__ void addBuoyancyKernel(float *velosityU, float* __restrict__ temperature, const int3 dim, float time_step, float ambiTemperature, float coeffBouyance, float gravity, bool* __restrict__ occ);

__global__ void setBoundaryKernel(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag);

__global__ void setBoundaryKernelDir(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag, float target);

__global__ void setBoundaryKernel2(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag);

__global__ void divergenceKernel(float* __restrict__ div, float* __restrict__ xu, float* __restrict__ xv, float* __restrict__ xw, const int3 dim, bool* __restrict__ occ);

__global__ void lin_solveKernel(float* __restrict__ dst, float* __restrict__ neighbor_src, float* __restrict__ src, const int3 dim, float a, bool sum, bool* __restrict__ occ);

__global__ void applyPressureKernel(float* __restrict__ u, float* __restrict__ v, float* __restrict__ w, float* __restrict__ p, const int3 dim, bool* __restrict__ occ);
