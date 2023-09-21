#include "fluid_kernels.cuh"


inline __host__ __device__ float3 operator-(float3 a, float3 b)
{
	return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline __host__ __device__ float3 operator-(float3 a, float b)
{
	return make_float3(a.x - b, a.y - b, a.z - b);
}

inline __host__ __device__ float3 operator*(float3 a, float b)
{
	return make_float3(a.x * b, a.y * b, a.z * b);
}

inline __host__ __device__ float3 operator*(float3 a, float3 b)
{
	return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline __host__ __device__ int3 make_int3(float3 a)
{
	return make_int3(int(a.x), int(a.y), int(a.z));
}

inline __host__ __device__ float3 floorr(const float3 v)
{
	return make_float3(floorf(v.x), floorf(v.y), floorf(v.z));
}

inline __device__ int get_voxel(int x, int y, int z, int3 d)
{
	__assume(x < 256);
	__assume(y < 256);
	__assume(z < 256);
	__assume(d.x < 256);
	__assume(d.y < 256);
	__assume(d.z < 256);

	return z + (x * d.y + y) * d.z;
}

inline __device__ float get_cell(int3 c, int3 d, float* __restrict__ src)
{
	__assume(c.x < 256);
	__assume(c.y < 256);
	__assume(c.z < 256);

	// 0 <= x < size
	if (c.x < 0) c.x = 0;
	else if (c.x > (d.x - 1)) c.x = d.x - 1;
	if (c.y < 0) c.y = 0;
	else if (c.y > (d.y - 1)) c.y = d.y - 1;
	if (c.z < 0) c.z = 0;
	else if (c.z > (d.z - 1)) c.z = d.z - 1;
	
	return src[get_voxel(c.x, c.y, c.z, d)];
}

inline __device__ float get_cellF(float3 p, int3 d, float* __restrict__ src)
{
	// bilinear interpolation
	float3 l = floorr(p);
	int3 rp = make_int3(l);
	float3 dif = p - l;
//	float sum = 0.f;

	float pretty[8];
	pretty[0] = abs( ((1.0f - 0.0f) - dif.x) * ((1.0f - 0.0f) - dif.y) * ((1.0f - 0.0f) - dif.z) ) * get_cell(make_int3(rp.x + 0, rp.y + 0, rp.z + 0), d, src);
	pretty[1] = abs( ((1.0f - 0.0f) - dif.x) * ((1.0f - 0.0f) - dif.y) * ((1.0f - 1.0f) - dif.z) ) * get_cell(make_int3(rp.x + 0, rp.y + 0, rp.z + 1), d, src);
	pretty[2] = abs( ((1.0f - 0.0f) - dif.x) * ((1.0f - 1.0f) - dif.y) * ((1.0f - 0.0f) - dif.z) ) * get_cell(make_int3(rp.x + 0, rp.y + 1, rp.z + 0), d, src);
	pretty[3] = abs( ((1.0f - 0.0f) - dif.x) * ((1.0f - 1.0f) - dif.y) * ((1.0f - 1.0f) - dif.z) ) * get_cell(make_int3(rp.x + 0, rp.y + 1, rp.z + 1), d, src);
	pretty[4] = abs( ((1.0f - 1.0f) - dif.x) * ((1.0f - 0.0f) - dif.y) * ((1.0f - 0.0f) - dif.z) ) * get_cell(make_int3(rp.x + 1, rp.y + 0, rp.z + 0), d, src);
	pretty[5] = abs( ((1.0f - 1.0f) - dif.x) * ((1.0f - 0.0f) - dif.y) * ((1.0f - 1.0f) - dif.z) ) * get_cell(make_int3(rp.x + 1, rp.y + 0, rp.z + 1), d, src);
	pretty[6] = abs( ((1.0f - 1.0f) - dif.x) * ((1.0f - 1.0f) - dif.y) * ((1.0f - 0.0f) - dif.z) ) * get_cell(make_int3(rp.x + 1, rp.y + 1, rp.z + 0), d, src);
	pretty[7] = abs( ((1.0f - 1.0f) - dif.x) * ((1.0f - 1.0f) - dif.y) * ((1.0f - 1.0f) - dif.z) ) * get_cell(make_int3(rp.x + 1, rp.y + 1, rp.z + 1), d, src);


	return pretty[0] + pretty[1] + pretty[2] + pretty[3] + pretty[4] + pretty[5] + pretty[6] + pretty[7];
}

inline __device__ bool isInside(int x, int y, int z, const int3 dim)
{
	// the grid size is N+2 including boundary
	//  1 <= x <= N; 
	return !(x < 1 || x >(dim.x - 2) ||
		y < 1 || y >(dim.y - 2) ||
		z < 1 || z >(dim.z - 2));
}

__global__ void advectionKernel(float* __restrict__ dst, float* __restrict__ src, float* __restrict__ xu, float* __restrict__ xv, float* __restrict__ xw, const int3 dim, float time_step)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;

	if (__builtin_expect(!isInside(x, y, z, dim), false))
		return;

	int vox = get_voxel(x, y, z, dim);
	float avg = (dim.x + dim.y + dim.z - 6) / 3.0f;

	float3 np = make_float3(float(x), float(y), float(z))
		- make_float3(xu[vox] * avg, xv[vox] * avg, xw[vox] * avg) * time_step;
		
	if (np.x < 0) np.x = 1;
	else if (np.x > (dim.x - 2)) np.x = dim.x - 2;

	if (np.y < 0) np.y = 1;
	else if (np.y > (dim.y - 2)) np.y = dim.y - 2;

	if (np.z < 0) np.z = 1;
	else if (np.z > (dim.z - 2)) np.z = dim.z - 2;

	dst[vox] = get_cellF(np, dim, src);
}

__global__ void addSourceKernel(float* __restrict__ dst, float* __restrict__ src, const int3 dim, float time_step, bool copy)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;

	//if (__builtin_expect(!isInside(x, y, z, dim), false))
	//    return;

	if (x >= dim.x ||
		y >= dim.y ||
		z >= dim.z)
		return;

	int vox = get_voxel(x, y, z, dim);
	if (copy) {
		if (src[vox] != 0.0f)
		{
			dst[vox] = src[vox];
		}
	}
	else
	{
		dst[vox] += src[vox];
	}
		//dst[vox] = /*dst[vox] +  time_step **/ src[vox];
}

//TD
__global__ void addBuoyancyKernel(float* __restrict__ velocityV, float* __restrict__ temperature, const int3 dim, float time_step, float ambiTemperature, float coeffBouyance, float gravity, bool* __restrict__ occ)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;

	if (__builtin_expect(!isInside(x, y, z, dim), false))
		return;

	int vox = get_voxel(x, y, z, dim);
	int voxU = get_voxel(x, y + 1, z, dim);
	int voxD = get_voxel(x, y - 1, z, dim);
	int voxR = get_voxel(x + 1, y, z, dim);
	int voxL = get_voxel(x - 1, y, z, dim);
	int voxF = get_voxel(x, y, z + 1, dim);
	int voxB = get_voxel(x, y, z - 1, dim);

	int neighborCount = 1;
	float avgTemp = temperature[vox];

	if (__builtin_expect(!occ[voxL], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxL];
	}
	if (__builtin_expect(!occ[voxR], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxR];
	}
	if (__builtin_expect(!occ[voxD], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxD];
	}
	if (__builtin_expect(!occ[voxU], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxU];
	}
	if (__builtin_expect(!occ[voxB], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxB];
	}
	if (__builtin_expect(!occ[voxF], true)) {
		neighborCount += 1;
		avgTemp += temperature[voxF];
	}

	avgTemp /= neighborCount;

	float diff = temperature[vox] - avgTemp;


	//float avgtemp = temperature[vox] - (temperature[voxU] + temperature[voxD]+
	//									temperature[voxL] + temperature[voxR]+
	//									temperature[voxB] + temperature[voxF]) / 6.0f;

//	float diff = (temperature[voxU] + avgtemp * 2.0f
//		+ temperature[voxR] + temperature[voxL] + temperature[voxF] + temperature[voxK]) / 7.f;

	//0.01:gravity term
	//avgtemp:bouyancy term
	velocityV[vox] += time_step * coeffBouyance * gravity * 0.01f * diff;
}

//TD
__global__ void setBoundaryKernel(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag)
{
	// neumann boundary condition is f' = a ;
	// we determin a = 0.

	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;


	if (x >= dim.x ||
		y >= dim.y ||
		z >= dim.z)
		return;

	int vox = get_voxel(x, y, z, dim);

	////////////////////////////////////////////////////////////////////////////////
	// out boundary 

	if (__builtin_expect(isInside(x, y, z, dim), true))
	{
		float dst_vox = dst[vox];
		// for X
		if (flag == 1)
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = -dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = -dst_vox;
		}
		else
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = dst_vox;
		}


		// for Y
		if (flag == 2)
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = -dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = -dst_vox;
		}
		else
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = dst_vox;
		}

		// for Z
		if (flag == 3)
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = -dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = -dst_vox;
		}
		else
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = dst_vox;
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// inside(object) boundary 

	if (occ != nullptr & occ[vox])
	{
		dst[vox] = 0.f;

		if (flag == 1)
		{ // x direction
			if(x > 0)
				dst[get_voxel(x - 1, y, z, dim)] = 0.f;
			if(x < dim.x - 1)
				dst[get_voxel(x + 1, y, z, dim)] = 0.f;
		}
		else if (flag == 2)
		{ // y direction
			if (y > 0)
				dst[get_voxel(x, y - 1, z, dim)] = 0.f;
			if (y < dim.y - 1)
				dst[get_voxel(x, y + 1, z, dim)] = 0.f;
		}
		else if (flag == 3)
		{ // z direction
			if (z > 0)
				dst[get_voxel(x, y, z - 1, dim)] = 0.f;
			if (z < dim.z - 1)
				dst[get_voxel(x, y, z + 1, dim)] = 0.f;
		}
	}


	////////////////////////////////////////////////////////////////////////////////
	// corner

	__syncthreads();


	if (__builtin_expect(z == 0, false))
	{
		// front side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
	}
	else if (__builtin_expect(z == (dim.z - 1), false))
	{
		// back side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
	}

}

//TD
__global__ void setBoundaryKernelDir(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag, float target)
{
	// neumann boundary condition is f' = a ;
	// we determin a = 0.


	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;


	if (x >= dim.x ||
		y >= dim.y ||
		z >= dim.z)
		return;

	int vox = get_voxel(x, y, z, dim);

	////////////////////////////////////////////////////////////////////////////////
	// out boundary 



	if (__builtin_expect(isInside(x, y, z, dim), true))
	{
		float dst_vox = dst[vox];

		// for X
		if (flag == 1)
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = dst_vox;
		}
		else
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = dst_vox;
		}


		// for Y
		if (flag == 2)
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = dst_vox;
		}
		else
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = dst_vox;
		}

		// for Z
		if (flag == 3)
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = dst_vox;
		}
		else
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = dst_vox;
		}

		////////////////////////////////////////////////////////////////////////////////
		// inside(object) boundary 

		if (flag == 0)
		{
			if (occ != nullptr & occ[vox])
			{
				float  l = dst[get_voxel(x - 1, y, z, dim)];
				float  r = dst[get_voxel(x + 1, y, z, dim)];

				float  u = dst[get_voxel(x, y + 1, z, dim)];
				float  d = dst[get_voxel(x, y - 1, z, dim)];

				float  f = dst[get_voxel(x, y, z - 1, dim)];
				float  b = dst[get_voxel(x, y, z + 1, dim)];

				dst[vox] = (l + r + u + d + f + b) / 6.f;
			}
		}
	} // end inside


	////////////////////////////////////////////////////////////////////////////////
	// corner

//	__syncthreads();

	if (__builtin_expect(z == 0, false))
	{
		// front side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
	}
	else if (__builtin_expect(z == (dim.z - 1), false))
	{
		// back side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
	}

}

//TD
__global__ void setBoundaryKernel2(float* __restrict__ dst, bool* __restrict__ occ, const int3 dim, int flag, float* __restrict__ target)
{
	// neumann boundary condition is f' = a ;
	// we determin a = 0.

	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;


	if (x >= dim.x ||
		y >= dim.y ||
		z >= dim.z)
		return;

	int vox = get_voxel(x, y, z, dim);

	////////////////////////////////////////////////////////////////////////////////
	// out boundary 

	if (__builtin_expect(isInside(x, y, z, dim), true))
	{
		float dst_vox = dst[vox];
		// for X
		if (flag == 1)
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = -dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = -dst_vox;
		}
		else
		{
			if (x == 1)
				dst[get_voxel(x - 1, y, z, dim)] = dst_vox;
			else if (x == (dim.x - 2))
				dst[get_voxel(x + 1, y, z, dim)] = dst_vox;
		}


		// for Y
		if (flag == 2)
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = -dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = -dst_vox;
		}
		else
		{
			if (y == 1)
				dst[get_voxel(x, y - 1, z, dim)] = dst_vox;
			else if (y == (dim.y - 2))
				dst[get_voxel(x, y + 1, z, dim)] = dst_vox;
		}

		// for Z
		if (flag == 3)
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = -dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = -dst_vox;
		}
		else
		{
			if (z == 1)
				dst[get_voxel(x, y, z - 1, dim)] = dst_vox;
			else if (z == (dim.z - 2))
				dst[get_voxel(x, y, z + 1, dim)] = dst_vox;
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// inside(object) boundary 

	if (occ != nullptr & occ[vox])
	{
		dst[vox] = 0.f;

		if (flag == 1)
		{ // x direction
			if (x > 0)
				dst[get_voxel(x - 1, y, z, dim)] = 0.f;
			if (x < dim.x - 1)
				dst[get_voxel(x + 1, y, z, dim)] = 0.f;
		}
		else if (flag == 2)
		{ // y direction
			if (y > 0)
				dst[get_voxel(x, y - 1, z, dim)] = 0.f;
			if (y < dim.y - 1)
				dst[get_voxel(x, y + 1, z, dim)] = 0.f;
		}
		else if (flag == 3)
		{ // z direction
			if (z > 0)
				dst[get_voxel(x, y, z - 1, dim)] = 0.f;
			if (z < dim.z - 1)
				dst[get_voxel(x, y, z + 1, dim)] = 0.f;
		}
	}


	////////////////////////////////////////////////////////////////////////////////
	// corner

	__syncthreads();


	if (__builtin_expect(z == 0, false))
	{
		// front side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z + 1, dim)]) / 3.f;
		}
	}
	else if (__builtin_expect(z == (dim.z - 1), false))
	{
		// back side
		if (y == 0)
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y + 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
		else if (y == (dim.y - 1))
		{
			if (x == 0)
				dst[vox] = (dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x + 1, y, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
			else if (x == (dim.x - 1))
				dst[vox] = (dst[get_voxel(x - 1, y, z, dim)] + dst[get_voxel(x, y - 1, z, dim)] + dst[get_voxel(x, y, z - 1, dim)]) / 3.f;
		}
	}

}

//TD
__global__ void divergenceKernel(float* __restrict__ div, float* __restrict__ xu, float* __restrict__ xv, float* __restrict__ xw, const int3 dim, bool* __restrict__ occ)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;

	if (__builtin_expect(!isInside(x, y, z, dim), false))
		return;

	//const float dx = dim.x > dim.y ? (dim.x > dim.z ? dim.x : dim.z) : (dim.y > dim.z ? dim.y : dim.z) - 2;

	div[get_voxel(x, y, z, dim)] = (
		(xu[get_voxel(x + 1, y, z, dim)] - xu[get_voxel(x - 1, y, z, dim)]) +
		(xv[get_voxel(x, y + 1, z, dim)] - xv[get_voxel(x, y - 1, z, dim)]) +
		(xw[get_voxel(x, y, z + 1, dim)] - xw[get_voxel(x, y, z - 1, dim)])) / 2.0f;
	//	(xw[get_voxel(x, y, z + 1, dim)] - xw[get_voxel(x, y, z - 1, dim)])) / (dx * 2.0f);
}

//TD
__global__ void lin_solveKernel(float* __restrict__ dst, float* __restrict__ neighbor_src, float* __restrict__ src, const int3 dim, float a, bool sum, bool* __restrict__ occ)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;

	if (__builtin_expect(!isInside(x, y, z, dim), false))
		return;

	int vox = get_voxel(x, y, z, dim);
	int voxU = get_voxel(x, y + 1, z, dim);
	int voxD = get_voxel(x, y - 1, z, dim);
	int voxR = get_voxel(x + 1, y, z, dim);
	int voxL = get_voxel(x - 1, y, z, dim);
	int voxF = get_voxel(x, y, z + 1, dim);
	int voxB = get_voxel(x, y, z - 1, dim);
	float src_v = src[vox];

	int neighborCount = 0;
	float temp = 0.f;

	//if (__builtin_expect(!occ[voxL], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxL];
	//}
	//if (__builtin_expect(!occ[voxR], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxR];
	//}
	//if (__builtin_expect(!occ[voxD], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxD];
	//}
	//if (__builtin_expect(!occ[voxU], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxU];
	//}
	//if (__builtin_expect(!occ[voxB], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxB];
	//}
	//if (__builtin_expect(!occ[voxF], true)) {
	//	neighborCount += 1;
	//	temp += neighbor_src[voxF];
	//}

	temp += neighbor_src[voxL];
	temp += neighbor_src[voxR];
	temp += neighbor_src[voxD];
	temp += neighbor_src[voxU];
	temp += neighbor_src[voxB];
	temp += neighbor_src[voxF];
	neighborCount = 6;
	//temp = neighbor_src[voxL] + neighbor_src[voxR] + neighbor_src[voxU] + neighbor_src[voxF] + neighbor_src[voxD] + neighbor_src[voxB];

	float c = (sum ? 1 : 0) - neighborCount * a;

	if (__builtin_expect(c == 0, false))
		return;


	temp = (src_v - a * temp) / c;

	dst[vox] = temp;
}

//TD
__global__ void applyPressureKernel(float* __restrict__ u, float* __restrict__ v, float* __restrict__ w, float* __restrict__ p, const int3 dim, bool* __restrict__ occ)
{
	volatile int x = blockDim.x * blockIdx.x + threadIdx.x;
	volatile int y = blockDim.y * blockIdx.y + threadIdx.y;
	volatile int z = blockDim.z * blockIdx.z + threadIdx.z;


	if (__builtin_expect(!isInside(x, y, z, dim), false))
		return;

	int cur = get_voxel(x, y, z, dim);

	//const float dx = dim.x > dim.y ? (dim.x > dim.z ? dim.x : dim.z) : (dim.y > dim.z ? dim.y : dim.z) - 2;

	u[cur] -= 0.5f * (p[get_voxel(x + 1, y, z, dim)] - p[get_voxel(x - 1, y, z, dim)]);
	v[cur] -= 0.5f * (p[get_voxel(x, y + 1, z, dim)] - p[get_voxel(x, y - 1, z, dim)]);
	w[cur] -= 0.5f * (p[get_voxel(x, y, z + 1, dim)] - p[get_voxel(x, y, z - 1, dim)]);
	//u[cur] -= (0.5f * dx) * (p[get_voxel(x + 1, y, z, dim)] - p[get_voxel(x - 1, y, z, dim)]);
	//v[cur] -= (0.5f * dx) * (p[get_voxel(x, y + 1, z, dim)] - p[get_voxel(x, y - 1, z, dim)]);
	//w[cur] -= (0.5f * dx) * (p[get_voxel(x, y, z + 1, dim)] - p[get_voxel(x, y, z - 1, dim)]);
}