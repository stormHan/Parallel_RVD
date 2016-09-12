
#ifndef H_CUDA_BASCIS_H
#define H_CUDA_BASCIS_H


__device__ double computeDistance(double3 p1, double3 p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x)
		+ (p1.y - p2.y) * (p1.y - p2.y)
		+ (p1.z - p2.z) * (p1.z - p2.z));
}



__device__
double3 add(double3 a, double3 b)
{
	a.x = a.x + b.x;
	a.y = a.y + b.y;
	a.z = a.z + b.z;

	return a;
}

__device__
double3 sub(double3 a, double3 b)
{
	a.x = a.x - b.x;
	a.y = a.y - b.y;
	a.z = a.z - b.z;

	return a;
}

__device__
double dot(double3 a, double3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ 
int sgn(double _d)
{
	return (_d > 0) ? 1 : (
		(_d < 0) ? -1 : 0
		);
}

__device__
double m_fabs(double x)
{
	if (x > 0) return x;
	else return -x;
}
#endif /* H_CUDA_BASCIS_H */