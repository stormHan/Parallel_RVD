
#ifndef H_CUDA_BASCIS_H
#define H_CUDA_BASCIS_H


__device__ double computeDistance(double3 p1, double3 p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x)
		+ (p1.y - p2.y) * (p1.y - p2.y)
		+ (p1.z - p2.z) * (p1.z - p2.z));
}

#endif /* H_CUDA_BASCIS_H */