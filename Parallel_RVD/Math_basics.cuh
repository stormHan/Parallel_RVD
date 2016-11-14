
#ifndef H_CUDA_BASCIS_H
#define H_CUDA_BASCIS_H


__device__ double computeDistance(double3 p1, double3 p2)
{
	return sqrt((p1.x - p2.x) * (p1.x - p2.x)
		+ (p1.y - p2.y) * (p1.y - p2.y)
		+ (p1.z - p2.z) * (p1.z - p2.z));
}

__device__
inline double sqr_d(double x)
{
	return x * x;
}

__device__ 
inline double distance2(double3 p1, double3 p2)
{
	return sqr_d(p2.x - p1.x)
		+ sqr_d(p2.y - p1.y) + sqr_d(p2.z - p1.z);
}
__device__
inline double distance2(double* p1, double* p2, int dim)
{
	double result = 0.0;
	for (int i = 0; i < dim; ++i)
	{
		result += sqr_d(p2[i] - p1[i]);
	}
	return result;
}

__device__
inline double maxd(double d1, double d2)
{
	return (d1 > d2) ? d1 : d2;
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

__device__ 
double computeTriangleArea(double3 p1, double3 p2, double3 p3)
{
	double a = computeDistance(p1, p2);
	double b = computeDistance(p2, p3);
	double c = computeDistance(p1, p3);

	//Heron's Formula to compute the area of the triangle
	double p = (a + b + c) / 2;
	if (p >= a && p >= b && p >= c)
		return sqrt(p * (p - a) * (p - b) * (p - c));
	else
		return 0.0;
}

__device__
void computeTriangleCentriod(
double3 p, double3 q, double3 r, double a, double b, double c,
double3& Vg, double& V
)
{
	double abc = a + b + c;
	double area = computeTriangleArea(p, q, r);
	V = area / 3.0 * abc;

	double wp = a + abc;
	double wq = b + abc;
	double wr = c + abc;

	double s = area / 12.0;
	Vg.x = s * (wp * p.x + wq * q.x + wr * r.x);
	Vg.y = s * (wp * p.y + wq * q.y + wr * r.y);
	Vg.z = s * (wp * p.z + wq * q.z + wr * r.z);
}

#endif /* H_CUDA_BASCIS_H */