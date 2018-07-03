#include "Mesh_2D.h"

double TriArea(double Xbeg, double Ybeg, double X_A, double Y_A,double X_B, double Y_B)
{
	double Area = 0.5 * ((X_A - Xbeg) * (Y_B - Ybeg) - (X_B - Xbeg) * (Y_A - Ybeg));
	return (Area > 0 ? Area : -Area);
}
void AllocateARK(double** &f,int const Qu,int const Qv)
{
	f = new double* [Qu];
	for(int i = 0;i < Qu;++i)
		f[i] = new double[Qv];
}
void DeallocateARK(double** &f,int const Qu,int const Qv)
{
	for(int i = 0;i < Qu;++i)
	{
		delete[] f[i];
		f[i] = nullptr;
	}
	delete[] f;
	f = nullptr;
}