#include <iostream>
#include <iomanip>
#include <omp.h>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

void LeastSquareDebug()
{
	for(int i = 0;i < Cells;++i)
	{
		Cell_2D *center = &CellArray[i], *neighbour = nullptr;
		for(int m = 0;m < DV_Qu;++m)
		for(int n = 0;n < DV_Qv;++n)
		{
			center->f.BarP[m][n] = -1.0;
			for(int Iface = 0;Iface < center->celltype;++Iface)
			{
			neighbour = center->Cell_C[Iface];
				if(neighbour->xc == center->xc && neighbour->yc > center->yc)
				{
					neighbour->f.BarP[m][n] = 1.5;
				}
				else if(neighbour->xc == center->xc && neighbour->yc < center->yc)
				{
					neighbour->f.BarP[m][n] = 0.0;
				}
				else if(neighbour->yc == center->yc && neighbour->xc > center->xc)
				{
					neighbour->f.BarP[m][n] = 5;
				}
				else if(neighbour->yc == center->yc && neighbour->xc < center->xc)
				{
					neighbour->f.BarP[m][n] = 2;
				}
				else
				{
					cout << "neither xc nor yc is the same;"<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
				}
			}
			double Sum_wdxdfBP = 0.0;
			double Sum_wdydfBP = 0.0;
			for(int Iface = 0;Iface < center->celltype;++Iface)
			{
				neighbour = center->Cell_C[Iface];
				Sum_wdxdfBP += center->wdx_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
				Sum_wdydfBP += center->wdy_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
			}
			center->f.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdfBP + center->LS_M[0][1]*Sum_wdydfBP;
			center->f.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdfBP + center->LS_M[1][1]*Sum_wdydfBP;
			if(CellArray[i].f.BarP_x[m][n] != 3.0*0.5*NL || CellArray[i].f.BarP_y[m][n] != 1.5*0.5*NL)
			{
				cout <<"CellArray : "<<i<<" ---------"<<endl;
				cout <<"m = "<<m <<"  "<<CellArray[i].f.BarP_x[m][n]
					<<"    "<<CellArray[i].f.BarP_y[m][n]<<endl;
				getchar();
			}
		}
	}
}
void LSCellMatrix(Cell_2D* const &Center,int k,
					const double &neighbour_xc, double const &neighbour_yc)
{
	double dx = neighbour_xc - Center->xc, dy = neighbour_yc - Center->yc;
//
	SetZero(dx);
	SetZero(dy);
//-----------------LeastSquare Left Matrix----------------------
	double DistanceP2P = dx*dx + dy*dy;
	Center->LS_M[0][0] += dx*dx/DistanceP2P;
	Center->LS_M[0][1] += dx*dy/DistanceP2P;
	Center->LS_M[1][0] += dy*dx/DistanceP2P;
	Center->LS_M[1][1] += dy*dy/DistanceP2P;
//-----------------LeastSquare Right Matrix---------------------
	Center->wdx_C[k] = dx/DistanceP2P;
	Center->wdy_C[k] = dy/DistanceP2P;
}

void InverseMatrix_2_2(double (&LS_M)[2][2])
{
	double a[2][2],A;
	A = LS_M[0][0] * LS_M[1][1] - LS_M[0][1] * LS_M[1][0];
	if(0.0 == A)
	{
		cout <<"Singular Maxtrix : " <<__FILE__<<"  "<<__LINE__<<"  "<<__func__<<endl;
		getchar();
		return;
	}
	a[0][0] = LS_M[1][1]/A;
	a[1][1] = LS_M[0][0]/A;
	a[0][1] = -LS_M[0][1]/A;
	a[1][0] = -LS_M[1][0]/A;
	for(int i = 0;i < 2;++i)
		for(int j = 0;j < 2;++j)
			LS_M[i][j] = a[i][j];
}

void Grad_LSMatrix()
{
	Cell_2D* neighbour = nullptr;
	Cell_2D* center = nullptr;
	for(int i = 0;i != Cells;++i)
	{
		center = &CellArray[i];
		for(int k = 0;k != center->celltype;++k)
		{
			neighbour = center->Cell_C[k];
			if(neighbour != nullptr)
			{
				LSCellMatrix(center,k,neighbour->xc,neighbour->yc);
			}
			else
			{
				cout << "CellArray : " << i <<"neighbour cell invalid : "<<endl;
				getchar();
			}
		}
		InverseMatrix_2_2(center->LS_M);
	}
	cout <<"LeastSquare Matrix Construction Done" << endl;
}
void Grad_VS_LS(Cell_2D *center)
{
	Cell_2D  *neighbour = nullptr;
	for(int m = 0;m < DV_Qu;++m)
	for(int n = 0;n < DV_Qv;++n)
	{
		double Sum_wdxdhBP = 0.0;
		double Sum_wdydhBP = 0.0;

		double Sum_wdxdfBP = 0.0;
		double Sum_wdydfBP = 0.0;

		double Sum_wdxdgBP = 0.0;
		double Sum_wdydgBP = 0.0;

		for(int Iface = 0;Iface < center->celltype;++Iface)
		{
			neighbour = center->Cell_C[Iface];

			#ifdef _ARK_ALLENCAHN_FLIP
			Sum_wdxdhBP += center->wdx_C[Iface]*(neighbour->h.BarP[m][n] - center->h.BarP[m][n]); 
			Sum_wdydhBP += center->wdy_C[Iface]*(neighbour->h.BarP[m][n] - center->h.BarP[m][n]);
			#endif
//
			Sum_wdxdfBP += center->wdx_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
			Sum_wdydfBP += center->wdy_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			Sum_wdxdgBP += center->wdx_C[Iface]*(neighbour->g.BarP[m][n] - center->g.BarP[m][n]);
			Sum_wdydgBP += center->wdy_C[Iface]*(neighbour->g.BarP[m][n] - center->g.BarP[m][n]);
			#endif
		}
		#ifdef _ARK_ALLENCAHN_FLIP
		center->h.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdhBP + center->LS_M[0][1]*Sum_wdydhBP;
		center->h.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdhBP + center->LS_M[1][1]*Sum_wdydhBP;
		#endif
		//
		center->f.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdfBP + center->LS_M[0][1]*Sum_wdydfBP;
		center->f.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdfBP + center->LS_M[1][1]*Sum_wdydfBP;
//isothermal flip	
		#ifndef _ARK_ISOTHERMAL_FLIP
		center->g.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdgBP + center->LS_M[0][1]*Sum_wdydgBP;
		center->g.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdgBP + center->LS_M[1][1]*Sum_wdydgBP;
		#endif
	}
}

double const

DeltaX = dx,

DeltaY = dy,

dxSq = DeltaX*DeltaX,

dySq = DeltaY*DeltaY,

_2dx = 2.0*DeltaX,

_2dy = 2.0*DeltaY;

int const

apsi = 4, 

bpsi = 1;

namespace APSI
{

int const

a = 4, 

b = 2,

c = 8,

d = 1,

_E = 6;

}

namespace BPSI
{

int const

a = 10,

b = 2, 

c = 20,

d = 1,

_E = 12;

}

namespace CPSI
{

int const

a = 1,

b = 0,

c = 2,

d = 0,

_E = 1;
}

namespace myPSI = CPSI;



void update_Phi_x(Cell_2D *cellptr)
{  
	cellptr->msq->Phi_x
	=
	(	apsi*
	    (
	      (cellptr->Cell_C[0]->msq->Phi) - (cellptr->Cell_C[2]->msq->Phi)
	    )
	  + bpsi*
	    (
	  		(cellptr->Cell_Diag[0]->msq->Phi) - (cellptr->Cell_Diag[1]->msq->Phi)
	  	+   (cellptr->Cell_Diag[3]->msq->Phi) - (cellptr->Cell_Diag[2]->msq->Phi)
	  	)
	)/(_2dx*6);
}
//
void update_Phi_y(Cell_2D *cellptr)
{
	cellptr->msq->Phi_y 
	=
	(	apsi*
	    (
	      (cellptr->Cell_C[1]->MsQ().Phi) - (cellptr->Cell_C[3]->MsQ().Phi)
	    )
	  + bpsi*
	    (
	  		(cellptr->Cell_Diag[0]->MsQ().Phi) - (cellptr->Cell_Diag[3]->MsQ().Phi)
	  	+   (cellptr->Cell_Diag[1]->MsQ().Phi) - (cellptr->Cell_Diag[2]->MsQ().Phi)
	  	)
	)/(_2dy*6);
}
double update_Phi_yy(Cell_2D *cellptr)
{	
	return
	(
	- myPSI::c * (cellptr->msq->Phi)
	- myPSI::b * (cellptr->Cell_C[0]->msq->Phi + cellptr->Cell_C[2]->msq->Phi)
	+ myPSI::a * (cellptr->Cell_C[1]->msq->Phi + cellptr->Cell_C[3]->msq->Phi)
	+ myPSI::d * 
	  (
		 cellptr->Cell_Diag[0]->msq->Phi + cellptr->Cell_Diag[1]->msq->Phi
	   + cellptr->Cell_Diag[2]->msq->Phi + cellptr->Cell_Diag[3]->msq->Phi
	  )
	)/( myPSI::_E*dySq );
}
double update_Phi_xx(Cell_2D *cellptr)
{
	return
	(
	- myPSI::c * (cellptr->msq->Phi)
	+ myPSI::a * (cellptr->Cell_C[0]->msq->Phi + cellptr->Cell_C[2]->msq->Phi)
	- myPSI::b * (cellptr->Cell_C[1]->msq->Phi + cellptr->Cell_C[3]->msq->Phi)
	+ myPSI::d * 
	  (
		 cellptr->Cell_Diag[0]->msq->Phi + cellptr->Cell_Diag[1]->msq->Phi
	   + cellptr->Cell_Diag[2]->msq->Phi + cellptr->Cell_Diag[3]->msq->Phi
	  )
	)/( myPSI::_E*dxSq );
}
void update_AbsPhi_x(Cell_2D *cellptr)
{

}
void update_AbsPhi_y(Cell_2D *cellptr)
{
	
}
void update_Laplacian_Phi(Cell_2D *cellptr)
{
	cellptr->msq->laplacianPhi
	=
	(	4*
		(
	  		(cellptr->Cell_C[0]->msq->Phi) + (cellptr->Cell_C[1]->msq->Phi)
	  	+   (cellptr->Cell_C[2]->msq->Phi) + (cellptr->Cell_C[3]->msq->Phi)
	  	)
	+	
	    (
	  		cellptr->Cell_Diag[0]->msq->Phi + cellptr->Cell_Diag[1]->msq->Phi
	  	+   cellptr->Cell_Diag[2]->msq->Phi + cellptr->Cell_Diag[3]->msq->Phi
	  	)
	-	20*(cellptr->msq->Phi)
	)/(6*MinL*MinL);
}
void Grad_Phi_6points(Cell_2D *center)
{
	update_Phi_x(center);
	update_Phi_y(center);
	update_Laplacian_Phi(center);
//	center->msq->laplacianPhi = update_Phi_xx(center) + update_Phi_yy(center);
//
	SetZero(center->msq->Phi_x);
	SetZero(center->msq->Phi_y);
}
void Grad_Phi_CD(Cell_2D *cellptr)
{
	cellptr->msq->Phi_x = 
		(cellptr->Cell_C[0]->msq->Phi - cellptr->Cell_C[2]->msq->Phi)/_2dx;
	cellptr->msq->Phi_y = 
		(cellptr->Cell_C[1]->msq->Phi - cellptr->Cell_C[3]->msq->Phi)/_2dy;
	cellptr->msq->laplacianPhi
	= 
	(
		cellptr->Cell_C[0]->msq->Phi
	  + 
	    cellptr->Cell_C[2]->msq->Phi
	  - 2*cellptr->msq->Phi
	)/dxSq 
	+
	(
		cellptr->Cell_C[1]->msq->Phi
	  + 
	    cellptr->Cell_C[3]->msq->Phi
	  - 2*cellptr->msq->Phi
	)/dySq ;
}
void Grad_Phi_LS(Cell_2D *center)
{
	Cell_2D  *neighbour = nullptr;
	double Sum_wdxdPhi = 0.0,Sum_wdydPhi = 0.0;
	for(int Iface = 0;Iface < center->celltype;++Iface)
	{
		neighbour = center->Cell_C[Iface];
//
		Sum_wdxdPhi += center->wdx_C[Iface] * ((neighbour->msq->Phi) - (center->msq->Phi));
		Sum_wdydPhi += center->wdy_C[Iface] * ((neighbour->msq->Phi) - (center->msq->Phi));
	}
	center->msq->Phi_x = center->LS_M[0][0]*Sum_wdxdPhi + center->LS_M[0][1]*Sum_wdydPhi;
	center->msq->Phi_y = center->LS_M[1][0]*Sum_wdxdPhi + center->LS_M[1][1]*Sum_wdydPhi;
	SetZero(center->msq->Phi_x);
	SetZero(center->msq->Phi_y);
	double L = sqrt(center->msq->SqPhixPhiy());
	if(L != 0)
	{
		center->msq->Phi_x /= L;
		center->msq->Phi_y /= L;
	}
//
}
void Grad_Phi_Zero(Cell_2D *center)
{
	(center->msq->Phi_x) = 0;
	(center->msq->Phi_y) = 0;
}