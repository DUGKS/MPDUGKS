#include <iostream>
#include <iomanip>
#include <omp.h>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

double const

DeltaX = dx,

DeltaY = dy,

dxSq = DeltaX*DeltaX,

dySq = DeltaY*DeltaY;

#include "GradSchemeBasic.h"

//------------------------------------Cartesian 6points---------------------------
// void Grad_VS_6points(Cell_2D *cellptr)
// {
// 	for(int m = 0;m < DV_Qu;++m)
// 	for(int n = 0;n < DV_Qv;++n)
// 	{
// 		#ifdef _ARK_ALLENCAHN_FLIP
// 		cellptr->h.BarP_x[m][n] = 
// 		(	apsi*
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->h.BarP[m][n]) - (cellptr->Cell_C[2]->h.BarP[m][n])
// 	  	  	)
// 	  	  + bpsi*
// 	  	  	(
// 	  			(cellptr->Cell_Diag[0]->h.BarP[m][n]) - (cellptr->Cell_Diag[1]->h.BarP[m][n])
// 	  		+   (cellptr->Cell_Diag[3]->h.BarP[m][n]) - (cellptr->Cell_Diag[2]->h.BarP[m][n])
// 	  	  	)
// 		)/(_2dx*6);
// 		cellptr->h.BarP_y[m][n] =
// 		(	apsi*
// 		    (
// 		      (cellptr->Cell_C[1]->h.BarP[m][n]) - (cellptr->Cell_C[3]->h.BarP[m][n])
// 		    )
// 		  + bpsi*
// 		    (
// 		  		(cellptr->Cell_Diag[0]->h.BarP[m][n]) - (cellptr->Cell_Diag[3]->h.BarP[m][n])
// 		  	+   (cellptr->Cell_Diag[1]->h.BarP[m][n]) - (cellptr->Cell_Diag[2]->h.BarP[m][n])
// 		  	)
// 		)/(_2dy*6);
// 		#endif
// 		//
// 		cellptr->f.BarP_x[m][n] = 
// 		(	apsi*
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->f.BarP[m][n]) - (cellptr->Cell_C[2]->f.BarP[m][n])
// 	  	  	)
// 	  	  + bpsi*
// 	  	  	(
// 	  			(cellptr->Cell_Diag[0]->f.BarP[m][n]) - (cellptr->Cell_Diag[1]->f.BarP[m][n])
// 	  		+   (cellptr->Cell_Diag[3]->f.BarP[m][n]) - (cellptr->Cell_Diag[2]->f.BarP[m][n])
// 	  	  	)
// 		)/(_2dx*6);
// 		cellptr->f.BarP_y[m][n] =
// 		(	apsi*
// 		    (
// 		      (cellptr->Cell_C[1]->f.BarP[m][n]) - (cellptr->Cell_C[3]->f.BarP[m][n])
// 		    )
// 		  + bpsi*
// 		    (
// 		  		(cellptr->Cell_Diag[0]->f.BarP[m][n]) - (cellptr->Cell_Diag[3]->f.BarP[m][n])
// 		  	+   (cellptr->Cell_Diag[1]->f.BarP[m][n]) - (cellptr->Cell_Diag[2]->f.BarP[m][n])
// 		  	)
// 		)/(_2dy*6);
// //isothermal flip	
// 		#ifndef _ARK_ISOTHERMAL_FLIP
// 		cellptr->g.BarP_x[m][n] = 
// 		(	apsi*
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->g.BarP[m][n]) - (cellptr->Cell_C[2]->g.BarP[m][n])
// 	  	  	)
// 	  	  + bpsi*
// 	  	  	(
// 	  			(cellptr->Cell_Diag[0]->g.BarP[m][n]) - (cellptr->Cell_Diag[1]->g.BarP[m][n])
// 	  		+   (cellptr->Cell_Diag[3]->g.BarP[m][n]) - (cellptr->Cell_Diag[2]->g.BarP[m][n])
// 	  	  	)
// 		)/(_2dx*6);
// 		cellptr->g.BarP_y[m][n] =
// 		(	apsi*
// 		    (
// 		      (cellptr->Cell_C[1]->g.BarP[m][n]) - (cellptr->Cell_C[3]->g.BarP[m][n])
// 		    )
// 		  + bpsi*
// 		    (
// 		  		(cellptr->Cell_Diag[0]->g.BarP[m][n]) - (cellptr->Cell_Diag[3]->g.BarP[m][n])
// 		  	+   (cellptr->Cell_Diag[1]->g.BarP[m][n]) - (cellptr->Cell_Diag[2]->g.BarP[m][n])
// 		  	)
// 		)/(_2dy*6);
// 		#endif
// 	}
// }
// void Grad_VS_4points(Cell_2D *cellptr)
// {
// 	for(int m = 0;m < DV_Qu;++m)
// 	for(int n = 0;n < DV_Qv;++n)
// 	{
// 		#ifdef _ARK_ALLENCAHN_FLIP
// 		cellptr->h.BarP_x[m][n] = 
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->h.BarP[m][n]) - (cellptr->Cell_C[2]->h.BarP[m][n])
// 	  	  	)/(_2dx);
// 		cellptr->h.BarP_y[m][n] =
// 		    (
// 		      (cellptr->Cell_C[1]->h.BarP[m][n]) - (cellptr->Cell_C[3]->h.BarP[m][n])
// 		    )/(_2dy);
// 		#endif
// 		//
// 		cellptr->f.BarP_x[m][n] = 
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->f.BarP[m][n]) - (cellptr->Cell_C[2]->f.BarP[m][n])
// 	  	  	)/(_2dx);
// 		cellptr->f.BarP_y[m][n] =
// 		    (
// 		      (cellptr->Cell_C[1]->f.BarP[m][n]) - (cellptr->Cell_C[3]->f.BarP[m][n])
// 		    )/(_2dy);
// //isothermal flip	
// 		#ifndef _ARK_ISOTHERMAL_FLIP
// 		cellptr->g.BarP_x[m][n] = 
// 	  	  	(
// 	  	    (cellptr->Cell_C[0]->g.BarP[m][n]) - (cellptr->Cell_C[2]->g.BarP[m][n])
// 	  	  	)/(_2dx);
// 		cellptr->g.BarP_y[m][n] =
// 		    (
// 		      (cellptr->Cell_C[1]->g.BarP[m][n]) - (cellptr->Cell_C[3]->g.BarP[m][n])
// 		    )/(_2dy);
// 		#endif
// 	}
// }
// void Grad_VS_LS(Cell_2D *center)
// {
// 	Cell_2D  *neighbour = nullptr;
// 	for(int m = 0;m < DV_Qu;++m)
// 	for(int n = 0;n < DV_Qv;++n)
// 	{
// 		#ifdef _ARK_ALLENCAHN_FLIP
// 		double Sum_wdxdhBP = 0.0;
// 		double Sum_wdydhBP = 0.0;
// 		#endif

// 		double Sum_wdxdfBP = 0.0;
// 		double Sum_wdydfBP = 0.0;

// 		#ifndef _ARK_ISOTHERMAL_FLIP
// 		double Sum_wdxdgBP = 0.0;
// 		double Sum_wdydgBP = 0.0;
// 		#endif

// 		for(int Iface = 0;Iface < center->celltype;++Iface)
// 		{
// 			neighbour = center->Cell_C[Iface];

// 			#ifdef _ARK_ALLENCAHN_FLIP
// 			Sum_wdxdhBP += center->wdx_C[Iface]*(neighbour->h.BarP[m][n] - center->h.BarP[m][n]); 
// 			Sum_wdydhBP += center->wdy_C[Iface]*(neighbour->h.BarP[m][n] - center->h.BarP[m][n]);
// 			#endif
// //
// 			Sum_wdxdfBP += center->wdx_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
// 			Sum_wdydfBP += center->wdy_C[Iface]*(neighbour->f.BarP[m][n] - center->f.BarP[m][n]);
// //isothermal flip
// 			#ifndef _ARK_ISOTHERMAL_FLIP
// 			Sum_wdxdgBP += center->wdx_C[Iface]*(neighbour->g.BarP[m][n] - center->g.BarP[m][n]);
// 			Sum_wdydgBP += center->wdy_C[Iface]*(neighbour->g.BarP[m][n] - center->g.BarP[m][n]);
// 			#endif
// 		}
// 		#ifdef _ARK_ALLENCAHN_FLIP
// 		center->h.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdhBP + center->LS_M[0][1]*Sum_wdydhBP;
// 		center->h.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdhBP + center->LS_M[1][1]*Sum_wdydhBP;
// 		#endif
// 		//
// 		center->f.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdfBP + center->LS_M[0][1]*Sum_wdydfBP;
// 		center->f.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdfBP + center->LS_M[1][1]*Sum_wdydfBP;
// //isothermal flip	
// 		#ifndef _ARK_ISOTHERMAL_FLIP
// 		center->g.BarP_x[m][n] = center->LS_M[0][0]*Sum_wdxdgBP + center->LS_M[0][1]*Sum_wdydgBP;
// 		center->g.BarP_y[m][n] = center->LS_M[1][0]*Sum_wdxdgBP + center->LS_M[1][1]*Sum_wdydgBP;
// 		#endif
// 	}
// }

void Grad_VS_4points(Cell_2D *cellptr)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	update_DVDF_Grad4points(cellptr,&Cell_2D::h);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	update_DVDF_Grad4points(cellptr,&Cell_2D::f);
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP
	update_DVDF_Grad4points(cellptr,&Cell_2D::g);
	#endif
}
void Grad_VS_6points(Cell_2D *cellptr)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	update_DVDF_Grad6points(cellptr,&Cell_2D::h);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	update_DVDF_Grad6points(cellptr,&Cell_2D::f);
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP
	update_DVDF_Grad6points(cellptr,&Cell_2D::g);
	#endif
}
void Grad_VS_LS(Cell_2D *center)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	update_DVDF_Grad_LS(center,&Cell_2D::h);
	#endif
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	update_DVDF_Grad_LS(center,&Cell_2D::f);
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP
	update_DVDF_Grad_LS(center,&Cell_2D::g);
	#endif
}
//!----------------------used to calculate pseudopotetial forece-------------------
void update_Psi_x(Cell_2D *cellptr)
{
	cellptr->msq->Fx = update_MQ_x(cellptr,&MacroQuantity::Psi);
}
void update_Psi_y(Cell_2D *cellptr)
{
	cellptr->msq->Fy = update_MQ_y(cellptr,&MacroQuantity::Psi);
}
void update_Psi_x(Face_2D *faceptr)
{
	faceptr->msqh->Fx = update_MQ_x(faceptr,&MacroQuantity::Psi);
}
void update_Psi_y(Face_2D *faceptr)
{
	faceptr->msqh->Fy = update_MQ_y(faceptr,&MacroQuantity::Psi);
}
//!
void update_Phi_x(Cell_2D *cellptr)
{  
	cellptr->msq->Phi_x = update_MQ_x(cellptr,&MacroQuantity::Phi);
}
//
void update_Phi_y(Cell_2D *cellptr)
{
	cellptr->msq->Phi_y = update_MQ_y(cellptr,&MacroQuantity::Phi);
}
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
// void update_AbsPhi_x(Cell_2D *cellptr)
// {

// }
// void update_AbsPhi_y(Cell_2D *cellptr)
// {
	
// }
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