#include <iostream>
#include <iomanip>
#include <omp.h>
#include "DUGKSDeclaration.h"
using std::ios;
using std::setiosflags;
using std::setprecision;
using std::cout;
using std::endl;

double const aTP = 4.0/3.0, bTP = 1.0/3.0;

int const ThreadNum = omp_get_max_threads();

double ResidualPer1k = 1.0;

double const PRECISION = 1.0E-16;

double const XDEBUG = 13.25,YDEBUG = 0.25;


//--------------------------DEBUG-------------------------------

extern void Output_L2Norm(double const &t,double &L2_uv, double &L2_p);

extern void Output_Flowfield(double const &t,int step);

extern void Output_SumRho(double t);

extern void Output_Residual(double t,double Residual);

extern void Output_MidX(int step);

//------------------------Boundary.cpp----------------------

extern void P_Inlet_4_Boundary();

extern void P_Outlet_5_Boundary();

extern void WallShadowC_fBP(Cell_2D &shadowCell);

extern void Wall_3_DS(Face_2D &face);

extern void Wall_3_NEE(Face_2D &face);

extern void Wall_3_BB(Face_2D &face);

extern void fluxCheck(Face_2D const* faceptr);

//----------------------------------------------------------

extern void Update_force(Cell_2D &cell);

extern void Update_force_h(Face_2D &face);

extern void MacroSource(Cell_2D *cellptr);

extern void Output_xcyc();

//----------------------------------DEBUG---------------------------------------
#ifdef _ZERO_NDEBUG_FLIP

extern void Output_UVP(double const &t);

extern void Output_fBh(Face_2D& face,double t);

extern void Output_gBh(Face_2D& face,double t);

extern void Output_f_hdt(Face_2D& face,double t);

extern void Output_g_hdt(Face_2D& face,double t);

extern void Output_fT(Cell_2D& face,double t);

extern void Output_gT(Cell_2D& face,double t);

extern void Output_fT_Append(Cell_2D &cell,double dt);

extern void Output_gT_Append(Cell_2D &cell,double dt);

extern void Output_f_hdt_Append(Face_2D &face,double dt);

extern void Output_g_hdt_Append(Face_2D &face,double dt);

// extern void Output_fFlux_Append(Cell_2D &cell,double dt);

// extern void Output_gFlux_Append(Cell_2D &cell,double dt);

extern void Output_phi_Bh(Face_2D &face,double t);

#endif

//-------------------------------------GradScheme.cpp---------------------------
extern void LeastSquareDebug();

extern void Grad_VS_LS(Cell_2D *center);

extern void Grad_VS_6points(Cell_2D *center);

extern void Grad_VS_4points(Cell_2D *center);
//------------------------------------------------------------------------------

void Update_phi_BP(Cell_2D& cell);

void Update_phi_h(Face_2D& face);

void Update_phiFlux_h(Face_2D& face);

void Update_Flux(Face_2D &face);

void Update_BoundFlux(Face_2D &face);

void Update_phi_T(Cell_2D& cell);

void Zero_PartialDerivatives(Cell_2D& cell);

void Flux_2D(Face_2D &face);

void Flux_2D_Limiter(Face_2D &face);

void Update_Residual(int step);

void UpdateL2Error(int step);

void Update_SumRho(int step);
//--------------------------------DmVn.cpp--------------------------
extern void Update_DVDF_Source(Cell_2D &cell);

extern void Update_DVDF_Source_h(Face_2D &face);

extern void IntegralShearStress();

//-------------------------------------------------------------------
 
void DUGKS2DSolver()
{

Output_Flowfield(0.0,0);
//Output_xcyc();
cout << "iteration start    ThreadNum : "<<ThreadNum<<endl;
omp_set_num_threads(ThreadNum);
step = 0;
#pragma omp parallel
{
//while(step < End_Step)
while(ResidualPer1k > RESIDUAL)
{
//---------------------------------------------------------
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		Update_phi_Eq(CellArray[n]);
		#ifdef _ARK_FORCE_FLIP
		Update_DVDF_Source(CellArray[n]);
		#endif
		Update_phi_BP(CellArray[n]);
	}
//-------------------Update-shadow fBP--------------------------------
	#ifdef _Wall_3_BCs_FLIP
	#pragma omp for schedule(guided)
	for(int n = 0;n < WallFaceNum;++n)
		WallShadowC_fBP(WallShadowCA[n]);
	#endif
	#ifdef _P_INLET_4_BCS_FLIP	
		P_Inlet_4_Boundary();
	#endif
	#ifdef _P_OUTLET_5_BCS_FLIP
		P_Outlet_5_Boundary();
	#endif
//-------------------------------Update Grad DVDF_BarPlus-------------------------------
	#ifdef _FLUX_SCHEME_UW_ARK
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
	//Grad_VS_LS(&CellArray[n]);
	Grad_VS_6points(&CellArray[n]);
	//Grad_VS_4points(&CellArray[n]);
	//Zero_PartialDerivatives(CellArray[n]);
	}
	#endif
//-------------------------------Flux-------------------------------------
//-------------------------------Interior Face-----------------------------	
	#ifdef _ARK_LIMITER_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D_Limiter(*InteriorFaceA[n]);
		}
		#pragma omp for schedule(guided)
		for(int n = 0;n < BoundFaceNum;++n)
			Flux_2D(*BoundFaceA[n]);
	#else
		#pragma omp for schedule(guided)
		for(int n = 0;n < InteriorFaceNum;++n)
		{
			Flux_2D(*InteriorFaceA[n]);
		}
		#ifdef _PERIODIC_12_8_BCs_FLIP
		#pragma omp for schedule(guided)
		for(int k = 0;k < PeriodicFaceNum;++k)
		{
			Flux_2D(*PeriodicFaceA[k]);
		}
		#endif
	#endif
//
	#ifdef _Wall_3_BCs_FLIP
		#pragma omp for schedule(guided)
		for(int n = 0;n < WallFaceNum;++n)
		{
			#ifdef _Wall_3_BCs_DS
			Wall_3_DS(*WallFaceA[n]);
			#endif
//
			#ifdef _Wall_3_BCs_NEE
			Wall_3_NEE(*WallFaceA[n]);
			#endif
//
			#ifdef _Wall_3_BCs_BB
			Flux_2D(*WallFaceA[n]);
			Wall_3_BB(*WallFaceA[n]);
			#endif
		}
	#endif
// //----------------------------------------Wall Face----------------------------------
// 	#ifdef _Wall_3_BCs_FLIP
// // 	for(int i_Wall = 0;i_Wall < WallFaceNum;++i_Wall)
// // 	{
		
// // 		#ifdef _NEE_BOUNDARY_SCHEME_FLIP
// // 		NonEquilibriumExtrapolation(*WallFaceA[i_Wall]);
// // 		#endif
// // //		
// // 		//DiffusiveScatter(*WallFaceA[i_Wall]);
// // //		
// // 		#ifdef _BB_BOUNDARY_SCHEME_FLIP
// // 		BounceBack(*WallFaceA[i_Wall]);
// // 		#endif
// // 	}
// 	for(int i = 0;i < WallFaceNum;++i)
// 		Update_BoundFlux(*WallFaceA[i]);
// 	#endif
	#pragma omp for schedule(guided) 
	LoopPS(Faces)
	{
		Update_phiFlux_h(FaceArray[n]);
	}
	#pragma omp for schedule(guided)
		for(int n = 0;n < WallFaceNum;++n)
		{
			fluxCheck(WallFaceA[n]);
		}
//--------------------auxilary DDF : phi tilde----------------------------
	#pragma omp for schedule(guided)
	LoopPS(Cells)
		Update_phi_T(CellArray[n]);
//---------------------------------macro variables---------------------
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		Update_MacroVar(CellArray[n]);
	}
//-------------------------pseudopotential force-------------------------
	#ifdef _ARK_ALLENCAHN_FLIP
	#pragma omp for schedule(guided)
	LoopPS(Cells)
	{
		MacroSource(&CellArray[n]);
	}
	#endif
//
	#pragma omp single
	{
		++step;
		if(step%ConvergenceControl == 0)
		{
			Update_SumRho(step);
			#ifdef _OUTPUT_L2NORM_ERROR_FLIP
			UpdateL2Error(step);
			#endif
		}
		if(step%ResidualControl == 0)
		{
			Update_Residual(step);
		}
	}
}
}
	IntegralShearStress();
	Output_Flowfield((step)*dt,step);
	Output_MidX(step);
}

void Zero_PartialDerivatives(Cell_2D &cell)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		cell.f.BarP_x[i][j] = 0.0;
		cell.f.BarP_y[i][j] = 0.0;
//isothermal flip
#ifndef _ARK_ISOTHERMAL_FLIP
		cell.g.BarP_x[i][j] = 0.0;
		cell.g.BarP_y[i][j] = 0.0;
#endif
	}
}
void Update_phi_BP(Cell_2D& cell)
{
	// for(int i = 0;i < DV_Qu;++i)
	// for(int j = 0;j < DV_Qv;++j)
	LoopVS(DV_Qu,DV_Qv)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.BarP[i][j] = cell.h.aBP*cell.h.Tilde[i][j] + cell.h.bBP*cell.h.Eq[i][j]
						  + cell.h.cBP*cell.h.So[i][j];
		#endif

		cell.f.BarP[i][j] = cell.f.aBP*cell.f.Tilde[i][j] + cell.f.bBP*cell.f.Eq[i][j];
		#ifdef _ARK_FORCE_FLIP
		cell.f.BarP[i][j] += cell.f.cBP*cell.f.So[i][j];
		#endif
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		cell.g.BarP[i][j] = cell.g.aBP*cell.g.Tilde[i][j] + cell.g.bBP*cell.g.Eq[i][j];
		#endif
	}
}
void Update_phi_h(Face_2D& face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[i][j] = face.h.ah*face.h.BhDt[i][j] + face.h.bh*face.h.EqhDt[i][j]
						 + face.h.ch*face.h.SohDt[i][j];
		#endif

		face.f.hDt[i][j] = face.f.ah*face.f.BhDt[i][j] + face.f.bh*face.f.EqhDt[i][j];
		#ifdef _ARK_FORCE_FLIP
		face.f.hDt[i][j] += face.f.ch*face.f.SohDt[i][j];
		#endif
		//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		face.g.hDt[i][j] = face.f.ah*face.g.BhDt[i][j] + face.f.bh*face.g.EqhDt[i][j];
		#endif
	}
}
void Update_phiFlux_h(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[i][j] *= face.xi_n_dS[i][j];
		#endif
		//
		face.f.hDt[i][j] *= face.xi_n_dS[i][j];
		//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		face.g.hDt[i][j] *= face.xi_n_dS[i][j];
		#endif
	}
}
inline
double VenkatakrishnanExpression(double a,double b)
{
	return (a*a + 2*a*b + Kh3)/(a*a + 2*b*b + a*b + Kh3);
} 
void VenkatakrishnanFluxLimiter(Cell_2D &cell,int const &i,int const &j)
{
	double GradfBPDotDelta[4],LfBP[4];
	double MaxfBP = cell.f.BarP[i][j],MinfBP = cell.f.BarP[i][j];
	double MinfBPLimiter = 1;
//
	// double MaxfBP = -1E+5,MinfBP = -MaxfBP, MaxgBP = -1E+5,MingBP = -MaxgBP;
	#ifndef _ARK_ISOTHERMAL_FLIP
	double GradgBPDotDelta[4],LgBP[4];
	double MaxgBP = cell.g.BarP[i][j],MingBP = cell.g.BarP[i][j];
	double MingBPLimiter = 1;
	#endif
//	
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		if(cell.Cell_C[iFace]->f.BarP[i][j] > MaxfBP)
			MaxfBP = cell.Cell_C[iFace]->f.BarP[i][j];
		if(cell.Cell_C[iFace]->f.BarP[i][j] < MinfBP)
			MinfBP = cell.Cell_C[iFace]->f.BarP[i][j];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		if(cell.Cell_C[iFace]->g.BarP[i][j] > MaxgBP)
			MaxgBP = cell.Cell_C[iFace]->g.BarP[i][j];
		if(cell.Cell_C[iFace]->g.BarP[i][j] < MingBP)
			MingBP = cell.Cell_C[iFace]->g.BarP[i][j];
#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		GradfBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.f.BarP_x[i][j]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.f.BarP_y[i][j];					        
		if(GradfBPDotDelta[iFace] > 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MaxfBP-cell.f.BarP[i][j],GradfBPDotDelta[iFace]);
		}
		else if(GradfBPDotDelta[iFace] < 0)
		{
			LfBP[iFace]
			=
			VenkatakrishnanExpression(MinfBP-cell.f.BarP[i][j],GradfBPDotDelta[iFace]);
		}
		else
			LfBP[iFace] = 1;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		GradgBPDotDelta[iFace] = (cell.Face_C[iFace]->xf - cell.xc)*cell.g.BarP_x[i][j]
						        + (cell.Face_C[iFace]->yf - cell.yc)*cell.g.BarP_y[i][j];
		if(GradgBPDotDelta[iFace] > 0)
		{
			LgBP[iFace]
			=  
			VenkatakrishnanExpression(MaxgBP-cell.g.BarP[i][j],GradgBPDotDelta[iFace]);
		}
		else if(GradgBPDotDelta[iFace] < 0)
		{
			LgBP[iFace]
			=
			VenkatakrishnanExpression(MingBP-cell.g.BarP[i][j],GradgBPDotDelta[iFace]);
		}
		else
			LgBP[iFace] = 1;
#endif
	}
	for(int iFace = 0;iFace < cell.celltype;++iFace)
	{
		if(LfBP[iFace] < MinfBPLimiter) MinfBPLimiter = LfBP[iFace];
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
		if(LgBP[iFace] < MingBPLimiter) MingBPLimiter = LgBP[iFace];
#endif
	}
	cell.fBPLimiter = MinfBPLimiter;
//isothermal
#ifndef _ARK_ISOTHERMAL_FLIP
	cell.gBPLimiter = MingBPLimiter;
#endif
}
void UW_Interior_phi_Bh_Limiter(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - hDt*(xi_u[QuIndex]) - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[j] - ptr_C->yc;
	VenkatakrishnanFluxLimiter(*ptr_C,i,j);
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[i][j] = ptr_C->h.BarP[i][j]
					  + (dx*ptr_C->h.BarP_x[i][j] + dy*ptr_C->h.BarP_y[i][j]);
	#endif
	//
	face.f.BhDt[i][j] = ptr_C->f.BarP[i][j]
					  + ptr_C->fBPLimiter*(dx*ptr_C->f.BarP_x[i][j] + dy*ptr_C->f.BarP_y[i][j]);
//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP
	face.g.BhDt[i][j] = ptr_C->g.BarP[i][j]
					  + ptr_C->gBPLimiter*(dx*ptr_C->g.BarP_x[i][j] + dy*ptr_C->g.BarP_y[i][j]);
	#endif
}
void UW_Interior_phi_Bh(Face_2D& face,Cell_2D* ptr_C,int const &i,int const &j)
{
	double dx = face.xf - hDt*(xi_u[QuIndex]) - ptr_C->xc;
	double dy = face.yf - hDt*xi_v[j] - ptr_C->yc;
	//
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[i][j] = ptr_C->h.BarP[i][j]
					  + (dx*ptr_C->h.BarP_x[i][j] + dy*ptr_C->h.BarP_y[i][j]);
	#endif
	//				  
	face.f.BhDt[i][j] = ptr_C->f.BarP[i][j]
					  + (dx*ptr_C->f.BarP_x[i][j] + dy*ptr_C->f.BarP_y[i][j]);
	//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP
	face.g.BhDt[i][j] = ptr_C->g.BarP[i][j]
					  + (dx*ptr_C->g.BarP_x[i][j] + dy*ptr_C->g.BarP_y[i][j]);
	#endif
}
void UW_Interior_phi_Bh(Face_2D& face,int const &i,int const &j)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	face.h.BhDt[i][j] = 0.5*(face.owner->h.BarP[i][j]+face.neigh->h.BarP[i][j]);
	#endif
	//				  
	face.f.BhDt[i][j] = 0.5*(face.owner->f.BarP[i][j]+face.neigh->f.BarP[i][j]);
	//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP
	face.g.BhDt[i][j] = 0.5*(face.owner->g.BarP[i][j]+face.neigh->g.BarP[i][j]);
	#endif
}
void CD_Interior_phi_Bh(Face_2D &face,int i,int j)
{
	#ifdef _ARK_ALLENCAHN_FLIP
	double _hBP_xF = 1000,_hBP_yF = 1000;
	#endif
	double _fBP_xF = 1000,_fBP_yF = 1000;
	#ifndef _ARK_ISOTHERMAL_FLIP
	double _gBP_xF = 1000,_gBP_yF = 1000;
	#endif
		if(EqualZero(face.Vx))
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			//_hBP_xF = 0.5*(face.owner->h.BarP_x[i][j] + face.neigh->h.BarP_x[i][j]);
			_hBP_xF = (0.5*(face.faceCells[3]->h.BarP[i][j] - face.faceCells[1]->h.BarP[i][j])
					+ 0.5*(face.faceCells[2]->h.BarP[i][j] - face.faceCells[0]->h.BarP[i][j]))
					/face._dx;
			_hBP_yF = (face.owner->h.BarP[i][j] - face.neigh->h.BarP[i][j])/face._dy;
			#endif
			//
			//_fBP_xF =  0.5*(face.owner->f.BarP_x[i][j] + face.neigh->f.BarP_x[i][j]);
			_fBP_xF = (0.5*(face.faceCells[3]->f.BarP[i][j] - face.faceCells[1]->f.BarP[i][j])
					+ 0.5*(face.faceCells[2]->f.BarP[i][j] - face.faceCells[0]->f.BarP[i][j]))
					/face._dx;
			_fBP_yF = (face.owner->f.BarP[i][j] - face.neigh->f.BarP[i][j])/face._dy;
			//
			#ifndef _ARK_ISOTHERMAL_FLIP
			//_gBP_xF =  0.5*(face.owner->g.BarP_x[i][j] + face.neigh->g.BarP_x[i][j]);
			_gBP_xF = (0.5*(face.faceCells[3]->g.BarP[i][j] - face.faceCells[1]->g.BarP[i][j])
					+ 0.5*(face.faceCells[2]->g.BarP[i][j] - face.faceCells[0]->g.BarP[i][j]))
					/face._dx;
			_gBP_yF = (face.owner->g.BarP[i][j] - face.neigh->g.BarP[i][j])/face._dy;
			#endif
		}
		else if(EqualZero(face.Vy))
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			//_hBP_yF = 0.5*(face.owner->h.BarP_y[i][j] + face.neigh->h.BarP_y[i][j]);
			_hBP_yF = (0.5*(face.faceCells[3]->h.BarP[i][j] - face.faceCells[1]->h.BarP[i][j])
					+ 0.5*(face.faceCells[2]->h.BarP[i][j] - face.faceCells[0]->h.BarP[i][j]))
					/face._dy;
			_hBP_xF = (face.owner->h.BarP[i][j] - face.neigh->h.BarP[i][j])/face._dx;
			#endif
			//
			//_fBP_yF = 0.5*(face.owner->f.BarP_y[i][j] + face.neigh->f.BarP_y[i][j]);
			_fBP_yF = (0.5*(face.faceCells[3]->f.BarP[i][j] - face.faceCells[1]->f.BarP[i][j])
					+ 0.5*(face.faceCells[2]->f.BarP[i][j] - face.faceCells[0]->f.BarP[i][j]))
					/face._dy;
			_fBP_xF = (face.owner->f.BarP[i][j] - face.neigh->f.BarP[i][j])/face._dx;
			//
			#ifndef _ARK_ISOTHERMAL_FLIP
			//_gBP_yF = 0.5*(face.owner->g.BarP_y[i][j] + face.neigh->g.BarP_y[i][j]);
			_gBP_yF = (0.5*(face.faceCells[3]->g.BarP[i][j] - face.faceCells[1]->g.BarP[i][j])
					+ 0.5*(face.faceCells[2]->g.BarP[i][j] - face.faceCells[0]->g.BarP[i][j]))
					/face._dy;
			_gBP_xF = (face.owner->g.BarP[i][j] - face.neigh->g.BarP[i][j])/face._dx;
			#endif
		}
		else
		{
			cout <<"wrong interface"<<endl;
			getchar();
		}
		#ifdef _ARK_ALLENCAHN_FLIP
		face.h.BhDt[i][j] = 0.5*(face.owner->h.BarP[i][j] + face.neigh->h.BarP[i][j])
				- hDt*(_hBP_xF*(xi_u[QuIndex]) + _hBP_yF*xi_v[j]);
		#endif
		//
		face.f.BhDt[i][j] = 0.5*(face.owner->f.BarP[i][j] + face.neigh->f.BarP[i][j])
				- hDt*(_fBP_xF*(xi_u[QuIndex]) + _fBP_yF*xi_v[j]);
		//
		#ifndef _ARK_ISOTHERMAL_FLIP
		face.g.BhDt[i][j] = 0.5*(face.owner->g.BarP[i][j] + face.neigh->g.BarP[i][j])
				- hDt*(_gBP_xF*(xi_u[QuIndex]) + _gBP_yF*xi_v[j]);
		#endif	
}
void Update_phi_fBh(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#if defined _FLUX_SCHEME_CD_ARK
			CD_Interior_phi_Bh(face,i,j);
		#elif defined _FLUX_SCHEME_UW_ARK
		if(face.xi_n_dS[i][j] > 0)
			UW_Interior_phi_Bh(face,face.owner,i,j);
		else if(face.xi_n_dS[i][j] < 0)
			UW_Interior_phi_Bh(face,face.neigh,i,j);
		else
			UW_Interior_phi_Bh(face,i,j);
		#else 
			exit(-1);
		#endif
	}
}
void Update_phi_fBh_Limiter(Face_2D &face)
{
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		#if defined _FLUX_SCHEME_CD_ARK
			CD_Interior_phi_Bh(face,i,j);
		#elif defined _FLUX_SCHEME_UW_ARK	
		if(face.xi_n_dS[i][j] >= 0)
			UW_Interior_phi_Bh_Limiter(face,face.owner,i,j);
		else
			UW_Interior_phi_Bh_Limiter(face,face.neigh,i,j);
		#else 
			exit(-1);
		#endif
	}
}
void Flux_2D_Limiter(Face_2D &face)
{	
	Update_phi_fBh_Limiter(face);
	Update_MacroVar_h(face);
	Update_phi_Eqh(face);
	Update_phi_h(face);
	Update_phiFlux_h(face);
//------------------------------------------DEBUG---------------------------------
}
void Flux_2D(Face_2D &face)
{	
	Update_phi_fBh(face);
	Update_MacroVar_h(face);
	Update_phi_Eqh(face);
	#ifdef _ARK_FORCE_FLIP
	Update_DVDF_Source_h(face);
	#endif
	Update_phi_h(face);
//	Update_phiFlux_h(face);
}
//-------------------------------------------------------------------------------
void Update_phi_T(Cell_2D &cell)
{
	// for(int i = 0;i < DV_Qu;++i)
	// for(int j = 0;j < DV_Qv;++j)
	LoopVS(DV_Qu,DV_Qv)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		double hFluxSum = 0.0;
		#endif
		double fFluxSum = 0;
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		gFluxSum = 0;
		#endif
		for(int k = 0;k < cell.celltype;++k)
		{
			#ifdef _ARK_ALLENCAHN_FLIP
			hFluxSum += cell.signFlux[k]*cell.Face_C[k]->h.hDt[i][j];
			#endif
			//
			fFluxSum += cell.signFlux[k]*cell.Face_C[k]->f.hDt[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			gFluxSum += cell.signFlux[k]*cell.Face_C[k]->g.hDt[i][j];
			#endif
		}
		#ifdef _ARK_ALLENCAHN_FLIP
		cell.h.Tilde[i][j] = aTP*cell.h.BarP[i][j] - bTP*cell.h.Tilde[i][j]
							+ cell.DtSlashVolume*hFluxSum;
		#endif
		//
		cell.f.Tilde[i][j] = aTP*cell.f.BarP[i][j] - bTP*cell.f.Tilde[i][j]
							+ cell.DtSlashVolume*fFluxSum;
//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP
		cell.g.Tilde[i][j] = aTP*cell.g.BarP[i][j] - bTP*cell.g.Tilde[i][j]
							+ cell.DtSlashVolume*gFluxSum;
		#endif
	}
}
void Update_SumRho(int step)
{
	SumRho = 0.0;
	SumT = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		SumRho += CellArray[i].MsQ().Rho;
		#ifndef _ARK_ISOTHERMAL_FLIP 
		SumT += CellArray[i].MsQ().T;
		#endif
	}
	Output_SumRho(step*dt);
}
void Update_Residual(int step)
{
//---------------------------density residual-------------------------
		// double SumRho = 0.0,SumdRho = 0.0;
		// double dRho = 0.0;
		// LoopPS(Cells)
		// {
		// 	dRho = CellArray[n].MsQ().Rho - CellArray[n].MsQ().Rho_1k;
		// 	SumdRho += dRho*dRho;
		// 	SumRho += CellArray[n].MsQ().Rho*CellArray[n].MsQ().Rho;
		// 	CellArray[n].MsQ().Rho_1k = CellArray[n].MsQ().Rho;
		// }
		// ResidualPer1k = sqrt(SumdRho/(SumRho + 1.0E-30));
//
//---------------------------velocity residual--------------------------
		double SumUV = 0.0,Sumdudv = 0.0;
		double du = 0.0, dv = 0.0;
		for(int i = 0;i < Cells;++i)
		{
			du = CellArray[i].MsQ().U - CellArray[i].MsQ().U_1k;
			dv = CellArray[i].MsQ().V - CellArray[i].MsQ().V_1k;
			Sumdudv += du*du + dv*dv;
			SumUV += CellArray[i].MsQ().SqUV();
			CellArray[i].MsQ().U_1k = CellArray[i].MsQ().U;
			CellArray[i].MsQ().V_1k = CellArray[i].MsQ().V;
		}
		ResidualPer1k = sqrt(Sumdudv/(SumUV + 1.0E-30));
		Output_Residual(step*dt,ResidualPer1k);
		if(step%writeFileControl == 0)
		{
			Output_Flowfield(step*dt, step);
			Output_MidX(step);
		}
		//Output_xcyc();
	#ifndef _ARK_NOHUP_FLIP
		cout << setiosflags(ios::scientific) << setprecision(12);
		cout <<step <<"    "<<step*dt<<"    "<<SumRho<<"    "<<SumT<<"    "<<ResidualPer1k<<'\n';
	#endif
}
void UpdateL2Error(int step)
{
	double L2_uv, L2_p=0;
	Output_L2Norm(step*dt,L2_uv,L2_p);
}