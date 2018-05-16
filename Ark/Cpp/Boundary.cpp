#include "DUGKSDeclaration.h"

void P_Inlet_4_Boundary()
{
	for(int k = 0;k < P_InletFaceNum;++k)
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_InletShadowCA[k].f.BarP[i][j] = P_InletShadowCA[k].Cell_C[0]->f.BarP[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			P_InletShadowCA[k].g.BarP[i][j] = P_InletShadowCA[k].Cell_C[0]->g.BarP[i][j];
			#endif
		}
}
void P_Outlet_5_Boundary()
{
	for(int k = 0;k < P_OutletFaceNum;++k)
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			P_OutletShadowCA[k].f.BarP[i][j] = P_OutletShadowCA[k].Cell_C[0]->f.BarP[i][j];
//isothermal flip
			#ifndef _ARK_ISOTHERMAL_FLIP
			P_OutletShadowCA[k].g.BarP[i][j] = P_OutletShadowCA[k].Cell_C[0]->g.BarP[i][j];
			#endif
		}
}
void ExtrapolationfBP(Cell_2D *shadowCell,Cell_2D const*neighb,Cell_2D const*nextNeighb)
{
	shadowCell->MsQ().Phi = 2*neighb->MsQ().Phi - nextNeighb->MsQ().Phi;
	shadowCell->MsQ().Rho = 2*neighb->MsQ().Rho - nextNeighb->MsQ().Rho;
	LoopVS(DV_Qu,DV_Qv)
	{
		#ifdef _ARK_ALLENCAHN_FLIP
		shadowCell->h.BarP[i][j] = 2*neighb->h.BarP[i][j] - nextNeighb->h.BarP[i][j];
		#endif
		//
		shadowCell->f.BarP[i][j] = 2*neighb->f.BarP[i][j] - nextNeighb->f.BarP[i][j];
		//
		#ifndef _ARK_ISOTHERMAL_FLIP
		shadowCell->g.BarP[i][j] = 2*neighb->g.BarP[i][j] - nextNeighb->g.BarP[i][j];
		#endif
	}
	
}
void WallShadowC_fBP(Cell_2D &shadowCell)
{
// used for GradfBP
	int const TOP = top,BOTTOM = bottom,LEFT = left,RIGHT = right;
	#ifdef _CARTESIAN_MESH_FLIP
		if(TOP == shadowCell.zone)
		{
			ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[3]);
		}
		else if(BOTTOM == shadowCell.zone)
		{
			ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[1]);
		}
		else if(RIGHT == shadowCell.zone)
		{
			ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[2]);
		}
		else if(LEFT == shadowCell.zone)
		{
			ExtrapolationfBP(&shadowCell,shadowCell.ShadowC,shadowCell.ShadowC->Cell_C[0]);
		}
	#endif
		//
		// #ifdef _ARK_ALLENCAHN_FLIP
		// shadowCell.MsQ()= cell->MsQ();
		// #endif
		// //
		// shadowCell.MsQ().U = -cell->MsQ().U;//Non-Equilibrium Extrapolation
		// shadowCell.MsQ().V = -cell->MsQ().V;
		// //
		// Update_phi_Eq(shadowCell);
	#ifndef _CARTESIAN_MESH_FLIP
		for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
		#ifdef _ARK_ALLENCAHN_FLIP
			shadowCell.h.BarP[i][j] = shadowCell.h.Eq[i][j]
			+ cell->h.aBP*(cell->h.Tilde[i][j] - cell->h.Eq[i][j]);
		#endif
			shadowCell.f.BarP[i][j] = shadowCell.f.Eq[i][j]
			+ cell->f.aBP*(cell->f.Tilde[i][j] - cell->f.Eq[i][j]);
	//isothermal flip
		#ifndef _ARK_ISOTHERMAL_FLIP	
			shadowCell.g.BarP[i][j] = shadowCell.g.Eq[i][j]
			+ cell->g.aBP*(cell->g.Tilde[i][j] - cell->g.Eq[i][j]);
		#endif
		}
	#endif
}
void Wall_3_NEE(Face_2D &face)
{
	Cell_2D &cell = *face.owner;
	face.MsQh().Rho = cell.MsQ().Rho;
	face.MsQh().Phi = cell.MsQ().Phi;
	face.MsQh().p = cell.MsQ().p;
	// face.MsQh().U = 0.0;
	// face.MsQh().V = 0.0;
	Update_phi_Eqh(face);
	//
	#ifdef _ARK_ALLENCAHN_FLIP
	double hNEq = 2.0*cell.h.tau/(2.0*cell.h.tau + ::dt);
	#endif
	//
	double fNEq = 2.0*cell.f.tau/(2.0*cell.f.tau + ::dt);
	//
	#ifndef _ARK_ISOTHERMAL_FLIP	
	double gNEq = 2.0*cell.g.tau/(2.0*cell.g.tau + ::dt);
	#endif
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
	#ifdef _ARK_ALLENCAHN_FLIP
		face.h.hDt[i][j] = face.h.EqhDt[i][j]
		+ hNEq*(cell.h.Tilde[i][j] - cell.h.Eq[i][j] + hDt*cell.h.So[i][j]);
		//
		face.h.hDt[i][j] *= face.xi_n_dS[i][j];
	#endif
		face.f.hDt[i][j] = face.f.EqhDt[i][j]
		+ fNEq*(cell.f.Tilde[i][j] - cell.f.Eq[i][j] + hDt*cell.f.So[i][j]);
		//
		face.f.hDt[i][j] *= face.xi_n_dS[i][j];
	//isothermal flip
	#ifndef _ARK_ISOTHERMAL_FLIP	
		face.g.hDt[i][j] = face.g.EqhDt[i][j]
		+ gNEq*(cell.g.Tilde[i][j] - cell.gEq[i][j] + hDt*cell.g.So[i][j]);
		//
		face.g.hDt[i][j] *= face.xi_n_dS[i][j];
	#endif
	}
}
