#include "DUGKSDeclaration.h"
//
using PhaseFieldAC::PhiL;
using PhaseFieldAC::PhiV;
using PhaseFieldAC::RhoL;
using PhaseFieldAC::RhoV;
using PhaseFieldAC::wI;
using PhaseFieldAC::Gx;
using PhaseFieldAC::Gy;

extern void Grad_Phi_6points(Cell_2D *cellptr);

extern void Grad_Phi_CD(Cell_2D *center);

inline double SourcePhi(double phi)
{
	return -4*(phi-PhiL)*(phi-PhiV)/(wI*(PhiL-PhiV));
}
inline double ChemicalPotential(double Phi,double laplacianPhi)
{
	return
	(
	2*PhaseFieldAC::Beta*(Phi-PhiL)*(Phi-PhiV)*(2*Phi-PhiL-PhiV)
	-PhaseFieldAC::Kappa*laplacianPhi
	);
}
void SourceMomentum(Cell_2D *cellptr)
{
	double F = ChemicalPotential(cellptr->MsQ().Phi,cellptr->MsQ().laplacianPhi);
	cellptr->msq->calcRho_xRho_y();
	cellptr->msq->calcFxFy(F);
//	cellptr->msq->Fy -= Gy*(cellptr->msq->Rho-RhoL);
	// cellptr->msq->Fx += Gx;
}
void MacroSource(Cell_2D *cellptr)
{
//
	Grad_Phi_6points(cellptr);
//	Grad_Phi_CD(cellptr);
//
	#ifdef _ARK_MOMENTUM_FLIP
	#ifdef _ARK_FORCE_FLIP
	SourceMomentum(cellptr);
	#endif
	#endif
//
	double L = sqrt(cellptr->MsQ().SqPhixPhiy());
	if(L!= 0)
	{
		double tmp = SourcePhi(cellptr->MsQ().Phi);
		(cellptr->MsQ().Phi_x) = (cellptr->MsQ().Phi_x)*tmp/L + cellptr->MsQ().dPhiU()/(RT*dt);
		(cellptr->MsQ().Phi_y) = (cellptr->MsQ().Phi_y)*tmp/L + cellptr->MsQ().dPhiV()/(RT*dt);
	}
	//
	cellptr->MsQ().prevPhiU = (cellptr->MsQ().Phi)*(cellptr->MsQ().U);
	cellptr->MsQ().prevPhiV = (cellptr->MsQ().Phi)*(cellptr->MsQ().V);
}
