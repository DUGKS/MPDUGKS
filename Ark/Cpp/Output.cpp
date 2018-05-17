#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <string>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::map;
using std::ios;
using std::setprecision;

int const Out_precision = 12;

int const IC = 198, JC = 197;

extern void TaylorGreenVortex(double t,double x,double y,double &u, double &v, double &p);

extern void TaylorCouetteAnalytical(double x,double y,double &u_A);

extern void AnalyticalForceDrivenTG(double x,double y,double &u_A, double &v_A,double &p_A);

extern void LayeredPoiseuilleAnalytical(double const xc,double const yc,double &u_A);

void OutputCase()
{
	ofstream OutFile_Case("../FlowField/Case.ark");
	if(!OutFile_Case)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"Case.ark open failed."<<endl;
		getchar();
		return;
	}
	OutFile_Case <<_MESHFILE_NAME_ARK<<"    =    Mesh File"<<endl
				 <<_MESHTYPE_ARK<<"    =    Mesh Type"<<endl
				 <<_FLUX_SCHEME_ARK<<"    =    Flux Scheme"<<endl
				 <<_BC_ARK<<"    =    Boundary Condition"<<endl
				 <<_QMODEL_ARK<<"    =    velocity model"<<endl
				 <<_FORCE_MODEL_ARK<<"    =    Force model"<<endl
				 <<_EOS_MODEL_ARK<<"    =    EoS model"<<endl
				 <<left<<fs<<"="<<fs<<"left"<<endl
				 <<right<<fs<<"="<<fs<<"right"<<endl
				 <<top<<fs<<"="<<fs<<"top"<<endl
				 <<bottom<<fs<<"="<<fs<<"bottom"<<endl
				 <<"#----------------------------------------"<<endl
				 <<Omega0<<"    =    Omega0//sutherland power"<<endl
				 <<Pr<<"    =    Prandtl"<<endl
				 <<nK<<"    =    nK//internal degrees of freedom"<<endl
				 <<Cv<<"    =    specific heat of constant volume"<<endl
				 <<Gamma<<"    =    Gamma//specific heat ratio"<<endl
				 <<"#----------------------------------------"<<endl
				 <<DV_Qu<<"    =    Qu"<<endl
				 <<DV_Qv<<"    =    Qv"<<endl
				 <<NL<<"    =    NL"<<endl
				 <<ChLength<<"    =    Character Length"<<endl
				 <<MinL<<"    =    MinL"<<endl
				 <<Kn<<"    =    Knudsen"<<endl
				 <<Ma<<"    =    Mach"<<endl
				 <<Re<<"    =    Re"<<endl
				 <<U0<<"    =    U0"<<endl
				 <<R0<<"    =    R0"<<endl
				 <<T0<<"    =    T0"<<endl
				 <<Mu0<<"    =    dynamic viscosity"<<endl
				 <<Nu0<<"    =    kinetic viscosity"<<endl
				 <<Tau0<<"    =    tau0"<<endl
				 <<Lambda0<<"    =    Lambda0"<<endl
				 <<CFL<<"    =    CFL"<<endl
				 <<dt<<"    =    dt"<<endl				 
				 <<dt/Tau0<<"    =    dt/tau0"<<endl
				 <<RESIDUAL<<"    =    RESIDUAL"<<endl
				 <<MaSpan<<"    =    MaSpan"<<endl
				 <<Eta<<"    =    Eta// T of D2V16"<<endl
				 <<"-------------Phase Field--------------"<<endl
				 <<PhaseFieldAC::Cn<<"    =    Cahn number"<<endl
				 <<PhaseFieldAC::Pe<<"    =    Peclet number"<<endl
				 <<PhaseFieldAC::Eo<<"    =    Eotvos(Bond) number"<<endl
				 <<PhaseFieldAC::Mo<<"    =    Morton number"<<endl
				 <<PhaseFieldAC::ReMP<<"    =    Reynolds number"<<endl
				 <<endl
				 <<PhaseFieldAC::M_Phi<<"    =    mobility"<<endl
				 <<PhaseFieldAC::W_InterFace<<"    =    inteface width"<<endl
				 <<PhaseFieldAC::Sigma<<"    =    Sigma(surface tension coefficient)"<<endl
				 <<PhaseFieldAC::Beta<<"    =    Beta"<<endl
				 <<PhaseFieldAC::Kappa<<"    =    Kappa"<<endl
				 <<PhaseFieldAC::RhoL<<"    =    liquid density"<<endl
				 <<PhaseFieldAC::RhoV<<"    =    vapor density"<<endl
				 <<PhaseFieldAC::NuL<<"    =    liquid viscosity"<<endl
				 <<PhaseFieldAC::NuV<<"    =    vapor viscosity"<<endl
				 <<PhaseFieldAC::TauMass<<"    =    mobility relaxition time";

	OutFile_Case.close();
//
}
void Output_Convergence()
{
	ostringstream oss_L2;
	oss_L2 <<"../Convergence/L2_uvp_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat"; 
	ofstream OutFile_L2(oss_L2.str().c_str());
	if(!OutFile_L2)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_L2 file open failed"<<endl;
		getchar();
	}
//---------------------------------------------------------------------
	ostringstream oss_SumRho;
	oss_SumRho << "../Convergence/SumRho_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumRho(oss_SumRho.str().c_str());
	if(!OutFile_SumRho)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumRho file open failed"<<endl;
		getchar();
	}
//---------------------------------------------------------------------
	ostringstream oss_Residual;
	oss_Residual <<"../Convergence/Residual_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_Residual(oss_Residual.str().c_str());
	if(!OutFile_Residual)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_Residual open failed"<<endl;
		getchar();
	}
}
void TaylorGreen_L2Norm(double const &t,double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp; 
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].MsQ().U - u_A;
		dv = CellArray[i].MsQ().V - v_A;
		dp = CellArray[i].MsQ().p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void TaylorCouette_L2Norm(double const &t,double &L2_uv)
{
	double u_A, du;
	double  Sumdudv = 0.0,Sumuv_A = 0.0;
	for(int i =0;i < Cells;++i)
	{
		TaylorCouetteAnalytical(CellArray[i].xc,CellArray[i].yc,u_A);
		du = sqrt(CellArray[i].MsQ().SqUV()) - u_A;
		Sumdudv += du*du;
		Sumuv_A += u_A*u_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
}
void LayeredPoiseuille_L2Norm(double const &t,double &L2_uv)
{
	double u_A, du;
	double  Sumdudv = 0.0,Sumuv_A = 0.0;
	for(int i =0;i < Cells;++i)
	{
		LayeredPoiseuilleAnalytical(CellArray[i].xc,CellArray[i].yc,u_A);
		du = CellArray[i].MsQ().U - u_A;
		Sumdudv += du*du;
		Sumuv_A += u_A*u_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
}
void ForceDrivenTaylorGreen_L2Norm(double &L2_uv, double &L2_p)
{
	double u_A, v_A, p_A, du, dv, dp;
	double  Sumdudv = 0.0,Sumdp = 0.0,Sumuv_A = 0.0,Sump_A = 0.0;
	for(int i = 0;i < Cells;++i)
	{
		AnalyticalForceDrivenTG(CellArray[i].xc,CellArray[i].yc,u_A,v_A,p_A);
		du = CellArray[i].MsQ().U - u_A;
		dv = CellArray[i].MsQ().V - v_A;
		dp = CellArray[i].MsQ().p - p_A;
		Sumdudv += du*du + dv*dv;
		Sumuv_A += u_A*u_A + v_A*v_A; 
		Sumdp   += dp*dp;
		Sump_A  += p_A*p_A;
	}
	L2_uv = sqrt(Sumdudv/Sumuv_A);
	L2_p  = sqrt(Sumdp/Sump_A);
}
void Output_L2Norm(double const &t,double &L2_uv, double &L2_p)
{	
	//ForceDrivenTaylorGreen_L2Norm(L2_uv,L2_p);
	//TaylorCouette_L2Norm(t,L2_uv);
	LayeredPoiseuille_L2Norm(t,L2_uv);
	ostringstream oss_L2;
	oss_L2 <<"../Convergence/L2_uvp_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_L2(oss_L2.str().c_str(),ofstream::app);
	if(!OutFile_L2)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_L2 file open failed"<<endl;
		getchar();
	}
	OutFile_L2 << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_L2 << t <<"    "<<L2_uv<<"    "<<L2_p<<'\n';
	OutFile_L2.close();
}
void Output_SumRho(double t)
{
	ostringstream oss_SumRho;
	oss_SumRho << "../Convergence/SumRho_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_SumRho(oss_SumRho.str().c_str(),ofstream::app);
	if(!OutFile_SumRho)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_SumRho file open failed"<<endl;
		getchar();
	}
	OutFile_SumRho << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_SumRho <<  t <<"    "<<SumRho<<endl;
	OutFile_SumRho.close();
}
void Output_Residual(double t,double Residual)
{
	ostringstream oss_Residual;
	oss_Residual <<"../Convergence/Residual_mu"<<Mu0<<"_Re"<<Re<<"_"<<NL<<_MESHFILE_NAME_ARK<<".dat";
	ofstream OutFile_Residual(oss_Residual.str().c_str(),ofstream::app);
	if(!OutFile_Residual)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_Residual open failed"<<endl;
		getchar();
	}
	OutFile_Residual << setiosflags(ios::scientific)<<setprecision(8);
	OutFile_Residual << t <<"    "<<Residual<<endl;
	OutFile_Residual.close();
}
void Output_Flowfield(double const &t,int step)
{
	ostringstream oss_FlowField;
	oss_FlowField <<"../FlowField/global/" << "step" << step <<"Ma"<< Ma<<"_"
					<<_MESHTYPE_ARK<<NL<<"_CFL"<<CFL<<"_T"<<t<<".dat";
	ofstream OutFile_FlowField(oss_FlowField.str().c_str());
	if(!OutFile_FlowField)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_FlowField open failed" << endl; 
		getchar();
		return;
	}
//	q<sub>x</sub>,q<sub>y</sub>,\
//	<Greek>t</Greek><sub>xx</sub>,<Greek>t</Greek><sub>xy</sub>,<Greek>t</Greek><sub>yy</sub>,
	ostringstream VarName,VarLocation,ZoneName,dataNE;
	VarName << "VARIABLES = X,Y,<Greek>r</Greek>,U,V,p,T,\
	<Greek>f</Greek>,\
	<Greek>f</Greek><sub>x</sub>,<Greek>f</Greek><sub>y</sub>\n";
	VarLocation <<"VarLocation=([1-2]=NODAL,[3-10]=CellCentered)\n";
	ZoneName<<"ZONE T = Time" << t <<"_"<<"Mu"<<Mu0<<"\n";
	dataNE<<"Nodes="<<Nodes<<", Elements="<<Cells<<", ZONETYPE=FEQuadrilateral\n";
	string tecformat[5]={VarName.str().c_str(),
						ZoneName.str().c_str(),
						dataNE.str().c_str(),
						"DATAPACKING=BLOCK\n",
						VarLocation.str().c_str()};
	OutFile_FlowField << tecformat[0]<<tecformat[1]<<tecformat[2]<<tecformat[3]<<tecformat[4];
	OutFile_FlowField << setiosflags(ios::scientific) << setprecision(12);
//	
	/*double *u_A = new double[Cells];
	double *v_A = new double[Cells];
	double *p_A = new double[Cells];
	for(int i = 0;i != Cells;++i)
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
	*/
//---------------------------------------------Node->X-------------------------------------------
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_FlowField << NodeArray[i].xN <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	OutFile_FlowField << endl;
//--------------------------------------------Node->Y-------------------------------------------
	for(int i = 0;i != Nodes;++i)
	{
		OutFile_FlowField << NodeArray[i].yN <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	OutFile_FlowField << endl;
//--------------------------------------------rho-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().Rho <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------u-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().U <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------v-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().V <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------p-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().p <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------T-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().T <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
/*
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().qx <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().qy <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[0][0] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[0][1] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].shearTau[1][1] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
*/
	// for(int i = 0;i != Cells;++i)
	// {
	// 	OutFile_FlowField << (*CellArray[i].pseudopsi) <<"   ";
	// 	if((i+1)%16 == 0)
	// 		OutFile_FlowField << "\n";
	// }
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().Phi <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().Phi_x <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << CellArray[i].MsQ().Phi_y <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
/*//--------------------------------------------u_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << u_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------v_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << v_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}
//--------------------------------------------p_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << p_A[i] <<"   ";
		if((i+1)%16 == 0)
			OutFile_FlowField << "\n";
	}*/
//--------------------------------------------relation----------------------------------------	
	for(int i = 0;i != Cells;++i)
	{
		OutFile_FlowField << MeshIndex(CellArray[i].cellNodes[0] , NodeArray)<< " "
						  << MeshIndex(CellArray[i].cellNodes[1] , NodeArray)<< " "
					 	  << MeshIndex(CellArray[i].cellNodes[2] , NodeArray)<< " " 
						  << MeshIndex(CellArray[i].cellNodes[3] , NodeArray)<< endl;
	}
	OutFile_FlowField.close();
	/*delete []u_A;
	delete []v_A;
	delete []p_A;*/
}
void Output_xcyc()
{
	ostringstream oss_xcyc;
	oss_xcyc <<"../FlowField/step" <<step<< "Nu"<< Nu0
					<<"_MeshCar"<<NL<<"_CFL"<<CFL<<"_xcyc.dat";
	ofstream OutFile_xcyc(oss_xcyc.str().c_str());
	if(!OutFile_xcyc)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"OutFile_xcyc open failed" << endl; 
		getchar();
		return;
	}
	ostringstream title,varName,dataIJK,VarLocation;
	title << "title = \"Step"<<step<<"_"<<"Nu"<<Nu0<<"\""<<"\n";
	varName << "variables = x,y,<Greek>r</Greek>,u,v,p,Fx,Fy\n";
	dataIJK << "zone I = "<<Nx<<", J = "<<Ny<<", F = BLOCK\n";
	string tec360Header[3]={title.str().c_str(),varName.str().c_str(),
							dataIJK.str().c_str()};
//
	OutFile_xcyc<<tec360Header[0]<<tec360Header[1]<<tec360Header[2];
	OutFile_xcyc << setiosflags(ios::scientific) << setprecision(Out_precision);
	// for(int i = 0;i != Cells;++i)
	// {
	// 	OutFile_xcyc << CellArray[i].xc<<"  "<<CellArray[i].yc<<"\n";
	// }
	// OutFile_xcyc <<"--------------P_Inlet-----------------"<<'\n';
	// for(int n = 0;n != P_InletFaceNum;++n)
	// 	OutFile_xcyc << P_InletShadowCA[n].xc<<"  "<<P_InletShadowCA[n].yc<<"\n";
	// OutFile_xcyc <<"--------------P_Outlet----------------"<<'\n';
	// for(int n = 0;n != P_OutletFaceNum;++n)
	// 	OutFile_xcyc << P_OutletShadowCA[n].xc<<"  "<<P_OutletShadowCA[n].yc<<"\n";
	// for(int n = 0;n != PeriodicFaceNum;++n)
	// 	OutFile_xcyc << PeriodicShadowCA[n].xc<<"  "<<PeriodicShadowCA[n].yc<<"\n";
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->xc<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//	
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->yc<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->Rho<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->U<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->V<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << CarCellArray[i][j]->msq->p<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	// for(int i = 1;i < Nxp1;++i)
	// for(int j = 1;j < Nyp1;++j)
	// {
	// 	OutFile_xcyc << *(CarCellArray[i][j]->pseudopsi)<<"  ";
	// 	if(j%10 == 0)
	// 		OutFile_xcyc<<"\n";
	// }
	// OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << (CarCellArray[i][j]->msq->Fx)<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		OutFile_xcyc << (CarCellArray[i][j]->msq->Fy)<<"  ";
		if(j%10 == 0)
			OutFile_xcyc<<"\n";
	}
	OutFile_xcyc<<"\n";
//
	OutFile_xcyc.close();
}
//----------------------------------------------DEBUG----------------------------------------
#ifdef _ZERO_NDEBUG_FLIP
void oss_XXX(ostringstream& oss_r,const string &folder,const string &suffix,double const &t)
{
	oss_r << "../FlowField/"<<folder<<"/Time" <<t<<"."<<suffix;
}
void FileOpen(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str());
	if(!OutFile_XXX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void FileOpenAppend(ofstream &OutFile_XXX,ostringstream &oss_XXX,string const &s)
{
	OutFile_XXX.open(oss_XXX.str().c_str(),ios::app);
	if(!OutFile_XXX)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<< ("OutFile_" + s + " open failed") << endl; 
		getchar();
		return;
	}
	OutFile_XXX << setiosflags(ios::scientific) << setprecision(Out_precision);
}
void printCorners(Cell_2D *cellptr)
{
	cout <<*cellptr->use<<"    "<<cellptr->xc<<"    "<<cellptr->yc<<endl;
}
void Output_UVP(double const &t)
{
	ostringstream oss_rho,oss_u,oss_v,oss_p,oss_uA,oss_vA,oss_pA;
	ofstream OutFile_rho,OutFile_u, OutFile_v, OutFile_p, OutFile_uA, OutFile_vA, OutFile_pA;
//--------------------------------------------------------------------------------
	double *u_A = new double[Cells];
	double *v_A = new double[Cells];
	double *p_A = new double[Cells];
	for(int i = 0;i != Cells;++i)
		TaylorGreenVortex(t,CellArray[i].xc,CellArray[i].yc,u_A[i], v_A[i], p_A[i]);
//---------------------------------------------------------------------------------
	oss_XXX(oss_rho,"UVP","rho",t);
	oss_XXX(oss_u,"UVP","u",t);
	oss_XXX(oss_v,"UVP","v",t);
	oss_XXX(oss_p,"UVP","p",t);
	oss_XXX(oss_uA,"UVP","uA",t);
	oss_XXX(oss_vA,"UVP","vA",t);
	oss_XXX(oss_pA,"UVP","pA",t);
	FileOpen(OutFile_rho,oss_rho,"rho");
	FileOpen(OutFile_u,oss_u,"u");
	FileOpen(OutFile_v,oss_v,"v");
	FileOpen(OutFile_p,oss_p,"p");
	FileOpen(OutFile_uA,oss_uA,"uA");
	FileOpen(OutFile_vA,oss_vA,"vA");
	FileOpen(OutFile_pA,oss_pA,"pA");
//--------------------------------------------rho-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_rho << CellArray[i].MsQ().Rho <<"\n";
	}
//--------------------------------------------u-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_u << CellArray[i].MsQ().U <<"\n";
	}
//--------------------------------------------v-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_v << CellArray[i].MsQ().V <<"\n";
	}
//--------------------------------------------p-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_p << CellArray[i].MsQ().p <<"\n";
	}
//--------------------------------------------u_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_uA << u_A[i] <<"\n";
	}
//--------------------------------------------v_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_vA << v_A[i] <<"\n";
	}
//--------------------------------------------p_A-------------------------------------------
	for(int i = 0;i != Cells;++i)
	{
		OutFile_pA << p_A[i] <<"\n";
	}
	OutFile_rho.close();
	OutFile_u.close();
	OutFile_v.close();
	OutFile_p.close();
	OutFile_uA.close();
	OutFile_vA.close();
	OutFile_pA.close();
	delete []u_A;
	delete []v_A;
	delete []p_A;
}
void Output_fBP(double const &t,int ii,int jj)
{
	ostringstream oss_fBP;
	oss_fBP <<"../FlowField/fBP/" << "Time" << t <<"_Mu"<< Mu0
					<<"_MeshCar"<<NL<<"-"<<NL<<"_CFL"<<CFL<<"_fBP"<<ii<<"_"<<jj<<".dat";
	ofstream OutFile_fBP(oss_fBP.str().c_str());
	if(!OutFile_fBP)
	{
		_PRINT_ERROR_MSG_FLIP
		cout <<"  "<<"OutFile_fBP open failed" << endl;
		getchar();
		return;
	}
	OutFile_fBP << setiosflags(ios::scientific) << setprecision(12);
	for(int n = 0;n != Cells;++n)
	{
		OutFile_fBP << CellArray[n].f.BarP[ii][jj] <<"\n";
	}
	OutFile_fBP.close();
}
void Output_fBh(Face_2D& face,double t)
{
	ostringstream oss_fBh;
	ofstream OutFile_fBh;
	oss_XXX(oss_fBh,"fBh","fBh",t);
	FileOpen(OutFile_fBh,oss_fBh,"fBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fBh<<face.f.BhDt[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_fBh <<face.owner->f.BarP[i][j]<<"    "
		// 				<<face.owner->f.Eq[i][j]<<"    "
		// 				<<face.f.BhDt[i][j] - face.owner->f.BarP[i][j]<<"    "
		// 				<<face.owner->f.Eq[i][j] - face.owner->f.BarP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_fBh <<face.neigh->f.BarP[i][j]<<"    "
		// 				<<face.neigh->f.Eq[i][j]<<"    "
		// 				<<face.f.BhDt[i][j] - face.neigh->f.BarP[i][j]<<"    "
		// 				<<face.neigh->f.Eq[i][j] - face.neigh->f.BarP[i][j]<<endl;
		// }
	}
	OutFile_fBh.close();
}
void Output_fh(Face_2D& face,double t)
{
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fh",t);
	FileOpen(OutFile_fh,oss_fh,"fh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_fh << face.f.hDt[i][j] <<"    "
					   << face.f.BhDt[i][j]<<"    "
					   << face.f.EqhDt[i][j]<<'\n';
		}
	OutFile_fh <<"Rho_h : "<<face.MsQh().Rho<<'\n'
			   <<"U_h : "<<face.MsQh().U<<'\n'
			   <<"V_h : "<<face.MsQh().V<<'\n'
			   <<"T_h: "<<face.MsQh().T<<'\n'
			   <<"qx_h : "<<face.MsQh().qx<<'\n'
			   <<"qy_h : "<<face.MsQh().qy<<'\n'
			   <<"xc : "<<face.xf<<'\n'
			   <<"yc : "<<face.yf<<'\n';
	OutFile_fh.close();
}
void Output_fh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fh;
	ofstream OutFile_fh;
	oss_XXX(oss_fh,"fh","fhAPP",dt);
	FileOpenAppend(OutFile_fh,oss_fh,"fh");
	OutFile_fh 	<< face.f.hDt[I][J] <<"    "<<face.f.ah<<"    "<<face.f.bh<<'\n';
				//<< face.f.BhDt[i][j]<<"    "
				//<< face.f.EqhDt[i][j]<<'\n';
	OutFile_fh.close();				
}
void Output_fT(Cell_2D &cell,double t)
{
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fT",t);
	FileOpen(OutFile_fT,oss_fT,"fT");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_fT <<cell.f.Tilde[i][j]<<"    "
				   <<cell.f.Eq[i][j]<<"    "
//				   <<cell.fFlux[i][j]<<"    "
				   <<cell.f.Tilde[i][j]-cell.f.Eq[i][j]
				   <<'\n';
	}
	OutFile_fT <<"Rho : "<<cell.MsQ().Rho<<'\n'
				<<"U : "<<cell.MsQ().U<<'\n'
				<<"V : "<<cell.MsQ().V<<'\n'
				<<"T: "<<cell.MsQ().T<<'\n'
				<<"Lambda: "<<cell.MsQ().Lambda<<'\n'
				<<"qx : "<<cell.MsQ().qx<<'\n'
				<<"qy : "<<cell.MsQ().qy<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_fT.close(); 
}
void Output_fT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_fT;
	ofstream OutFile_fT;
	oss_XXX(oss_fT,"fT","fTAPP",dt);
	FileOpenAppend(OutFile_fT,oss_fT,"fT");
	OutFile_fT  <<cell.f.Tilde[I][J]<<"    "//<<cell.aBP<<"    "<<cell.bBP
				<<cell.Cell_C[0]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[1]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[2]->f.Tilde[I][J]<<"    "
				<<cell.Cell_C[3]->f.Tilde[I][J]<<"    "
				<<'\n';
	OutFile_fT.close();
}

// void Output_fFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_fFlux;
// 	ofstream OutFile_fFlux;
// 	oss_XXX(oss_fFlux,"fFlux","fFluxAPP",dt);
// 	FileOpenAppend(OutFile_fFlux,oss_fFlux,"fFlux");
// 	OutFile_fFlux  <<cell.fFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->fh[I][J]<<"    "
// 				    <<cell.Face_C[1]->fh[I][J]<<"    "
// 				    <<cell.Face_C[2]->fh[I][J]<<"    "
// 				    <<cell.Face_C[3]->fh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_fFlux.close();
// }
// void Output_gFlux_Append(Cell_2D &cell,double dt)
// {
// 	int I = IC,J = JC;
// 	ostringstream oss_gFlux;
// 	ofstream OutFile_gFlux;
// 	oss_XXX(oss_gFlux,"gFlux","gFluxAPP",dt);
// 	FileOpenAppend(OutFile_gFlux,oss_gFlux,"gFlux");
// 	OutFile_gFlux  <<cell.gFlux[I][J]<<"    "
// 				    <<cell.Face_C[0]->gh[I][J]<<"    "
// 				    <<cell.Face_C[1]->gh[I][J]<<"    "
// 				    <<cell.Face_C[2]->gh[I][J]<<"    "
// 				    <<cell.Face_C[3]->gh[I][J]<<"    "
// 				    <<'\n';
// 	OutFile_gFlux.close();
// }
#ifndef _ARK_ISOTHERMAL_FLIP
void Output_gBh(Face_2D& face,double t)
{
	ostringstream oss_gBh;
	ofstream OutFile_gBh;
	oss_XXX(oss_gBh,"gBh","gBh",t);
	FileOpen(OutFile_gBh,oss_gBh,"gBh");
	//
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		OutFile_gBh<<face.g.BhDt[i][j]<<'\n';
		// if(face.xi_n_dS[i][j] >= 0)
		// {
		// 	OutFile_gBh <<face.owner->g.BarP[i][j]<<"    "
		// 				<<face.owner->gEq[i][j]<<"    "
		// 				<<face.g.BhDt[i][j] - face.owner->g.BarP[i][j]<<"    "
		// 				<<face.owner->gEq[i][j] - face.owner->g.BarP[i][j]<<endl;
		// }
		// else
		// {
		// 	OutFile_gBh <<face.neigh->g.BarP[i][j]<<"    "
		// 				<<face.neigh->gEq[i][j]<<"    "
		// 				<<face.g.BhDt[i][j] - face.neigh->g.BarP[i][j]<<"    "
		// 				<<face.neigh->gEq[i][j] - face.neigh->g.BarP[i][j]<<endl;
		// }
	}
	OutFile_gBh.close();
}
void Output_gT(Cell_2D &cell,double t)
{
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gT",t);
	FileOpen(OutFile_gT,oss_gT,"gT");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gT <<cell.g.Tilde[i][j]<<"    "
//					   <<cell.f.Eq[i][j]<<"    "
//					   <<cell.fFlux[i][j]<<"    "
//					   <<cell.g.Tilde[i][j]-cell.f.Eq[i][j]-cell.DtSlashVolume*cell.fFlux[i][j]
					   <<endl;
		}
	OutFile_gT <<"Rho : "<<cell.Rho<<'\n'
				<<"U : "<<cell.U<<'\n'
				<<"V : "<<cell.V<<'\n'
				<<"T: "<<cell.T<<'\n'
				<<"Lambda: "<<cell.Lambda<<'\n'
				<<"qx : "<<cell.qx<<'\n'
				<<"qy : "<<cell.qy<<'\n'
				<<"xc : "<<cell.xc<<'\n'
				<<"yc : "<<cell.yc<<'\n';
	OutFile_gT.close(); 
}
void Output_gh_Append(Face_2D& face,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","ghAPP",dt);
	FileOpenAppend(OutFile_gh,oss_gh,"gh");
	OutFile_gh 	<< face.g.hdt[I][J] <<"    "<<'\n';
				//<< face.f.BhDt[i][j]<<"    "
				//<< face.f.EqhDt[i][j]<<'\n';
	OutFile_gh.close();				
}
void Output_gh(Face_2D& face,double t)
{
	ostringstream oss_gh;
	ofstream OutFile_gh;
	oss_XXX(oss_gh,"gh","gh",t);
	FileOpen(OutFile_gh,oss_gh,"gh");
	for(int i = 0;i < DV_Qu;++i)
		for(int j = 0;j < DV_Qv;++j)
		{
			OutFile_gh << face.g.hDt[i][j] <<"    "
					   << face.g.BhDt[i][j]<<"    "
					   << face.g.EqhDt[i][j]<<'\n';
		}
	OutFile_gh.close();
}
void Output_gT_Append(Cell_2D &cell,double dt)
{
	int I = IC,J = JC;
	ostringstream oss_gT;
	ofstream OutFile_gT;
	oss_XXX(oss_gT,"gT","gTAPP",dt);
	FileOpenAppend(OutFile_gT,oss_gT,"gT");
	OutFile_gT  <<cell.g.Tilde[i][j]<<"    "
				// <<cell.g.BarP[i][j]<<"    "
				// <<cell.g.Eq[i][J]<<"    "
				<<cell.Cell_C[0]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[1]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[2]->g.Tilde[i][j]<<"    "
				<<cell.Cell_C[3]->g.Tilde[i][j]
				<<'\n';
	OutFile_gT.close();
}
#endif
void Output_phi_Bh(Face_2D &face,double t)
{
	ostringstream oss_fBh_L;
	ostringstream oss_fBh_R;
	ostringstream oss_gBh_L;
	ostringstream oss_gBh_R;
	ofstream OutFile_fBh_L;
	ofstream OutFile_fBh_R;
	ofstream OutFile_gBh_L;
	ofstream OutFile_gBh_R;
	oss_XXX(oss_fBh_L,"fBh","fBhL",t);
	oss_XXX(oss_fBh_R,"fBh","fBhR",t);
	oss_XXX(oss_gBh_L,"fBh","gBhL",t);
	oss_XXX(oss_gBh_R,"fBh","gBhR",t);
	FileOpen(OutFile_fBh_L,oss_fBh_L,"fBh_L");
	FileOpen(OutFile_fBh_R,oss_fBh_R,"fBh_R");
	FileOpen(OutFile_gBh_L,oss_gBh_L,"gBh_L");
	FileOpen(OutFile_gBh_R,oss_gBh_R,"gBh_R");
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		if(xi_u[QuIndex]>=0)
		{
			OutFile_fBh_R<<face.f.BhDt[i][j]<<'\n';
			#ifndef _ARK_ISOTHERMAL_FLIP
			OutFile_gBh_R<<face.g.BhDt[i][j]<<'\n';
			#endif
		}
		else
		{
			OutFile_fBh_L<<face.f.BhDt[i][j]<<'\n';
			#ifndef _ARK_ISOTHERMAL_FLIP
			OutFile_gBh_L<<face.g.BhDt[i][j]<<'\n';
			#endif
		}
	}
	OutFile_fBh_L.close();
	OutFile_fBh_R.close();
	OutFile_gBh_L.close();
	OutFile_gBh_R.close();
}
#endif