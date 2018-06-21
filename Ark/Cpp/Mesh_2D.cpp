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
Cell_2D::Cell_2D():msq(new MacroQuantity()),use(new int(1))
{
}
Cell_2D::Cell_2D(const Cell_2D &rhs)
{
//
	#ifdef _ARK_ALLENCAHN_FLIP
	h = rhs.h;
	#endif
	//
	#ifdef _ARK_MOMENTUM_FLIP
	f = rhs.f;
	#endif
	//
	#ifndef _ARK_ISOTHERMAL_FLIP
	g = rhs.g;
	#endif
//
	msq = rhs.msq;
//
	use   = rhs.use;
	++*use;
}
Cell_2D& Cell_2D::operator=(const Cell_2D &rhs)
{
    if(use != rhs.use)
    {
		if(--*use == 0)
		{	
			delete msq;
			delete use;
		}
		#ifdef _ARK_ALLENCAHN_FLIP
		h = rhs.h;
		#endif
		//!momentum
		#ifdef _ARK_MOMENTUM_FLIP
		f = rhs.f;
		#endif
		//
		#ifndef _ARK_ISOTHERMAL_FLIP
		g = rhs.g;
		#endif
		msq = rhs.msq;
		//
		use   = rhs.use;
		++*use;
    }
//
	return *this;
}
Cell_2D::~Cell_2D()
{
	if(--*use == 0)
	{
		delete msq;
		delete use;
	}
}
Cell_2D::DVDF::DVDF()
{
	AllocateARK(BarP,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(BarP_x,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(BarP_y,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(Tilde,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(Eq,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(So,Qu,Qv);
}
Cell_2D::DVDF::DVDF(const DVDF& rhs)
{
	BarP = rhs.BarP;
	BarP_x = rhs.BarP_x;
	BarP_y = rhs.BarP_y;
	Tilde  = rhs.Tilde;
	Eq  = rhs.Eq;
	So  = rhs.So;
//
	token = rhs.token;
	++*token;
}
Cell_2D::DVDF& Cell_2D::DVDF::operator=(const DVDF& rhs)
{
	if(token != rhs.token)
	{
		if(--*token == 0)
		{
			DeallocateARK(BarP,Cell_2D::Qu,Cell_2D::Qv);
			DeallocateARK(BarP_x,Cell_2D::Qu,Cell_2D::Qv);
			DeallocateARK(BarP_y,Cell_2D::Qu,Cell_2D::Qv);
			DeallocateARK(Tilde,Cell_2D::Qu,Cell_2D::Qv);
			DeallocateARK(Eq,Cell_2D::Qu,Cell_2D::Qv);
			DeallocateARK(So,Qu,Qv);
		//
			delete token;
		}
		BarP = rhs.BarP;
		BarP_x = rhs.BarP_x;
		BarP_y = rhs.BarP_y;
		Tilde  = rhs.Tilde;
		Eq  = rhs.Eq;
		So  = rhs.So;
//
		token = rhs.token;
		++*token;
	}
	return *this;
}
Cell_2D::DVDF::~DVDF()
{
	if(--*token == 0)
	{
		DeallocateARK(BarP,Cell_2D::Qu,Cell_2D::Qv);
		DeallocateARK(BarP_x,Cell_2D::Qu,Cell_2D::Qv);
		DeallocateARK(BarP_y,Cell_2D::Qu,Cell_2D::Qv);
		DeallocateARK(Tilde,Cell_2D::Qu,Cell_2D::Qv);
		DeallocateARK(Eq,Cell_2D::Qu,Cell_2D::Qv);
		DeallocateARK(So,Qu,Qv);
	}
}
void Cell_2D::DVDF::setxBP()
{
	aBP = (2.0*tau-::hDt)/(2.0*tau + ::dt);
	bBP = 1.0 - aBP;
	cBP = tau*bBP;
}
void Cell_2D::Factor()
{
//
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	f.setxBP();
	#endif
	//!energy
	#ifndef _ARK_ISOTHERMAL_FLIP
	g.setxBP();
	#endif
}
void Cell_2D::FactorAC()
{
	#ifdef _ARK_ALLENCAHN_FLIP
	h.setxBP();
	#endif
}
void Cell_2D::SetVolume()
{
	if(3 == celltype) 
	{
		xc = (cellNodes[0]->xN + cellNodes[1]->xN + cellNodes[2]->xN)/3.0;
		yc = (cellNodes[0]->yN + cellNodes[1]->yN + cellNodes[2]->yN)/3.0;
		volume = TriArea(
						 cellNodes[0]->xN,cellNodes[0]->yN,
						 cellNodes[1]->xN,cellNodes[1]->yN,
						 cellNodes[2]->xN,cellNodes[2]->yN
						);
	}
	else if(4 == celltype)
	{
//		
		double xc_A,yc_A,xc_B,yc_B,volume_A,volume_B;
//		
		xc_A = (cellNodes[0]->xN + cellNodes[1]->xN + cellNodes[2]->xN)/3.0;
		yc_A = (cellNodes[0]->yN + cellNodes[1]->yN + cellNodes[2]->yN)/3.0;
		volume_A = TriArea(
							cellNodes[0]->xN,cellNodes[0]->yN,
							cellNodes[1]->xN,cellNodes[1]->yN,
							cellNodes[2]->xN,cellNodes[2]->yN
						  );
//
		xc_B = (cellNodes[0]->xN + cellNodes[3]->xN + cellNodes[2]->xN)/3.0;
		yc_B = (cellNodes[0]->yN + cellNodes[3]->yN + cellNodes[2]->yN)/3.0;
		volume_B = TriArea(
							cellNodes[0]->xN,cellNodes[0]->yN,
							cellNodes[3]->xN,cellNodes[3]->yN,
							cellNodes[2]->xN,cellNodes[2]->yN
						  );
//		
		volume = volume_A + volume_B;
//
		xc = (volume_A*xc_A + volume_B*xc_B)/volume;
		yc = (volume_A*yc_A + volume_B*yc_B)/volume;
	}
	DtSlashVolume = dt/volume;
}
//
Face_2D::Face_2D():msqh(new MacroQuantity()),use(new int(1))
{
	AllocateARK(xi_n_dS,Qu,Qv);
}
Face_2D::~Face_2D()
{
//
	if(--*use == 0)	
	{
		DeallocateARK(xi_n_dS,Qu,Qv);
//
		delete msqh;
		delete use;
		//cout << "~Face_2D()" <<'\n';
	}
}
Face_2D::DVDF::DVDF()
{
	AllocateARK(hDt,Face_2D::Qu,Face_2D::Qv);
	AllocateARK(BhDt,Face_2D::Qu,Face_2D::Qv);
	AllocateARK(EqhDt,Face_2D::Qu,Face_2D::Qv);
	AllocateARK(SohDt,Face_2D::Qu,Face_2D::Qv);
}
Face_2D::DVDF::~DVDF()
{
	DeallocateARK(hDt,Qu,Qv);
	DeallocateARK(BhDt,Qu,Qv);
	DeallocateARK(EqhDt,Qu,Qv);
	DeallocateARK(SohDt,Qu,Qv);
	//cout << "Face_2D::DVDF::~DVDF()" <<'\n';
}
void Face_2D::DVDF::setxh()
{
	ah = 2.0*tauh/(2.0*tauh+::hDt);
	bh = 1.0 - ah;
	ch = tauh*bh;
}
void Face_2D::Factor()
{
	//!momentum
	#ifdef _ARK_MOMENTUM_FLIP
	f.setxh();
	#endif
	//energy
	#ifndef _ARK_ISOTHERMAL_FLIP
	g.setxh();
	#endif
}
void Face_2D::FactorAC()
{
	#ifdef _ARK_ALLENCAHN_FLIP
	h.setxh();
	#endif
}
void Face_2D::SetArea()
{
	xf = 0.5*(faceNodes[0]->xN + faceNodes[1]->xN);
	yf = 0.5*(faceNodes[0]->yN + faceNodes[1]->yN);
	Area = sqrt(
				(faceNodes[1]->yN - faceNodes[0]->yN) * (faceNodes[1]->yN - faceNodes[0]->yN)
				 +
				(faceNodes[1]->xN - faceNodes[0]->xN) * (faceNodes[1]->xN - faceNodes[0]->xN)
			   );
}
void Face_2D::SetNormalV()
{
	double dy = (faceNodes[1]->yN - faceNodes[0]->yN), 
	 	   dx = (faceNodes[1]->xN - faceNodes[0]->xN);
	SetZero(dx);SetZero(dy); 	 	   
	Vx = dy/Area;Vy = -dx/Area;
}