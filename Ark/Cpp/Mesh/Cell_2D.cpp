#include "Mesh_2D.h"

void Cell_2D::assign(const Cell_2D &rhs)
{
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
}
Cell_2D::Cell_2D():msq(new MacroQuantity()),use(new int(1))
{
}
Cell_2D::Cell_2D(const Cell_2D &rhs)
{
	assign(rhs);
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
		assign(rhs);
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
//!----------------------------DVDF----------------------------
//!---------------------------private--------------------------
void Cell_2D::DVDF::allocate()
{
	AllocateARK(BarP,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(BarP_x,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(BarP_y,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(Tilde,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(Eq,Cell_2D::Qu,Cell_2D::Qv);
	AllocateARK(So,Qu,Qv);
}
void Cell_2D::DVDF::deallocate()
{
	DeallocateARK(BarP,Cell_2D::Qu,Cell_2D::Qv);
	DeallocateARK(BarP_x,Cell_2D::Qu,Cell_2D::Qv);
	DeallocateARK(BarP_y,Cell_2D::Qu,Cell_2D::Qv);
	DeallocateARK(Tilde,Cell_2D::Qu,Cell_2D::Qv);
	DeallocateARK(Eq,Cell_2D::Qu,Cell_2D::Qv);
	DeallocateARK(So,Qu,Qv);
}
void Cell_2D::DVDF::assign(const DVDF& rhs)
{
	BarP = rhs.BarP;
	BarP_x = rhs.BarP_x;
	BarP_y = rhs.BarP_y;
	Tilde  = rhs.Tilde;
	Eq  = rhs.Eq;
	So  = rhs.So;
}
//----------------------------------constructors-------------------------
Cell_2D::DVDF::DVDF():token(new int(1))
{
	allocate();
}
Cell_2D::DVDF::DVDF(const DVDF& rhs)
{
	assign(rhs);
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
			deallocate();
		//
			delete token;
		}
		assign(rhs);
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
		deallocate();
	}
}
//-------------------------------------------------------------------
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