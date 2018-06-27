#ifndef _ARK_D2Q9_H
#define _ARK_D2Q9_H
//
#define QuIndex j
//
#include <cmath>

int const

DV_Qu = 1,

DV_Qv = 9;

double const

MaxU = sqrt(3.0*RT),

Eta = 0,

MaSpan = sqrt(3.0*RT);

namespace D2Q9{

double const xi_u[DV_Qv] = {0,1,1,0,-1,-1,-1,0,1};

double const xi_v[DV_Qv] = {0,0,1,1,1,0,-1,-1,-1};

int const _BB[DV_Qv] = {0,5,6,7,8,1,2,3,4};

}

#endif
// const double RT = 1.0/3.0,s_RT = sqrt(3.0*RT);
// //
// const typexi xi[Q]={{0,0},{s_RT,0},{s_RT,s_RT},{0,s_RT},
// 					{-s_RT,s_RT},{-s_RT,0},{-s_RT,-s_RT},{0,-s_RT},{s_RT,-s_RT}};
// //
// const double omega[Q]={4.0/9.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0,
// 						1.0/9.0, 1.0/36.0};					
// #endif
