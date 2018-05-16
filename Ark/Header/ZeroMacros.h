#ifndef _ZERO_MACRO_ARK
#define _ZERO_MACRO_ARK

#ifndef _PRINT_ERROR_MSG_FLIP
#define _PRINT_ERROR_MSG_FLIP  cout<<"File : "<<__FILE__<<"  Line : "\
<<__LINE__<<"  fun : "<<__func__<<'\n';
#endif

#ifndef _PRINT_SPLITLINE_ARK
#define _PRINT_SPLITLINE_ARK	cout<<"----------------------------------"<<'\n';
#endif

#define LoopPS(MESHNUM) for(int n=0;n<MESHNUM;++n)

#define LoopVS(NUMQU,NUMQV)		\
for(int i=0;i<NUMQU;++i)		\
for(int j=0;j<NUMQV;++j)

#define nl '\n'

#define fs "    "

#define Info cout<<setiosflags(ios::scientific)<<setprecision(12)

double const

infinitesimal = 1.0E-14,

PI = 3.141592653589793;

namespace ARK{

int const 

ND = 2,

cellN = 4,

cellF = 4,

faceN = 2,

faceF = 0;

}

#endif