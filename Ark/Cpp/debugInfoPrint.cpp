#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include "DUGKSDeclaration.h"

using std::cout;
using std::endl;
using std::ofstream;
using std::ostringstream;
using std::string;
using std::ios;
using std::setprecision;
using std::setiosflags;
//
void printSplitLine(char c)
{
	cout <<"----------------------------------"<<c<<endl;
}
void printErrorMessage(int line,const char *file,const char *func)
{
	cout<<"File : "<<file<<"  Line : "<<line<<"  fun : "<<func<<'\n';
}