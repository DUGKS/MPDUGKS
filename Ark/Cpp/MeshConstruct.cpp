#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>
#include "MeshConstructFunc.h"

using std::cout;
using std::cin;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;

void EraseParenthese(string &string_line,int &body)
{
	//only process lines whose first char is a left bracket;
	if('(' != string_line.front()) 
	{
		body = -1;
		return;
	}
	//pop space
	//while(' ' == string_line.back() || '\n' == string_line.back() || '\r' == string_line.back())
	while(std::isspace(string_line.back()))
		string_line.pop_back();
	//
	if(string_line.back() == '(') ++body;
	auto It = remove_if(string_line.begin(),string_line.end(),
									[](const char& c){return (c == '(' || c == ')');});
	//
	if(It == string_line.end()) body = -1;
	string_line.erase(It,string_line.end());
}
int StringToHex(int* const &ptrHexLine, istringstream &iss_line,int &count)
{
	for(;iss_line >> std::hex >>ptrHexLine[count];++count);
	return 0;
}
void HeadProcess(const int &index,int* const &ptrHexLine,int const &count)
{
	Allocate_Mesh(index,ptrHexLine,count);
}
void BodyProcess(const int &index,ifstream& InFile_Mesh,int* const &ptrHexLine)
{
	if(index == 10)
	{		
		ConstructNodes(ptrHexLine,InFile_Mesh);
	}
	else if(index == 12)
	{
		ConstructCells(ptrHexLine,InFile_Mesh);
	}
	else if(index == 13)
	{
		ConstructFaces(ptrHexLine,InFile_Mesh);
	}
	else if(index == 18)
	{
		ConstructPeriodicFaces(ptrHexLine,InFile_Mesh);
	}
	else
	{
		cout << "Invalid index during body processing" << endl;
		getchar();
	}
}
int MeshConstruct(const string &s)
{
	ifstream InFile_Mesh;
	const char *phome = std::getenv("HOME");
	string shome(phome);
	shome += "/Mesh/";
	InFile_Mesh.open(shome + s +".cas");
	if(!InFile_Mesh)
	{
		cout <<"----------------fatal error!!!---------------------"<<endl;
		cout << "Mesh file open failed!!!  "<<endl;
		cout <<"name of mesh file : "<<shome + s +".cas"<<endl;
		cout <<__FILE__ <<" : "<<__LINE__<<"  "<<__func__<<endl; 
		exit(0);
	}
//
	int *ptrHexLine = new int[MeshPerLine];
	int line = 0, count = 0, index = 0,body = 0,FaceCount = 0;
	string string_line;
	while(getline(InFile_Mesh,string_line))
	{
		if(string_line.size() == 0) continue;
		cout << string_line << endl;
		body = 0;index = 0;	count = 0;
		for(int i = 0;i != MeshPerLine;++i)
			ptrHexLine[i] =  0;
		EraseParenthese(string_line,body);
		if(-1 == body) continue;		
		istringstream iss_line(string_line);
		if(!(iss_line >> index) || !(index == 10 || index == 12 || index == 13 || index == 18))
			continue;
		StringToHex(ptrHexLine,iss_line,count);
		if(0 == body)
			HeadProcess(index,ptrHexLine,count);
		else if(1 == body)
		{
			BodyProcess(index,InFile_Mesh,ptrHexLine);
		}
		else
		{
			cout << "unknown body indicator"<<endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(0);
		}
	}
	cout << "Faces " << Faces << " Nodes " << Nodes << " Cells " << Cells <<endl;
	delete []ptrHexLine;
	return 0;
}

void AllocateFaces(int BoundFaceNum, Face_2D** &ptrBoundFace, vector<Face_2D*> &vecBoundFace)
{
	if (0 == BoundFaceNum) return;
	ptrBoundFace  = new Face_2D*[BoundFaceNum];
	for(int i = 0;i != BoundFaceNum;++i)
		ptrBoundFace[i] = vecBoundFace[i];
}
void FacesClassify()
{
	vector<Face_2D*> InteriorVec,WallVec,PeriodicVec,P_InletVec,P_OutletVec,V_InletVec,
					 SymmetryVec,P_FarfieldVec,BoundVec;//2,3,12-8,4,5,10,7,9;
	for(int i = 0;i != Faces; ++i)
	{
		if(2 == FaceArray[i].bc_type)
			InteriorVec.push_back(&FaceArray[i]);
		else if(12 == FaceArray[i].bc_type || 8 == FaceArray[i].bc_type)
		{
			PeriodicVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(3 == FaceArray[i].bc_type)
		{
			WallVec.push_back(&FaceArray[i]);
			BoundVec.push_back(FaceArray+i);
		}
		else if(4 == FaceArray[i].bc_type)
		{
			P_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(5 == FaceArray[i].bc_type)
		{
			P_OutletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(7 == FaceArray[i].bc_type)
		{
			SymmetryVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(9 == FaceArray[i].bc_type)
		{
			P_FarfieldVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else if(10 == FaceArray[i].bc_type)
		{
			V_InletVec.push_back(FaceArray+i);
			BoundVec.push_back(FaceArray+i);
		}
		else
		{
			cout << "FaceArray : " << i <<" unknown bc_type : " << FaceArray[i].bc_type
			<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
			getchar();
		}
	}
	#ifdef _ZERO_NDEBUG_FLIP
	if(InteriorFaceNum != InteriorVec.size())
	{
		cout << "FaceArray : " << InteriorFaceNum <<" interior Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(WallFaceNum != WallVec.size())
	{
		cout << "FaceArray : " << WallFaceNum <<" Wall Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(PeriodicFaceNum != PeriodicVec.size())
	{
		cout << "FaceArray : " << PeriodicFaceNum <<" Periodic Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_InletFaceNum != P_InletVec.size())
	{
		cout << "FaceArray : " << P_InletFaceNum <<" Pressure Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_OutletFaceNum != P_OutletVec.size())
	{
		cout << "FaceArray : " << P_OutletFaceNum <<" Pressure Outlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(V_InletFaceNum != V_InletVec.size())
	{
		cout << "FaceArray : " << V_InletFaceNum <<" Velocity Inlet Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(SymmetryVec.size() != SymmetryFaceNum)
	{
		cout << "FaceArray : " << SymmetryFaceNum <<" Symmetry Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(P_FarfieldVec.size() != P_FarfieldFaceNum)
	{
		cout << "FaceArray : " << P_FarfieldFaceNum <<" Pressure Farfield Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	if(BoundVec.size() != BoundFaceNum)
	{
		cout << "FaceArray : " << BoundFaceNum <<" Boundary Faces doesn't equal : "
		<<"  "<<__FILE__<<"  "<<__LINE__<<"  "<<__func__<< endl;
		getchar();
	}
	cout <<"Faces : "<<Faces<<'\n'
		 <<"BoundFaceNum : "<<BoundFaceNum<<'\n'
		 <<"InFaceNum : "<<InteriorFaceNum<<'\n'
		 <<"WallFaceNum : "<<WallFaceNum<<'\n'
		 <<"PeriodicFaceNum : "<<PeriodicFaceNum<<'\n'
		 <<"P_InletFaceNum : "<<P_InletFaceNum<<'\n'
		 <<"P_OutletFaceNum : "<<P_OutletFaceNum<<'\n'
		 <<"SymmetryFaceNum : "<<SymmetryFaceNum<<'\n'
		 <<"P_FarfieldFaceNum : "<<P_FarfieldFaceNum<<'\n'
		 <<"V_InletFaceNum : "<<V_InletFaceNum<<endl;
	if(InteriorFaceNum + WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != Faces)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "Faces != InteriorFaceNum + WallFaceNum"<<endl;
		getchar();
	}
	if(WallFaceNum + PeriodicFaceNum + P_InletFaceNum + P_OutletFaceNum
		+ SymmetryFaceNum + P_FarfieldFaceNum + V_InletFaceNum != BoundFaceNum)
	{
		_PRINT_ERROR_MSG_FLIP
		cout << "BoundFaceNum != ..."<<endl;
		getchar();
	}
	#endif
	AllocateFaces(InteriorFaceNum,InteriorFaceA,InteriorVec);
	AllocateFaces(PeriodicFaceNum,PeriodicFaceA,PeriodicVec);
	AllocateFaces(WallFaceNum,WallFaceA,WallVec);
	AllocateFaces(P_InletFaceNum,P_InletFaceA,P_InletVec);
	AllocateFaces(P_OutletFaceNum,P_OutletFaceA,P_OutletVec);
	AllocateFaces(SymmetryFaceNum,SymmetryFaceA,SymmetryVec);
	AllocateFaces(P_FarfieldFaceNum,P_FarfieldFaceA,P_FarfieldVec);
	AllocateFaces(V_InletFaceNum,V_InletFaceA,V_InletVec);
	AllocateFaces(BoundFaceNum,BoundFaceA,BoundVec);
}

//index of diagonal cell--index of neighbour cell-------index of the neighbour faces
/*
---------------				----------------				 ------
| D1 |   | D0 |				|    | C1 |    |			     | F1 |     
---------------				----------------			----------------
|    | C |    |				| C2 |  C | C0 |			  F2 |  C | F0  
---------------				----------------			----------------
| D2 |   | D3 |				|    | C3 |    |			     | F3 |     
---------------				----------------				 ------
*/
void sortFacesInThisCell(Cell_2D &cell)
{
	Face_2D *Facetmp[4] = {nullptr,nullptr,nullptr,nullptr};
	for(int i = 0;i < cell.celltype;++i)
	{
		if(fabs(cell.Face_C[i]->xf - cell.xc) > infinitesimal)
		{
			if(cell.Face_C[i]->xf > cell.xc)
			{
				Facetmp[0] = cell.Face_C[i];
			}
			else
			{
				Facetmp[2] = cell.Face_C[i];
			}
		}
		else
		{
			if(cell.Face_C[i]->yf > cell.yc)
			{
				Facetmp[1] = cell.Face_C[i];
			}
			else
			{
				Facetmp[3] = cell.Face_C[i];
			}
		}
	}
	for(int i = 0;i < cell.celltype;++i)
	{
		cell.Face_C[i] = Facetmp[i];
	}
}
void NeighbourCellConstruct()
{
	for(int n = 0;n < Cells;++n)
	{
		#ifdef _CARTESIAN_MESH_FLIP
			sortFacesInThisCell(CellArray[n]);
		#endif
		for(int i = 0;i < CellArray[n].celltype;++i)
		{
		CellArray[n].Cell_C[i] = ((CellArray[n].Face_C[i] -> owner ==  (&CellArray[n])) ? 
						CellArray[n].Face_C[i] -> neigh : CellArray[n].Face_C[i] -> owner);
		CellArray[n].signFlux[i] = ((CellArray[n].Face_C[i] -> owner ==  (&CellArray[n])) ? -1 : 1);
		}
	}
}
void DiagonalCellConstruct()
{
	int ne = 0,nw = 0,se = 0,sw = 0;
	for(int n = 0;n < Cells;++n)
	{
	//---------------------------------Diagonal 0 3------------------------
		if(nullptr == CellArray[n].Cell_C[0]->ShadowC)
		{
			CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[0]->Cell_C[1];
			CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[0]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[0] = CellArray[n].Cell_C[1]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[0] = PeriodicShadowC_NE;
				++ne;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[3] = CellArray[n].Cell_C[3]->Cell_C[0];
			}
			else
			{
				CellArray[n].Cell_Diag[3] = PeriodicShadowC_SE;
				++se;
			}
		}
	//---------------------------------Diagonal 1 2------------------------
		if(nullptr == CellArray[n].Cell_C[2]->ShadowC)
		{
			CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[2]->Cell_C[1];
			CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[2]->Cell_C[3];
		}
		else
		{
			if(nullptr == CellArray[n].Cell_C[1]->ShadowC)
			{
				CellArray[n].Cell_Diag[1] = CellArray[n].Cell_C[1]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[1] = PeriodicShadowC_NW;
				++nw;
			}
			if(nullptr == CellArray[n].Cell_C[3]->ShadowC)
			{
				CellArray[n].Cell_Diag[2] = CellArray[n].Cell_C[3]->Cell_C[2];
			}
			else
			{
				CellArray[n].Cell_Diag[2] = PeriodicShadowC_SW;
				++sw;
			}
		}
	}
}
void CarfaceCellsConstruct()
{
	LoopPS(Faces)
	{
		Face_2D &face = FaceArray[n];
		if(EqualZero(FaceArray[n].Vx))
		{
			if(0 != face._dx)
			{
				cout <<"Collision : Vx == 0 && face dx != 0 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
			face._dx = 2*(face.faceNodes[1]->xN - face.faceNodes[0]->xN);
			if(doubleEqual(FaceArray[n].Vy,1))
			{
				if(face.neigh != face.owner->Cell_C[1])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[0];
				face.faceCells[1] = face.owner->Cell_Diag[0];
				face.faceCells[2] = face.owner->Cell_C[2];
				face.faceCells[3] = face.owner->Cell_Diag[1];
			}
			else if(doubleEqual(FaceArray[n].Vy,-1))
			{
				if(face.neigh != face.owner->Cell_C[3])
				{
					_PRINT_SPLITLINE_ARK
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					cout <<"face xf : "<<face.xf<<"    "<<"face.yf : "<<face.yf<<endl;
					cout <<"face Vx : "<<face.Vx<<"    "<<"face.Vy : "<<face.Vy<<endl;
					cout <<face.neigh->xc - face.owner->xc <<"    "
						 <<face.neigh->yc - face.owner->yc<<endl;
					cout <<"face neigh : "<<face.neigh<<"    "<<"face.owner : "<<face.owner<<endl;
					for(int k = 0;k < 4;++k)
					{
						cout <<"Cell_C["<<k<<"] : "<<face.owner->Cell_C[k]<<"    ";
					}
					cout << endl;
					_PRINT_ERROR_MSG_FLIP
					_PRINT_SPLITLINE_ARK
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[2];
				face.faceCells[1] = face.owner->Cell_Diag[2];
				face.faceCells[2] = face.owner->Cell_C[0];
				face.faceCells[3] = face.owner->Cell_Diag[3];
			}
			else
			{
				cout <<"Vy != 1 && Vy != -1 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
		}
		else if(EqualZero(FaceArray[n].Vy))
		{
			if(0 != face._dy)
			{
				cout <<"Collision : Vx == 0 && face dx != 0 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
			face._dy = 2*(face.faceNodes[1]->yN - face.faceNodes[0]->yN);
			if(doubleEqual(FaceArray[n].Vx,1))
			{
				if(face.neigh != face.owner->Cell_C[0])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[3];
				face.faceCells[1] = face.owner->Cell_Diag[3];
				face.faceCells[2] = face.owner->Cell_C[1];
				face.faceCells[3] = face.owner->Cell_Diag[0];
			}
			else if(doubleEqual(FaceArray[n].Vx,-1))
			{
				if(face.neigh != face.owner->Cell_C[2])
				{
					cout <<"Collision : neigh and owner doesn't match "<<endl;
					_PRINT_ERROR_MSG_FLIP
					getchar();
					exit(-1);
				}
				face.faceCells[0] = face.owner->Cell_C[1];
				face.faceCells[1] = face.owner->Cell_Diag[1];
				face.faceCells[2] = face.owner->Cell_C[3];
				face.faceCells[3] = face.owner->Cell_Diag[2];
			}
			else
			{
				cout <<"Vx != 1 && Vx != -1 "<<endl;
				_PRINT_ERROR_MSG_FLIP
				getchar();
				exit(-1);
			}
		}
		else
		{
			cout <<"Attempting to construct faceCells for non Cartesian Mesh "<<endl;
			_PRINT_ERROR_MSG_FLIP
			getchar();
			exit(-1);
		}
	}
}
void AllocateCarCellArray()
{
	const double MinDx = dx, MinDy = dy;
//
	CarCellArray = new Cell_2D** [Nxp2];
	for(int i = 0;i < Nxp2;++i)
		CarCellArray[i] = new Cell_2D* [Nyp2];
	for(int i = 1;i < Nxp1;++i)
	for(int j = 1;j < Nyp1;++j)
	{
		CarCellArray[i][j]=0;
	}
//
	for(int n = 0;n < Cells;++n)
	{
		int i = std::round((CellArray[n].xc - X_Beg + 0.5*MinDx)/MinDx);
		int j = std::round((CellArray[n].yc - Y_Beg + 0.5*MinDy)/MinDy);
		CarCellArray[i][j] = &CellArray[n];
	}
	#ifdef _PERIODIC_12_8_BCs_FLIP
	for(int k = 0;k < PeriodicFaceNum;++k)
	{
		int i = std::round((PeriodicShadowCA[k].xc - X_Beg + 0.5*MinDx)/MinDx);
		int j = std::round((PeriodicShadowCA[k].yc - Y_Beg + 0.5*MinDy)/MinDy);
		CarCellArray[i][j] = &PeriodicShadowCA[k];
	}
	#endif
	#ifdef _Wall_3_BCs_FLIP
	for(int k = 0;k < WallFaceNum;++k)
	{
		int i = std::round((WallShadowCA[k].xc - X_Beg + 0.5*MinDx)/MinDx);
		int j = std::round((WallShadowCA[k].yc - Y_Beg + 0.5*MinDy)/MinDy);
		CarCellArray[i][j] = &WallShadowCA[k];
	}
	#endif
	CarCellArray[0][0] = PeriodicShadowC_SW;
	CarCellArray[0][Nyp1] = PeriodicShadowC_NW;
	CarCellArray[Nxp1][0] = PeriodicShadowC_SE;
	CarCellArray[Nxp1][Nyp1] = PeriodicShadowC_NE;
}
void DeallocateCarCellArray()
{
	for(int i = 0;i < Nxp2;++i)
	{
		delete[] CarCellArray[i];
		CarCellArray[i] = nullptr;
	}
	delete[] CarCellArray;
}
void SetFace_dxdy()
{
	for(int n = 0;n < Faces;++n)
	{
		FaceArray[n]._dx = FaceArray[n].owner->xc - FaceArray[n].neigh->xc;
		FaceArray[n]._dy = FaceArray[n].owner->yc - FaceArray[n].neigh->yc; 
		SetZero(FaceArray[n]._dx);
		SetZero(FaceArray[n]._dy);
	}
}

void setXiDotdS()
{
	for(int n = 0;n < Faces;++n)
	for(int i = 0;i < DV_Qu;++i)
	for(int j = 0;j < DV_Qv;++j)
	{
		FaceArray[n].xi_n_dS[i][j] = FaceArray[n].Area*((xi_u[QuIndex])*FaceArray[n].Vx + xi_v[j]*FaceArray[n].Vy);
	}
}