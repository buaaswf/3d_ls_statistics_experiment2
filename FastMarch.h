/**************************************************************************
	ORIGINAL AUTHOR: 
		Emud Mokhberi (emud@ucla.edu)
	MODIFIED BY:
	
	CONTRIBUTORS:
		

-----------------------------------------------
	
 ***************************************************************
 ******General License Agreement and Lack of Warranty ***********
 ****************************************************************

 This software is distributed for noncommercial use in the hope that it will 
 be useful but WITHOUT ANY WARRANTY. The author(s) do not accept responsibility
 to anyone for the consequences of using it or for whether it serves any 
 particular purpose or works at all. No guarantee is made about the software 
 or its performance.

 You are allowed to modify the source code, add your name to the
 appropriate list above and distribute the code as long as 
 this license agreement is distributed with the code and is included at 
 the top of all header (.h) files.

 Commercial use is strictly prohibited.
***************************************************************************/

/*
	FastMarch2D : A class for resetting the grid values to be the signed distance from the interface
	Inputs: grid and cell size
			
	The class maintains its own special grid of type FMContainer and in addition to two lists, one
	of ClosePoints which is used to determine the initial set of close points and a second list for
	maintaining the min heap which is used to determine the signed distance values
	
	The class works by first resetting all the negative signed distance values and then setting all
	of the positive signed ditance values 

	Public Functions:
	Set				- initialized the value of the FMContainer grid for the specified index
	Reinitialize	- Performs the fastmarching method on the FMContainer grid. It is assumed
					  that the values to be reset are in the FMContainer grid and that is where
					  the updated grid values will be when the function is done.
	
	Private Functions:
	PopHeap			- Called by FastMarch, it pops the current closest grid value from the heap 
					  and cosidered that cell value to be done. It then restores the heap to be a min heap
	FastMarch		- Called by ReinitHalf. While there are close points in the heap, it pops the heap, 
					  determines if the threshold for the fastmarching has been met (we only need to restore
					  the signed distance function to within a finite number of cell away from the
					  interface) and if it hasn't, it updates the neighboring cells if they are not
					  considered done points
	AddToHeap		- Called by FindPhi. Adds a former far point which has become a close point to the 
				      heap and updates the heap to remain a min heap
	UpdateHeap		- Called by FindPhi. This function is called if a close point's (which is in the heap) 
					  value is changed as a result of a neighboring cell being 'done'. It updates the heap 
					  to remain a min heap
	CheckMax2		- Called by FindPhi. This function is used if a close point has two neighboring done
					  points. It is used to make sure that the distance between the two neighboring done
					  points is not greater than the cell size. If it is, it only uses the smaller of the
					  two values in updating the signed distance of the close point
	CheckMax3		- Called by FindPhi. This function is used if a close point has three neighboring done
					  points. It is used to make sure that the distance between the three neighboring done
					  points is not greater than the cell size. If it is, it only uses the two smaller
					  values in updating the signed distance of the close point
	CheckFront		- Called by FindPhi. This function checks to see if the cell in front of the current
					  close point is done
	CheckBehind		- Called by FindPhi. This function checks to see if the cell behind of the current
					  close point is done. This function also handles the case where a close points has
					  done points on both sides along one axis
	FindPhi			- Called by FastMarch. This function updates the value os the specified cell using 
					  the values of neighboring 'done' cells. If the cell is already a close point, it 
					  calls UpdateHeap. Otherwise, it changes the cell to be a close point and calls
					  AddToHeap.
	SetBoundary		- Called by ReinitHalf. Sets the boundary grid done flags to -1
	ReinitHalf		- Called by Reinitialize. Resets either the exterior or interior of the interface
					  to be a signed distance function.
	Initialize		- Called by ReinitHalf. Steps through the grid along each axis and determines which 
					  cells are 'done', and which are 'close'. It uses linear interpolation to set the 
					  values of those cells adjacent to the interface.
	InitHeap		- Called by ReinitHalf. After Initialize has run and determined the first set of
					  'close' points, this function will create a min heap of those points. This is
					  done by stepping through the ClosePoints list and calling FindPhi for each cell
					  that isn't already in the heap
	AddClose		- Adds cell to ClosePoints list.
					  

	Created by Emud Mokhberi: UCLA : 09/04/04
*/

#ifndef FastMarch_H
#define FastMarch_H

#include "global.h"
#include <vector>
using std::vector;

class FastMarch
{
public:
	FastMarch(float x, float y, float z);
	~FastMarch() {freeMem();}
	struct FMContainer {
		BYTE DoneFlag;
		int HeapPosition;
		float value;	// distance value
		float grad[3];	// the gradient of the distance transform
	};


	inline float operator[] (int index) { FMContainer* temp=NULL; if (mapDist.Lookup(index, temp)) return temp->value; else {ASSERT(FALSE); return float(gPara.FMM_LIMIT_P);} }
	inline float operator() (int i, int j, int k) { FMContainer* temp=NULL; if (mapDist.Lookup(GI(i,j,k), temp)) return temp->value; else {ASSERT(FALSE); return float(gPara.FMM_LIMIT_P);} }

	inline BOOL	grad(int index, float*& v)
	{
		FMContainer* temp = NULL;
		if (mapDist.Lookup(index, temp))
		{
			v = temp->grad;
			return TRUE;
		} else return FALSE;
	}
	inline BOOL	grad(int i, int j, int k, float v[3]) {FMContainer* temp=NULL; int index = GI(i,j,k); if (!mapDist.Lookup(index, temp)) return FALSE; else {if (temp->DoneFlag == 10) {memcpy(v,temp->grad,3*sizeof(float)); return TRUE;} else return FALSE;}}

    void Reinitialize(BYTE* gridTYpe);

	POSITION sPos, cPos;
	void gradHalf(bool bFst=false);
	void smooth(BYTE* type);

protected:
	inline float dist(int index) 
	{// the fast version
		FMContainer* temp=NULL; 
		if (mapDist.Lookup(index, temp))
			return temp->value;
		else 
			return float(gPara.FMM_LIMIT_P);
	}

private:
	int PopHeap(FMContainer*& data);
    void March();
    void AddToHeap(FMContainer* data, int index);
	void UpdateHeap(FMContainer* data, int index);
    inline void CheckMax2(int& a, bool& flag, float& phi1,
						  const float &phi2, const float &inv2);
    inline void CheckMax3(int& a, bool& flag, float& phi1, 
                          const float &phi2, const float &phi3, 
						  const float &inv2, const float &inv3);
    inline void CheckFront(float& phi, int& a, bool& flag, int index);
	inline void CheckBehind(float& phi, int& a, bool& flag, int index);
    void FindPhi(FMContainer* data, int index, int x, int y, int z);
    
	void ReinitPositiveHalf(BYTE* gridTYpe);
	void ReinitNegativeHalf(BYTE* gridTYpe);

	void InitPositive(BYTE* gridTYpe);
	void InitNegative(BYTE* gridTYpe);

    inline void AddClose(float f, int index);
	inline void freeMem();
	void InitHeap();
	inline void CalcGradient(FMContainer* data, int x, int y, int z);

	// The real unit length in xyz-direction of each voxel, in mm, and their squares.
	float fx, fy, fz, fxx, fyy, fzz, fxxINV, fyyINV, fzzINV, sumINV;

	int heapSize;
	int closeSize;
	CMap<int, int, FMContainer*, FMContainer*&> mapDist;
	vector<int> FMHeap;
	vector<int> ClosePoints;

	float	CurrenLimit;
};

#endif