
#include "FastMarch.h"
#include "global.h"

#define FM_INIT_SIZE	30000000
FastMarch::FMContainer* cd = NULL;
float sq = 0;

FastMarch::FastMarch(float x, float y, float z)
	: fx(x), fy(y), fz(z), fxx(x*x), fyy(y*y), fzz(z*z),
	fxxINV(1/fxx), fyyINV(1/fyy), fzzINV(1/fzz), sumINV(fxxINV+fyyINV+fzzINV)
{
	mapDist.InitHashTable(FM_INIT_SIZE);
    
	heapSize = 0;
	closeSize = 0;

	CurrenLimit = 0;
}

void FastMarch::Reinitialize(BYTE* gridType)
{
	freeMem();

	// initialize the initial layer
	FOR_EACH if (gridType[i] & VXL_TISSUE_BORDER)
	{
		FMContainer* data = new FMContainer;
		data->DoneFlag = 1;
		data->HeapPosition = -1;
		data->value = 0.f;		// the distance value of the mucosa layer should be zero.
		memset(data->grad, 0, 3*sizeof(float));
		mapDist.SetAt(i, data);
	}

	FMHeap.resize(FM_INIT_SIZE);
	ClosePoints.resize(FM_INIT_SIZE);

    // Negative Phi first
	CurrenLimit = gPara.FMM_LIMIT_N;
    ReinitNegativeHalf(gridType);

	// swap the sign of each negative distance
	POSITION pos = mapDist.GetStartPosition();
	while (pos)
	{
		mapDist.GetNextAssoc(pos, i, cd);
		if ( gridType[i]&(VXL_AIR|VXL_MUCOSA) && !(gridType[i]&VXL_TISSUE_BORDER) )
			cd->value = -cd->value;
	}

    // Then Positive Phi
	CurrenLimit = gPara.FMM_LIMIT_P;
	ReinitPositiveHalf(gridType);

	ClosePoints.resize(0);
	FMHeap.resize(0);

	// 
	pos = mapDist.GetStartPosition();
	sPos = pos; // save the start position, used for multi-thread processing
	int cnt = 1+mapDist.GetCount()/2, ii = 0;
	while (ii < cnt)
	{
		mapDist.GetNextAssoc(pos, i, cd);
		ii ++;
	} cPos = pos;
}

void FastMarch::gradHalf(bool bFst)
{
	POSITION s, e;
	if (bFst)	{s = sPos; e = cPos;}
	else		{s = cPos; e = NULL;}
	int ig, ii, jj, kk;
	while (s != e)
	{
		mapDist.GetNextAssoc(s, ig, cd);
		if (cd->DoneFlag == 1)
		{
			GIJK(ig, ii, jj, kk);
			CalcGradient(cd, ii, jj, kk);
		}
	}

}

void FastMarch::smooth(	BYTE* type )
{
	// build kernel
	float x, y, z, w;

	// smoothing......
	FMContainer* cd = NULL, *cd0 = NULL;
	int ig = 0, igC = 0;
	POSITION pos = mapDist.GetStartPosition();
	while (pos != NULL)
	{
		mapDist.GetNextAssoc(pos, ig, cd);
		if (type[ig] & VXL_MUCOSA)
		{
			GIJK(ig, i, j, k);
			w = sqrt(cd->grad[0]*cd->grad[0] + cd->grad[1]*cd->grad[1] + cd->grad[2]*cd->grad[2]);
			if (w < 0.00001)
			{
				x = y = z = 0.;

				igC = GI(i+1, j, k); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				igC = GI(i-1, j, k); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				igC = GI(i, j+1, k); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				igC = GI(i, j-1, k); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				igC = GI(i, j, k+1); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				igC = GI(i, j, k-1); mapDist.Lookup(igC, cd0); x+=cd0->grad[0]; y+=cd0->grad[1]; z+=cd0->grad[2];
				
				x /= 6; y /= 6; z /= 6;
				w = sqrt(x*x + y*y + z*z);
				cd->grad[0] = x/w; cd->grad[1] = y/w; cd->grad[2] = z/w;
			}
		}
	}
}

/////////////////////////////////////////////////////////////
// Using sobel operator to compute the gradient for each voxel in 
// the mucosa volume area based on the distance transform values
// 
/////////////////////////////////////////////////////////////
inline void FastMarch::CalcGradient(FMContainer* data, int x, int y, int z)
{
	static float scale = 1./12.;

	data->grad[0]=(  3*(dist(PPZ)-dist(MMZ)) - 3*(dist(MPZ)-dist(PMZ)) + 6*(dist(PZZ)-dist(MZZ))
					+3*(dist(PZP)-dist(MZM)) + 3*(dist(PZM)-dist(MZP)) + 2*(dist(PPP)-dist(MMM))
					+2*(dist(PMP)-dist(MPM)) - 2*(dist(MMP)-dist(PPM)) - 2*(dist(MPP)-dist(PMM)) ) * scale;
	data->grad[1]=(  3*(dist(PPZ)-dist(MMZ)) + 6*(dist(ZPZ)-dist(ZMZ)) + 3*(dist(MPZ)-dist(PMZ))
					+3*(dist(ZPP)-dist(ZMM)) - 3*(dist(ZMP)-dist(ZPM)) + 2*(dist(PPP)-dist(MMM))
					-2*(dist(PMP)-dist(MPM)) - 2*(dist(MMP)-dist(PPM)) + 2*(dist(MPP)-dist(PMM)) ) * scale;
	data->grad[2]=(  6*(dist(ZZP)-dist(ZZM)) + 3*(dist(ZPP)-dist(ZMM)) + 3*(dist(ZMP)-dist(ZPM))
					+3*(dist(PZP)-dist(MZM)) - 3*(dist(PZM)-dist(MZP)) + 2*(dist(PPP)-dist(MMM))
					+2*(dist(PMP)-dist(MPM)) + 2*(dist(MMP)-dist(PPM)) + 2*(dist(MPP)-dist(PMM)) ) * scale;

	float norm = sqrt(data->grad[0]*data->grad[0] + data->grad[1]*data->grad[1] + data->grad[2]*data->grad[2]);
	if (norm > 0.00001)
	{data->grad[0] /= norm; data->grad[1] /= norm; data->grad[2] /= norm;}
	data->DoneFlag = 10;	// indicates the gradient of this voxel has been successfully calculated.
}

inline void FastMarch::freeMem()
{
	POSITION pos = mapDist.GetStartPosition();
	int n = 0;
	FMContainer* data = NULL;
	while (pos)
	{
		data = NULL;
		mapDist.GetNextAssoc(pos, n, data);
		if (data != NULL)
			delete data;
	}
}

void FastMarch::ReinitPositiveHalf(BYTE* gridType) {
    heapSize = 0;
	closeSize = 0;
	InitPositive(gridType);
	InitHeap();
	March();
}

void FastMarch::ReinitNegativeHalf(BYTE* gridType)
{
    heapSize = 0;
	closeSize = 0;
	InitNegative(gridType);
	InitHeap();
	March();
}

void FastMarch::InitPositive(BYTE* gridType) {
    static int ci, pi;	// current index and precursor index
	static BYTE flag;	// flag of voxel with current index
	int ti = 0;

	// build the positive narrow band along z-direction
	FOR_YX
		ci = GI(i,j,0);
		for (k=1; k<NZ; k++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_TISSUE|VXL_MUCOSA) )
			{
				ti = (flag==VXL_TISSUE) ? ci : pi;
				AddClose(fz, ti);
			}
			ci = pi;
		}
	END_TWO

	// build the positive narrow band along y-direction
	FOR_ZX
		ci = GI(i,0,k);
		for (j=1; j<NY; j++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_TISSUE|VXL_MUCOSA) )
			{
				ti = (flag==VXL_TISSUE) ? ci : pi;
				AddClose(fy, ti);
			}
			ci = pi;
		}
	END_TWO


	// build the positive narrow band along x-direction
	FOR_ZY
		ci = GI(0,j,k);
		for (i=1; i<NX; i++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_TISSUE|VXL_MUCOSA) )
			{
				ti = (flag==VXL_TISSUE) ? ci : pi;
				AddClose(fx, ti);
			}
			ci = pi;
		}
	END_TWO
}

void FastMarch::InitNegative(BYTE* gridType) {
    static int ci, pi;	// current index and precursor index
	static BYTE flag;	// flag of voxel with current index
	int ti = 0;

	// build the negative narrow band along z-direction
	FOR_YX
		ci = GI(i,j,0);
		for (k=1; k<NZ; k++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag^gridType[pi])==VXL_TISSUE_BORDER && ( 
				(flag&gridType[pi])==VXL_MUCOSA/* ||
				(flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_AIR) */)
				)
			{
				ti = (flag==VXL_MUCOSA) ? ci : pi;
				AddClose(fz, ti);
			}
			ci = pi;
		}
	END_TWO

	// build the negative narrow band along y-direction
	FOR_ZX
		ci = GI(i,0,k);
		for (j=1; j<NY; j++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag^gridType[pi])==VXL_TISSUE_BORDER && ( 
				(flag&gridType[pi])==VXL_MUCOSA/* ||
				(flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_AIR) */)
				)
			{
				ti = (flag==VXL_MUCOSA) ? ci : pi;
				AddClose(fy, ti);
			}
			ci = pi;
		}
	END_TWO


	// build the negative narrow band along x-direction
	FOR_ZY
		ci = GI(0,j,k);
		for (i=1; i<NX; i++)
		{
			flag = gridType[ci];
			pi = GI(i,j,k);
			if ( (flag^gridType[pi])==VXL_TISSUE_BORDER && ( 
				(flag&gridType[pi])==VXL_MUCOSA/* ||
				(flag|gridType[pi])==(VXL_TISSUE_BORDER|VXL_AIR) */)
				)
			{
				ti = (flag==VXL_MUCOSA) ? ci : pi;
				AddClose(fx, ti);
			}
			ci = pi;
		}
	END_TWO
}

inline void FastMarch::AddClose(float f, int index)
{
	if (mapDist.Lookup(index, cd)){ 
		// the distance value has been previously calculated
		cd->value = min(cd->value, f);
	} else {
		// the distance value has not been previously calculated
		cd = new FMContainer;
		cd->DoneFlag = 0;
		cd->HeapPosition = -1;
		cd->value = f;
		mapDist.SetAt(index, cd);

		//Add index to list for initialization of close band
		ClosePoints[closeSize] = index;
		closeSize++;
	}
}

void FastMarch::InitHeap() {
    static int x,y,z;
	for(i = 0; i < closeSize; i++) {
		VERIFY(mapDist.Lookup(ClosePoints[i], cd));
		if (cd->HeapPosition == -1 && cd->DoneFlag == 0)
		{
			GIJK(ClosePoints[i], x, y, z);
			FindPhi(cd, ClosePoints[i], x, y, z);
		}
	}
}

void FastMarch::FindPhi(FMContainer* data, int index, int x, int y, int z) {
	static float phiX, phiY, phiZ, b, quotient, phi;
    static int a;
    static bool flagX, flagY, flagZ;

    phiX = phiY = phiZ = 0.;
    a = 0;
    flagX = flagY = flagZ = 0;

	//Find The phiS
	if (inDomain(x+1,y,z)) CheckFront (phiX, a, flagX, GI(x+1,  y,  z));
	if (inDomain(x-1,y,z)) CheckBehind(phiX, a, flagX, GI(x-1,  y,  z));
	if (inDomain(x,y+1,z)) CheckFront (phiY, a, flagY, GI(x  ,y+1,  z));
	if (inDomain(x,y-1,z)) CheckBehind(phiY, a, flagY, GI(x  ,y-1,  z));
    if (inDomain(x,y,z+1)) CheckFront (phiZ, a, flagZ, GI(x  ,y  ,z+1));
	if (inDomain(x,y,z-1)) CheckBehind(phiZ, a, flagZ, GI(x  ,y  ,z-1));

	// new value should be no less than that of the accepted neighbors
	if (a == 3)
	{
		if     ((phiX >= phiY) && (phiX >= phiZ))	CheckMax3(a, flagX, phiX, phiY, phiZ, fyyINV, fzzINV);
		else if((phiY >= phiX) && (phiY >= phiZ))	CheckMax3(a, flagY, phiY, phiX, phiZ, fxxINV, fzzINV);
		else										CheckMax3(a, flagZ, phiZ, phiX, phiY, fxxINV, fyyINV);
	}
	if(a == 2) {
		if(!flagX) {
			if(phiY >= phiZ)CheckMax2(a, flagY, phiY, phiZ, fzzINV);
			else			CheckMax2(a, flagZ, phiZ, phiY, fyyINV);
		}
		else if(!flagY){
			if(phiX >= phiZ)CheckMax2(a, flagX, phiX, phiZ, fxxINV);
			else			CheckMax2(a, flagZ, phiZ, phiX, fzzINV);
		}
		else {
			if(phiX >= phiY)CheckMax2(a, flagX, phiX, phiY, fxxINV);
			else			CheckMax2(a, flagY, phiY, phiX, fyyINV);
		}
	}
	if (a == 0) return; // ?? need further fixing...

	sq = sumINV;
	if (!flagX) sq -= fxxINV;
	if (!flagY) sq -= fyyINV;
	if (!flagZ) sq -= fzzINV;

	b = phiX*fxxINV + phiY*fyyINV + phiZ*fzzINV;
	quotient = square(b) - 
		sq * (square(phiX)*fxxINV + square(phiY)*fyyINV + square(phiZ)*fzzINV - 1);
	if(quotient < 0.) TRACE(" 0 \n");
	else {
		phi = b + sqrt(quotient);
		phi /= sq;
		if (data == NULL)
		{// allocate memory for new trial voxel
			data = new FMContainer;
			data->DoneFlag = 0;
			data->HeapPosition = -1;
			memset(data->grad, 0, 3*sizeof(float));
			mapDist.SetAt(index, data);
		}

		data->value = phi;
		if(data->HeapPosition==-1)	AddToHeap (data, index);
	    else						UpdateHeap(data, index); 
	}
}

inline void FastMarch::CheckFront(float& phi, int& a, bool& flag, int index) 
{
	if (mapDist.Lookup(index, cd) && cd->DoneFlag == 1) {
		phi = cd->value;
		flag = 1;
		a++;
	}
}

inline void FastMarch::CheckBehind(float& phi, int& a, bool& flag, int index)
{
	if(mapDist.Lookup(index, cd) && cd->DoneFlag == 1) {
		if(!flag) { phi = cd->value; a++; }
		else phi = min(cd->value, phi);
        flag = 1;
	}
}

inline void FastMarch::CheckMax2(int& a, bool &flag, float& phi1, const float &phi2, 
								 const float &inv2) {
	if(square((phi1-phi2)*inv2) > 1.) { phi1 = 0; a = 1; flag = 0; }
}

inline void FastMarch::CheckMax3(int& a, bool& flag, float& phi1, 
                                 const float &phi2, const float &phi3,
								 const float &inv2, const float &inv3) {
	if((square((phi1-phi2)*inv2) + square((phi1-phi3)*inv3)) > 1.)
    {   phi1 = 0; a = 2; flag = 0; }
}

void FastMarch::AddToHeap(FMContainer* data, int index) {
	static FMContainer* dy;
	FMHeap[heapSize] = index;
	data->HeapPosition = heapSize;
	int y, x = heapSize;
	for(x; x > 0; x = y) {
		y = (x-1)/2;
		VERIFY(mapDist.Lookup(FMHeap[y], dy));
		if (data->value < dy->value ){
			FMHeap[x] = FMHeap[y];
			dy->HeapPosition = x;
			FMHeap[y] = index;
			data->HeapPosition = y;
		}
		else
			break;
	}
	heapSize++;
}

void FastMarch::UpdateHeap(FMContainer* data, int index) {
	static FMContainer* dy;
	int y, x = data->HeapPosition;
	for(x; x > 0; x = y) {
		y = (x-1)/2;
		VERIFY(mapDist.Lookup(FMHeap[y], dy));
		if (data->value < dy->value ) {
			FMHeap[x] = FMHeap[y];
			dy->HeapPosition = x;
			FMHeap[y] = index;
			data->HeapPosition = y;
		}
		else
			break;
	}
}

void FastMarch::March() {
	static int x, y, z, index;
	FMContainer* data = NULL;
	for(index = PopHeap(data); index != -1; index = PopHeap(data)) {
		mapDist.Lookup(index, cd);
		if(cd->value > CurrenLimit) return;
		GIJK(index, x, y, z);

		index = GI(x-1,y,z);
		if (inDomain(x-1, y, z))
		{
			if (!mapDist.Lookup(index, cd)) FindPhi(NULL, index, x-1, y, z);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x-1, y, z);
		}

		index = GI(x+1,y,z);
		if (inDomain(x+1, y, z))
		{
			if (!mapDist.Lookup(index, cd))	FindPhi(NULL, index, x+1, y, z);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x+1, y, z);
		}

		index = GI(x,y-1,z);
		if (inDomain(x, y-1, z))
		{
			if (!mapDist.Lookup(index, cd))	FindPhi(NULL, index, x, y-1, z);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x, y-1, z);
		}

		index = GI(x,y+1,z);
		if (inDomain(x, y+1, z))
		{
			if (!mapDist.Lookup(index, cd))	FindPhi(NULL, index, x, y+1, z);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x, y+1, z);
		}
        
		index = GI(x,y,z-1);
		if (inDomain(x, y, z-1))
		{
			if (!mapDist.Lookup(index, cd))	FindPhi(NULL, index, x, y, z-1);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x, y, z-1);
		}

		index = GI(x,y,z+1);
		if (inDomain(x, y, z+1))
		{
			if (!mapDist.Lookup(index, cd))	FindPhi(NULL, index, x, y, z+1);
			else if (cd->DoneFlag == 0)		FindPhi(cd  , index, x, y, z+1);
		}
	}
}

int FastMarch::PopHeap(FMContainer*& data) {
	static FMContainer *d1, *di, *dlc, *drc;
	if(heapSize == 0)
		return -1;
	int index = FMHeap[0];
	VERIFY(mapDist.Lookup(index, data));
	data->DoneFlag = 1;
	heapSize--;
	FMHeap[0] = FMHeap[heapSize];
	VERIFY(mapDist.Lookup(FMHeap[heapSize], d1));
	d1->HeapPosition = 0;
	for(i = 0; i < (heapSize-1); i = j) {
		int lc = 2*i+1;
		int rc = 2*i+2;
		mapDist.Lookup(FMHeap[i], di);
		float current = di->value;
		float lv, rv;
		if(lc < heapSize) {
			mapDist.Lookup(FMHeap[lc], dlc);
			lv = dlc->value;
			di = dlc;
			if(rc < heapSize) {
				mapDist.Lookup(FMHeap[rc], drc);
				rv = drc->value;
				if(lv > rv) {
					lc = rc;
					lv = rv;
					di = drc;
				}
			}
			if(current > lv) {
				FMHeap[i] = FMHeap[lc];
				di->HeapPosition = i;
				FMHeap[lc] = FMHeap[heapSize];
				d1->HeapPosition = lc;
				j = lc;
			}
			else
				break;
		}
		else
			break;
	}
	return index;
}