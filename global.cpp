
#include "global.h"
using namespace std;
int i =0, j = 0, k = 0;
int NX = 512, NY = 512, NZ = 512, dj = 0, dk = 0, size = 0;
Float djInv = 1., dkInv = 1.;
float FX = 1., FY = 1., FZ = 1., RZ = 1.;

GlobalPara gPara;

void EndSwap(int *objp)
{
  unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
  std::reverse(memp, memp + sizeof(int));
}
void RAW_HEADER::load(CFile& file)
{
	file.Read(&nX, sizeof(int));
	file.Read(&nY, sizeof(int));
	file.Read(&nZ, sizeof(int));

	file.Read(&fX, sizeof(float));
	file.Read(&fY, sizeof(float));
	file.Read(&fZ, sizeof(float));
}

void SEG_HEADER::load(CFile& file)//--for test aim henry modified 2013-4-11
{   int xT = 0;
	file.Read(reinterpret_cast<char*>(&xT), sizeof(int));
    EndSwap( &xT );
	//nX = xT;

	int yT = 0;
	file.Read(reinterpret_cast<char*>(&yT), sizeof(int));
    EndSwap( &yT );
	//nY = yT;

	int zT = 0;
	file.Read(reinterpret_cast<char*>(&zT), sizeof(int));
    EndSwap( &zT );
	//nZ = zT;

	float ff = 0.0f;
	file.Read(reinterpret_cast<char*>(&ff) , sizeof(float));
	//endswap( &ff );
    f = ff;
}

void SEG_HEADER::write(CFile& file)
{
	file.Write(&nX, sizeof(int));
	file.Write(&nY, sizeof(int));
	file.Write(&nZ, sizeof(int));
	file.Write(&f , sizeof(float));
}



void GEO_HEADER::load(CFile& file)
{
	file.Read(&bOkGrad, sizeof(BOOL));
	file.Read(&bOkCurv, sizeof(BOOL));
	file.Read(&bOkLICC, sizeof(BOOL));
	file.Read(&bOkSICN, sizeof(BOOL));
	file.Read(&num	  , sizeof(int ));
	file.Read(&maxGrad, sizeof(Float));
	file.Read(&minGrad, sizeof(Float));
	file.Read(&maxCurv, sizeof(Float));
	file.Read(&minCurv, sizeof(Float));
	file.Read(&maxSI  , sizeof(Float));
	file.Read(&minSI  , sizeof(Float));
	file.Read(&maxCurvedness, sizeof(Float));
	file.Read(&minCurvedness, sizeof(Float));
}
void GEO_HEADER::write(CFile& file)
{
	file.Write(&bOkGrad, sizeof(BOOL));
	file.Write(&bOkCurv, sizeof(BOOL));
	file.Write(&bOkLICC, sizeof(BOOL));
	file.Write(&bOkSICN, sizeof(BOOL));
	file.Write(&num	   , sizeof(int ));
	file.Write(&maxGrad, sizeof(Float));
	file.Write(&minGrad, sizeof(Float));
	file.Write(&maxCurv, sizeof(Float));
	file.Write(&minCurv, sizeof(Float));
	file.Write(&maxSI  , sizeof(Float));
	file.Write(&minSI  , sizeof(Float));
	file.Write(&maxCurvedness, sizeof(Float));
	file.Write(&minCurvedness, sizeof(Float));
}

//--henry added new function for test 20121102
void GEO_HEADER::write_new(ofstream &out )
{
	out.write((char*)&bOkGrad, sizeof(BOOL));
	out.write((char*)&bOkCurv, sizeof(BOOL));
	out.write((char*)&bOkSICN, sizeof(BOOL));
	out.write((char*)&num	   , sizeof(int ));
	out.write((char*)&maxCurv, sizeof(Float));
	out.write((char*)&minCurv, sizeof(Float));
	out.write((char*)&maxSI  , sizeof(Float));
	out.write((char*)&minSI  , sizeof(Float));
	out.write((char*)&maxCurvedness, sizeof(Float));
	out.write((char*)&minCurvedness, sizeof(Float));
}

void PLP_HEADER::load(CFile& file)
{
	file.Read(&nIPC, sizeof(int));
	file.Read(&nTP, sizeof(int));
	file.Read(dtcIPC, MAX_NUM_TP*sizeof(int));
}

void PLP_HEADER::write(CFile& file)
{
	file.Write(&nIPC, sizeof(int));
	file.Write(&nTP, sizeof(int));
	file.Write(dtcIPC, MAX_NUM_TP*sizeof(int));
}


////////////////////////////////////////////////////////////////////////////////////
// class Candidate

int	Candidate::canID = 0;
CStdioFile Candidate::featFiles[FEAT_NUM+2];
bool Candidate::isFeatFilesOpen = false;
void Candidate::openFeatFiles(const CString& sPath)
{
	if (isFeatFilesOpen)
		return;

	CString fName;
	int n = 0;
	for (; n<FEAT_NUM+1; n++)
	{
		fName = sPath + "\\" + FEAT_NAME + FEAT_EXT[n];
		featFiles[n].Open(fName, CFile::typeText|CFile::modeCreate|CFile::modeNoTruncate|CFile::modeWrite);
		featFiles[n].SeekToEnd();
	}
	fName = sPath + "\\" + FEAT_NAME + ".cid";
	featFiles[n].Open(fName, CFile::typeText|CFile::modeCreate|CFile::modeNoTruncate|CFile::modeWrite);
	featFiles[n].SeekToEnd();
	isFeatFilesOpen = true;
}
void Candidate::closeFeatFiles()
{
	int n = 0;
	for (n=0; n<FEAT_NUM+2; n++)
		featFiles[n].Close();
	isFeatFilesOpen = false;
}

void Candidate::writeFeatursMat(Candidate* c)
{
	int n = 0;
	for (n=0; n<FEAT_NUM; n++)
	{// each feature
		CString s; s.Format("%f\n", c->feature(FEAT_TYPE(n)));
		featFiles[n].SeekToEnd(); featFiles[n].WriteString(s);
	}
	
	// update "features.tf"
	CString s = c->isTruth() ? "1\n" : "0\n";
	featFiles[n].SeekToEnd(); featFiles[n].WriteString(s);

	// update "features.cid"
	s.Format("%d\n", c->cid());
	featFiles[n+1].SeekToEnd(); featFiles[n+1].WriteString(s);
}


float Candidate::volEstimation(const IndexMap& mpMV, float* mv)
{
	GEO_TYPE tCV = CV_SMOOTH;
	CArray<float, float> cv; int im;
	vector<int>::iterator it;
	for (it=vxls.begin(); it<vxls.end(); it++)
	{
		mpMV.Lookup(*it, im); im *= GEO_NUM;
		cv.Add(mv[im+tCV]);
	}
	float cvm, cvv, cvs, cvk;
	calMeanVarSkewKurt(cv.GetData(), cv.GetCount(), cvm, cvv, cvs, cvk);

	// maximum radius
	float rad = 1. / (cvm-sqrt(cvv));
	return 4.*M_PI*pow(rad,3)/3.;
}

void Candidate::setVxlCanID(const IndexMap& mpMV, float* MV)
{
	int im = 0;
	vector<int>::iterator it;
	for (it=vxls.begin(); it<vxls.end(); it++)
	{
		if (mpMV.Lookup(*it, im))
		{
			im *= GEO_NUM;
			MV[im + CAN_ID] = ID;
		}
	}
}

void Candidate::setupLocalFrame()
{
	float vx[3];	// x-axis of the local frame
	float vy[3];	// y-axis of the local frame
					// noml is the z-axis of the local frame
	float org[3];	// the origin of the local frame
	
	// determine the origin
	int x, y, z; GIJK(center, x, y, z);
	org[0] = x, org[1] = y, org[2] = z;
	
	// determine the non-optimized x-axis
	float fx, fy, fz;
	fx = fabs(noml[0]);	fy = fabs(noml[1]);	fz = fabs(noml[2]);
	if (fx<fy && fx<fz)
	{
		vx[0] = 0;	vx[1] = noml[2];	vx[2] = -noml[1];
	}
	else if (fy<fx && fy<fz)
	{
		vx[0] = noml[2];	vx[1] = 0;	vx[2] = -noml[0];
	}
	else
	{
		vx[0] = noml[1];	vx[1] = -noml[0];	vx[2] = 0;
	}

	// determine the y-axis
	vtrCrossProduct(noml, vx, vy);
	// normalize vx and vy
	vtrUnitize(vx); vtrUnitize(vy);
	// setup the transient matrix from original frame to the non-optimized frame
	mtx.SetData(vx, vy, noml, org);
	mtx.GetInverse(mtx);

	// set up the matrix used to transform vectors
	Point3d pt(0,0,0);
	mV.SetData(vx, vy, noml, pt.getData());
	mV.GetInverse(mV);

	// setup the bounding boxes under the original and non-optimized frames
	cIterator it;
	float v[3];
	lfBox[0] = lfBox[1] = lfBox[2] = F_MAX;
	lfBox[3] = lfBox[4] = lfBox[5] = -F_MAX;
	ofBox[0] = ofBox[1] = ofBox[2] = I_P_INFINITE;
	ofBox[3] = ofBox[4] = ofBox[5] = I_N_INFINITE;
	for (it=begin(); it<end(); it++)
	{
		GIJK(*it, x, y, z); 

		if (ofBox[0] > x) ofBox[0] = x;
		if (ofBox[1] > y) ofBox[1] = y;
		if (ofBox[2] > z) ofBox[2] = z;

		if (ofBox[3] < x) ofBox[3] = x;
		if (ofBox[4] < y) ofBox[4] = y;
		if (ofBox[5] < z) ofBox[5] = z;

		v[0] = x, v[1] = y, v[2] = z;
		mtx.ApplyPoint(v);
		
		if (lfBox[0] > v[0]) lfBox[0] = v[0];
		if (lfBox[1] > v[1]) lfBox[1] = v[1];
		if (lfBox[2] > v[2]) lfBox[2] = v[2];

		if (lfBox[3] < v[0]) lfBox[3] = v[0];
		if (lfBox[4] < v[1]) lfBox[4] = v[1];
		if (lfBox[5] < v[2]) lfBox[5] = v[2];
	}

	// determine the optimized x,y-directions using PCA
	int dx = 512, id, i = 0;
	IndexMap iMap;
	for (it=begin(); it<end(); it++)
	{
		GIJK(*it, x, y, z); 
		v[0] = x, v[1] = y, v[2] = z;
		mtx.ApplyPoint(v);
		x = Round(v[0]-lfBox[0]), y = Round(v[1]-lfBox[1]);
		id = x + dx*y;
		if (!iMap.Lookup(id, id))
			iMap.SetAt(id, id);
	}
	int ptNum = iMap.GetCount();
	Float *fxx = new Float[ptNum];
	Float *fyy = new Float[ptNum];
	POSITION pos = iMap.GetStartPosition();
	while (pos != NULL)
	{
		iMap.GetNextAssoc(pos, id, id);
		fxx[i] = id%dx;	fyy[i] = id/dx;	i ++;
	}
	Float ex[3], ey[3];
	principalDir2D(fxx, fyy, ptNum, ex, ey);
	delete [] fxx, delete [] fyy;

	vx[0] = ex[0], vx[1] = ex[1], vx[2] = 0; // optimized x-direction in the non-optimized frame
	vy[0] = ey[0], vy[1] = ey[1], vy[2] = 0; // optimized y-direction in the non-optimized frame
	mV.InversePoint(vx);	// optimized x-direction in the original frame
	mV.InversePoint(vy);	// optimized y-direction in the original frame
	vtrUnitize(vx); vtrUnitize(vy);
	mtx.SetData(vx, vy, noml, org); 
	mtx.GetInverse(mtx);	// the transient matrix from the original frame to the optimized frame

	// re-set the matrix to transform vectors
	mV.SetData(vx, vy, noml, pt.getData());
	mV.GetInverse(mV);

	// re-calculate the local frame box
	lfBox[0] = lfBox[1] = lfBox[2] = F_MAX;
	lfBox[3] = lfBox[4] = lfBox[5] = -F_MAX;
	for (it=begin(); it<end(); it++)
	{
		GIJK(*it, x, y, z); 
		v[0] = x, v[1] = y, v[2] = z;
		mtx.ApplyPoint(v);

		if (lfBox[0] > v[0]) lfBox[0] = v[0];
		if (lfBox[1] > v[1]) lfBox[1] = v[1];
		if (lfBox[2] > v[2]) lfBox[2] = v[2];

		if (lfBox[3] < v[0]) lfBox[3] = v[0];
		if (lfBox[4] < v[1]) lfBox[4] = v[1];
		if (lfBox[5] < v[2]) lfBox[5] = v[2];
	}

	// axis-ratio
	Float w = lfBox[3]-lfBox[0], h = lfBox[4]-lfBox[1], l = lfBox[5]-lfBox[2];
	feat[AXIS_RATIO] = max(max(w,h),l) / min(min(w,h),l);
}

bool Candidate::isInShowBox(int x, int y, int z)
{
	static float xyR = 2, zR = 1.5;
	float dx = xyR*(lfBox[3]-lfBox[0]), dy = xyR*(lfBox[4]-lfBox[1]), dz = zR*(lfBox[5]-lfBox[2]);
	float showBox[6] = { lfBox[0]-dx, lfBox[1]-dy, lfBox[2]-dz, 
						 lfBox[3]+dx, lfBox[4]+dy, lfBox[5]+dz };
	float v[3]; v[0] = x, v[1] = y, v[2] = z;
	mtx.ApplyPoint(v);
	if (v[0]<showBox[3] && v[0]>showBox[0] &&
		v[1]<showBox[4] && v[1]>showBox[1] &&
		v[2]<showBox[5] && v[2]>showBox[2] )
		return true;
	else return false;
}

//////////////////////////////////////////////////////////////////////////
//
#define InCan(i,j,k)	(vxlCan.Lookup(ig=GI((i),(j),(k)), ig))		// in reliable region
#define InCan0(i,j,k)	(mpCanVxl.Lookup(ig=GI((i),(j),(k)), ig))	// in whole candidate

static void h2c(Float h, Float& r, Float& g, Float& b, Float& a)
{// mapping Hounsfield value to color, according to Perry's design in 2004.
	h = max(-1000, h); h = min(h, 1000);

	if (h < 200)
	{
		// r-mapping
//		if (h<-64) r = 0; else if (h<-16) r = 0.68-0.03*(-16.f-h)/48.f; else r = 1; // according to the fig.2 in Pickhardt's paper
		if (h<-64) r = 0; else if (h<-16) r = 0.68-0.03*(-16.f-h)/48.f; else r = 1;

		// g-mapping
//		if (h<-436)g = 0; else if (h<-16) g = 0.73-0.23*(-16.-h)/420.f; else g = 0;
		if (h<-436)g = 0; else if (h<-16) g = 0.73-0.23*(-16.-h)/420.f; else g = 0;

		// b-mapping
//		if (h<-920)b = 0; else if (h<-128)b = 0.95-0.40*(-128-h)/792.f; else b = 0;
		if (h<-920)b = 0; else if (h<-128)b = 0.95-0.40*(-128-h)/792.f; else b = 0;
	} else 
	{// the white component in Perry's chart.
		r = g = b = 1.;
	}

	// a-mapping
	if (h<-920)a = 0; else if (h<-16) a = 0.15-0.15*(-16.-h)/904.f; else a = 0.3+0.05*(h+16)/1016.;
}

/*
static void blend(Float& dr, Float& dg, Float& db, Float& da,
				  Float  sr, Float  sg, Float  sb, Float  sa)
{// simulate the volume rendering using glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
	dr = min(1, sa*sr+dr*(1-sa));
	dg = min(1, sa*sg+dg*(1-sa));
	db = min(1, sa*sb+db*(1-sa));
	da = min(1, sa*sa+da*(1-sa));
}
*/
static void blend(Float& dr, Float& dg, Float& db, Float& da,
				  Float  sr, Float  sg, Float  sb, Float  sa)
{// simulate the volume rendering using glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
	dr = sa*sr+dr*(1-sa);
	dg = sa*sg+dg*(1-sa);
	db = sa*sb+db*(1-sa);
	da = sa*sa+da*(1-sa);
}

//////////////////////////////////////////////////////////////////////////
// Generate the image patterns
//////////////////////////////////////////////////////////////////////////
void Candidate::genBiopsyImgs(RAW_DATA_TYPE* cleanRaw, RAW_DATA_TYPE* dirtyRaw, 
							  Grid2dF& gA, Grid2dF& gC, Grid2dF& gS, 
							  Grid2d4v& cA, Grid2d4v& cC, Grid2d4v& cS,
							  Grid2dU& tA, Grid2dU& tC, Grid2dU& tS)
{
	// The min CT value in objective volume of interest (this candidate, or the extended box)
	float FOUT = -1000; // F_MAX;

	float dx, dy, dz, vl[3], vg[3], g[8], gd[8], v, vd;
	int nx, ny, nz, x, y, z, ig;

	static float ext = 1.;	// extend 1 mm
	dx = ext/FX, dy = ext/FY, dz = ext/FZ;
	float box[6] = { lfBox[0]-dx, lfBox[1]-dy, lfBox[2]-dz, 
					 lfBox[3]+dx, lfBox[4]+dy, lfBox[5]+dz };
	dx = dy = dz = 0.25;	// sampling interval to generate the images

	nx = 1+int( (box[3]-box[0])/dx );
	ny = 1+int( (box[4]-box[1])/dy );
	nz = 1+int( (box[5]-box[2])/dz );

	// allocate memory
	gA.init(nx, ny); gC.init(ny, nz); gS.init(nz, nx);
	cA.init(nx, ny); cC.init(ny, nz); cS.init(nz, nx);
	tA.init(nx, ny); tC.init(ny, nz); tS.init(nz, nx);

	IndexMap vxlCan;
	for (ig=0; ig<vxls.size(); ig++)
		vxlCan.SetAt(vxls[ig], ig);

	Float sr, sg, sb, sa;
	for (z=0; z<nz; z++) { for (y=0; y<ny; y++) { for (x=0; x<nx; x++) {
		vl[0] = box[0] + x*dx; vl[1] = box[1] + y*dy; vl[2] = box[2] + z*dz;
		mtx.InversePoint(vl, vg);
		i = int(vg[0]); j = int(vg[1]); k = int(vg[2]);

		if (inDomain(i  ,j  ,k  )) {g[0] = cleanRaw[GI(i  ,j  ,k  )]; gd[0] = dirtyRaw[GI(i  ,j  ,k  )];} else {g[0] = FOUT; gd[0] = FOUT;}
		if (inDomain(i+1,j  ,k  )) {g[1] = cleanRaw[GI(i+1,j  ,k  )]; gd[1] = dirtyRaw[GI(i+1,j  ,k  )];} else {g[1] = FOUT; gd[1] = FOUT;}
		if (inDomain(i+1,j+1,k  )) {g[2] = cleanRaw[GI(i+1,j+1,k  )]; gd[2] = dirtyRaw[GI(i+1,j+1,k  )];} else {g[2] = FOUT; gd[2] = FOUT;}
		if (inDomain(i  ,j+1,k  )) {g[3] = cleanRaw[GI(i  ,j+1,k  )]; gd[3] = dirtyRaw[GI(i  ,j+1,k  )];} else {g[3] = FOUT; gd[3] = FOUT;}
		if (inDomain(i  ,j  ,k+1)) {g[4] = cleanRaw[GI(i  ,j  ,k+1)]; gd[4] = dirtyRaw[GI(i  ,j  ,k+1)];} else {g[4] = FOUT; gd[4] = FOUT;}
		if (inDomain(i+1,j  ,k+1)) {g[5] = cleanRaw[GI(i+1,j  ,k+1)]; gd[5] = dirtyRaw[GI(i+1,j  ,k+1)];} else {g[5] = FOUT; gd[5] = FOUT;}
		if (inDomain(i+1,j+1,k+1)) {g[6] = cleanRaw[GI(i+1,j+1,k+1)]; gd[6] = dirtyRaw[GI(i+1,j+1,k+1)];} else {g[6] = FOUT; gd[6] = FOUT;}
		if (inDomain(i  ,j+1,k+1)) {g[7] = cleanRaw[GI(i  ,j+1,k+1)]; gd[7] = dirtyRaw[GI(i  ,j+1,k+1)];} else {g[7] = FOUT; gd[7] = FOUT;}

		v  = trilerp( g[0], g[1], g[3], g[2],  g[4],  g[5],  g[7],  g[6], vg[0]-i, vg[1]-j, vg[2]-k);
		vd = trilerp(gd[0],gd[1],gd[3],gd[2], gd[4], gd[5], gd[7], gd[6], vg[0]-i, vg[1]-j, vg[2]-k);

		// Simple linear accumulation to generate gray image
		gA[nx*y + x] += v; gC[ny*z + y] += v; gS[nz*x + z] += v;

		// mapping to a colorful image ... 
		h2c(vd, sr, sg, sb, sa);
		blend(cA(x,y,0), cA(x,y,1), cA(x,y,2), cA(x,y,3), sr, sg, sb, sa);
		blend(cC(y,z,0), cC(y,z,1), cC(y,z,2), cC(y,z,3), sr, sg, sb, sa);
		blend(cS(z,x,0), cS(z,x,1), cS(z,x,2), cS(z,x,3), sr, sg, sb, sa);

		// determine pixel type
		i = clamp(Round(vg[0]), 0, NX-1);
		j = clamp(Round(vg[1]), 0, NY-1);
		k = clamp(Round(vg[2]), 0, NZ-1);
		if (InCan(i, j, k))
		{// flag plp pixels
			tA(x,y) |= type[ig]; tC(y,z) |= type[ig]; tS(z,x) |= type[ig];
		}
	}}}
	for (z=0; z<nz; z++) { for (y=0; y<ny; y++){
		if (tC(y,z)&(PLP_SEED|PLP_GROW|PLP_EXT)) ;
		else tC(y,z) = PLP_OUT;
	}}
	for (z=0; z<nz; z++) { for (x=0; x<nx; x++){
		if (tS(z,x)&(PLP_SEED|PLP_GROW|PLP_EXT)) ;
		else tS(z,x) = PLP_OUT;
	}}
	for (y=0; y<ny; y++) { for (x=0; x<nx; x++){
		if (tA(x,y)&(PLP_SEED|PLP_GROW|PLP_EXT)) ;
		else tA(x,y) = PLP_OUT;
	}}
	gA.rescale(); gC.rescale(); gS.rescale();

	//
	int i = 0, num = cA.getSize(nx, ny);
	float max = -F_MAX;
	Point3d e;
	for (i=0; i<num; i++)
	{
		cA.getElm3d(i, e);
		if (max < e[0]) max = e[0];
		if (max < e[1]) max = e[1];
		if (max < e[2]) max = e[2];
	} max = 1./max;
	for (i=0; i<num; i++)
	{
		cA.getElm3d(i, e);
		e *= max;
		cA.setElm3d(i, e);
	}

	num = cC.getSize(nx, ny);
	max = -F_MAX;
	for (i=0; i<num; i++)
	{
		cC.getElm3d(i, e);
		if (max < e[0]) max = e[0];
		if (max < e[1]) max = e[1];
		if (max < e[2]) max = e[2];
	} max = 1./max;
	for (i=0; i<num; i++)
	{
		cC.getElm3d(i, e);
		e *= max;
		cC.setElm3d(i, e);
	}

	num = cS.getSize(nx, ny);
	max = -F_MAX;
	for (i=0; i<num; i++)
	{
		cS.getElm3d(i, e);
		if (max < e[0]) max = e[0];
		if (max < e[1]) max = e[1];
		if (max < e[2]) max = e[2];
	} max = 1./max;
	for (i=0; i<num; i++)
	{
		cS.getElm3d(i, e);
		e *= max;
		cS.setElm3d(i, e);
	}
}

void Candidate::genBiopsyFeats(RAW_DATA_TYPE* cleanRaw, RAW_DATA_TYPE* dirtyRaw, const CString& path)
{
	GIJK(center,i,j,k);

	Grid2dF			*gA = new Grid2dF,	// Axial gray biopsy image
					*gC = new Grid2dF,	// Coronal gray biopsy image
					*gS = new Grid2dF;	// Saggital gray biopsy image
	Grid2d4v		*cA = new Grid2d4v,	// Axial colorful biopsy image
					*cC = new Grid2d4v,	// Coronal colorful biopsy image
					*cS = new Grid2d4v;	// Saggital colorful biopsy image
	Grid2dU			*tA = new Grid2dU,	// pixel types of axial gray/colorful biopsy image
					*tC = new Grid2dU,	// pixel types of coronal gray/colorful biopsy image
					*tS = new Grid2dU;	// pixel types of saggital gray/colorful biopsy image
	
	// generate gray and colorful images
	genBiopsyImgs(cleanRaw, dirtyRaw, *gA, *gC, *gS, *cA, *cC, *cS, *tA, *tC, *tS);

	// process gray images
	proGrayImg(*gA, *gC, *gS, *tA, *tC, *tS);

	// release memory of gray images
	delete gA; delete gC; delete gS;
	// clear the flags
	tA->Clear(GRY_GROW|GRY_BGD);
	tC->Clear(GRY_SEED|GRY_GROW|GRY_BGD);
	tS->Clear(GRY_SEED|GRY_GROW|GRY_BGD);

 	// process colorful images
	proClrImg(*cA, *cC, *cS, *tA, *tC, *tS);

 	// release memory
 	delete cA; delete cC; delete cS;
 	delete tA; delete tC; delete tS;
}


//////////////////////////////////////////////////////////////////////////
// new version built during drafting the paper
void Candidate::proGrayImg(Grid2dF& gA, Grid2dF& gC, Grid2dF& gS,
					Grid2dU& tA, Grid2dU& tC, Grid2dU& tS)
{
	int nx, ny, nz;	gA.GetSize(nx, ny);	nz = gC.GetNy();

	genCCASeeds(gA, tA, PLP_SEED);
	genCCASeeds(gC, tC, PLP_SEED|PLP_GROW);
	genCCASeeds(gS, tS, PLP_SEED|PLP_GROW);

	// extracting highlighted patch in axial image, light patches in saggital and coronal images
	gA.connCompFromSeeds(0.75, 1.15, tA, GRY_SEED, GRY_GROW, true); tA.removeHoles(GRY_SEED, GRY_GROW);
	gC.connCompFromSeeds(0.95, 1.45, tC, GRY_SEED, GRY_GROW, true); tC.removeHoles(GRY_SEED, GRY_GROW);
	gS.connCompFromSeeds(0.95, 1.45, tS, GRY_SEED, GRY_GROW, true); tS.removeHoles(GRY_SEED, GRY_GROW);
	genCCASeeds(gA, tA, GRY_GROW); tA.Clear(GRY_GROW);
	genCCASeeds(gC, tC, GRY_GROW); tC.Clear(GRY_GROW);
	genCCASeeds(gS, tS, GRY_GROW); tS.Clear(GRY_GROW);

	// extracting bright area in saggital and coronal patterns
	gC.connCompFromSeeds(0.95, 2.85, tC, GRY_SEED, GRY_GROW);
	gS.connCompFromSeeds(0.95, 2.85, tS, GRY_SEED, GRY_GROW);

	genGrayFeats(tA, tC, tS, gA, gC, gS);
}

void Candidate::genGrayFeats(Grid2dU& tA, Grid2dU& tC, Grid2dU& tS,
							 Grid2dF& gA, Grid2dF& gC, Grid2dF& gS)
{
	unsigned char flag = 0x80;
	int num = 0, nx, ny, i, x, y, n;

	// 
	unsigned char HA = GRY_SEED,	// highlighted patch in axial image
		LCS = GRY_SEED,				// light patch in coronal and sagittal images
		BCS = GRY_GROW;				// bright patch in coronal and sagittal images

	// handle axial image
	num = tA.GetSize(nx, ny);

	vector<int> cl, nl, *pcl, *pnl, *temp;
	Float mh = 0.f, mb = 0.f, mx = 0, my = 0;
	int nh = 0;
	for (i=0; i<num; i++)
	{//
		if (tA[i] & HA)
		{
			nh ++;	mh += gA[i];	tA[i] |= flag;
			x = i%nx; y = i/nx; mx += x; my += y;
			cl.push_back(i);
		}
	} mh /= nh; mx /= nh; my /= nh;
	
	pcl = &cl, pnl = &nl;
	int LS = 4, step = 0, nb = 0;
	while (step < LS)
	{
		pnl->clear();
		for (i=0; i<pcl->size(); i++)
		{
			n = (*pcl)[i];	x = n%nx;	y = n/nx;
			if (x-1>=0) { if (!(tA[n- 1]&flag)) {mb += gA[n- 1]; tA[n- 1]|=flag; pnl->push_back(n- 1);} }
			if (x+1<nx) { if (!(tA[n+ 1]&flag)) {mb += gA[n+ 1]; tA[n+ 1]|=flag; pnl->push_back(n+ 1);} }
			if (y-1>=0) { if (!(tA[n-nx]&flag)) {mb += gA[n-nx]; tA[n-nx]|=flag; pnl->push_back(n-nx);} }
			if (y+1<ny) { if (!(tA[n+nx]&flag)) {mb += gA[n+nx]; tA[n+nx]|=flag; pnl->push_back(n+nx);} }
		}
		nb += pnl->size();
		temp = pcl; pcl = pnl; pnl = temp;
		step ++;
	} mb /= nb;
	feat[AXI_HR] = mh/mb; // highlighting ratio

	int nBin = 60;
	Float delta = 360.0f/nBin, ang;
	Grid1dI Bin(nBin);
	for (y=0; y<ny; y++){ for (x=0; x<nx; x++){
		if (tA(x,y) & HA)
		{
			ang = RAD2DEG(angle(x-mx, y-my));
			Bin[int(ang/delta)] ++;
		}
	}}
	Float mn;
	Bin.meanvar(mn, feat[AXI_DL]); // disk-likeness

	// handle coronal image
	num = tC.GetSize(nx, ny);

	Float mlight = 0, mbright = 0, cxl = 0, cyl = 0, cxb = 0, cyb = 0;
	int nlight = 0, nbright = 0;

	for (y=0; y<ny; y++){ for (x=0; x<nx; x++){
		if (tC(x, y) & LCS)
		{
			nlight ++;		mlight += gC(x, y);
			cxl += x;		cyl += y;
		}
		else if (tC(x, y) & BCS)
		{
			nbright ++;		mbright += gC(x, y);
			cxb += x;		cyb += y;
		}
	}}
	feat[COR_LR] = (mlight/nlight) / (mbright/nbright); // lightness ratio
	feat[COR_GCX] = cxl/(nlight*nx);	// normalized cx of light patch
	feat[COR_GCY] = cyl/(nlight*ny);	// normalized cy of light patch
	feat[COR_BCX] = cxb/(nbright*nx);	// normalized cx of bright patch
	feat[COR_BCY] = cyb/(nbright*ny);	// normalized cy of bright patch
	if (nbright == 0)
	{// assume they are FPs
		feat[COR_LR] = 1;
		feat[COR_BCX] = 1;
		feat[COR_BCY] = 1;
	}

	// handle sagittal image
	num = tS.GetSize(nx, ny);

	mlight = 0, mbright = 0, cxl = 0, cyl = 0, cxb = 0, cyb = 0;
	nlight = 0, nbright = 0;

	for (y=0; y<ny; y++){ for (x=0; x<nx; x++){
		if (tS(x, y) & LCS)
		{
			nlight ++;		mlight += gS(x, y);
			cxl += x;		cyl += y;
		}
		else if (tS(x, y) & BCS)
		{
			nbright ++;		mbright += gS(x, y);
			cxb += x;		cyb += y;
		}
	}}
	feat[SAG_LR] = (mlight/nlight) / (mbright/nbright); // lightness ratio
	feat[SAG_GCX] = cxl/(nlight*nx);	// normalized cx of light patch
	feat[SAG_GCY] = cyl/(nlight*ny);	// normalized cy of light patch
	feat[SAG_BCX] = cxb/(nbright*nx);	// normalized cx of bright patch
	feat[SAG_BCY] = cyb/(nbright*ny);	// normalized cy of bright patch
	if (nbright == 0)
	{// assume they are FPs
		feat[SAG_LR] = 1;
		feat[SAG_BCX] = 1;
		feat[SAG_BCY] = 1;
	}
}

void Candidate::proClrImg(Grid2d4v& cA, Grid2d4v& cC, Grid2d4v& cS,
				   Grid2dU& tA, Grid2dU& tC, Grid2dU& tS)
{
	cA.clamp(); cC.clamp(); cS.clamp();

	// output colorful images
	genClrFeats(cA, tA, CLR_ROI);
}

//////////////////////////////////////////////////////////////////////////
bool Candidate::RED(Float H, Float S, Float I)
{
	if ((H<1./6.&&H >= 0) || (H>5./6. && H<=1))
	{
		if (S > 0.23)
		{
			if (I>0.16 && I<0.78)
				return true;
		}
	}
	return false;
}
bool Candidate::WHITE(Float H, Float S, Float I)
{
	if (S < 0.23) return true;
	else if (I > 0.80) return true;
	else return false;
}
bool Candidate::GREEN(Float H, Float S, Float I)
{
	return true;
}

//////////////////////////////////////////////////////////////////////////
//
void Candidate::genClrFeats(Grid2d4v& c, Grid2dU& t, unsigned char or)
{// TEMP_12 ->
	Color3d cfr(1,0,0), cfw(1,1,1), e;
	int i, x, y, n, nx, ny;
	n = c.getSize(nx, ny);
	float fr, fw, fd, sum_fd = 0, sum_fr = 0, sum_fw = 0, cx = .5*nx, cy = .5*ny;
	for (y=0; y<ny; y++) { for (x=0; x<nx; x++) {
		i = x + nx*y;
		if (t[i] & or)
		{
			c.getElm3d(i, e);
			fr = 1./(e.similarity(cfr)+0.0001);
			fw = 1./(e.similarity(cfw)+0.0001);
			fd = 1./(sqrt((x-cx)*(x-cx) + (y-cy)*(y-cy))+.01);
			sum_fr += fr * fd;
			sum_fw += fw * fd;
			sum_fd += fd;
		}
	}}
	feat[TEMP_00] = sum_fr / sum_fd;
	feat[TEMP_01] = sum_fw / sum_fd;

	// HSI color model
	Float sigma = 0;	// the ...
	cx = cy = 0;		// the center of the PLP_SEED patch
	n = 0;
	for (y=0; y<ny; y++) { for (x=0; x<nx; x++) {
		if (t[x + nx*y] & PLP_SEED)
		{
			n ++;
			cx += x; cy += y;
		}
	}} cx /= n; cy /= n;
	Float maxx = -1, minx = nx+1, maxy = -1, miny = ny+1;
	for (y=0; y<ny; y++) { for (x=0; x<nx; x++) {
		if (t[x + nx*y] & or)
		{
			if (maxx < x) maxx = x;
			if (minx > x) minx = x;
			if (maxy < y) maxy = y;
			if (miny > y) miny = y;
		}
	}}
	sigma = 0.618*min(maxx-minx, maxy-miny);

	Float H, S, I, w, sum_H = 0, sum_S = 0, sum_I = 0, sum_w = 0;
	for (y=0; y<ny; y++) { for (x=0; x<nx; x++) {
		i = x + nx*y;
		if (t[i] & or)
		{// (cx, cy) should take the center of the PLP_SEED patch, but what about the sigma?
			c.getElm3d(i, e);
			RGB2HSI(e[0], e[1], e[2], H, S, I);
			w = exp(-((x-cx)*(x-cx) + (y-cy)*(y-cy))/(2*sigma*sigma)) / (sqrt(2*PI)*sigma);
			sum_H += w*H;
			sum_S += w*S;
			sum_I += w*I;
			sum_w += w;
		}
	}}
	feat[CLR_AH] = sum_H / sum_w;
	feat[CLR_AS] = sum_S / sum_w;
	feat[CLR_AI] = sum_I / sum_w;
}

//////////////////////////////////////////////////////////////////////////
void Candidate::genCCASeeds(Grid2dF& g, Grid2dU& t, unsigned char s)
{
	int i, j, n;
	n = g.GetSize(i, j);
	for (i=0; i<n; i++)
		if (t[i]&s && g[i]>F_ZERO)
			t[i] |= GRY_SEED;
}

//////////////////////////////////////////////////////////////////////////
//
void Candidate::save(CFile& file)
{
	file.Write(&ID, sizeof(int));
	file.Write(&center, sizeof(int));
	file.Write(&TPID, sizeof(int));
	int numFeat = FEAT_NUM;
	file.Write(&numFeat, sizeof(int));
	file.Write(feat, numFeat*sizeof(float));
	file.Write(&CID, sizeof(int));

	// original/local frames
	file.Write(noml, 3*sizeof(float));
	file.Write(ofBox, 6*sizeof(int));
	file.Write(lfBox, 6*sizeof(Float));
	mV.save(&file);
	mtx.save(&file);


	// 
	file.Write(&num, sizeof(int));
	int* n = new int[num]; int nn = 0;
	for (cIterator it=begin(); it<end(); it++) n[nn++] = *it;
	file.Write(n, num*sizeof(int));
	delete [] n;
	file.Write(&type[0], num*sizeof(BYTE));

	// plp voxels
	int nPV = (int)plpVxls.size();
	file.Write(&nPV, sizeof(int));
	n = new int[nPV]; nn = 0;
	for (cIterator it=plpVxls.begin(); it<plpVxls.end(); it++) n[nn++] = *it;
	file.Write(n, nPV*sizeof(int));
	delete [] n;
}

void Candidate::load(CFile& file)
{
	file.Read(&ID, sizeof(int));
	file.Read(&center, sizeof(int));
	file.Read(&TPID, sizeof(int));
	int numFeat = FEAT_NUM;
	file.Read(&numFeat, sizeof(int)); ASSERT(numFeat == FEAT_NUM);
	file.Read(feat, numFeat*sizeof(float));
	file.Read(&CID, sizeof(int));

	file.Read(noml, 3*sizeof(float));
	file.Read(ofBox, 6*sizeof(int));
	file.Read(lfBox, 6*sizeof(Float));
	mV.load(&file);
	mtx.load(&file);

	file.Read(&num, sizeof(int));
	int* n = new int[num];
	file.Read(n, num*sizeof(int));
	for (int nn=0; nn<num; nn++) vxls.push_back(n[nn]);
	delete [] n;
	type.resize(num);
	file.Read(&type[0], num*sizeof(BYTE));

	//
	int nPV = 0;
	file.Read(&nPV, sizeof(int));
	n = new int[nPV];
	file.Read(n, nPV*sizeof(int));
	for (int ii=0; ii<nPV; ii++) plpVxls.push_back(n[ii]);
	delete [] n;
}

void Candidate::saveCAD(CStdioFile& file, float zf)
{

/*	// Save all the voxels in the detection
	int x, y, z, id, nv; 
	IndexMap v;
	for (cIterator it=begin(); it<end(); it++)
	{
		GIJK(*it, x, y, z); z = int(z*zf);
		id = GI(x, y, z);
		if (!v.Lookup(id, nv))
		{
			v.SetAt(id, *it);
		}
	}

	GIJK(center, x, y, z);
	CString sMsg, temp;
	sMsg.Format("%d\t%d\t%d\t1\t%d\t", x, y, int(z*zf), v.GetSize());
	POSITION pos = v.GetStartPosition();
	while (pos != NULL) {
		v.GetNextAssoc(pos, id, nv); GIJK(nv, x, y, z);
		temp.Format("%d\t%d\t%d\t", x, y, int(z*zf));
		sMsg += temp;
	}
	sMsg.TrimRight(); sMsg += "\n";
	file.WriteString(sMsg);
*/

	// just save the center of the detection
	int x, y, z, id, nv; 
	GIJK(center, x, y, z);
	CString sMsg, temp;
	sMsg.Format("%d\t%d\t%d", x, y, int(z*zf));
	sMsg.TrimRight(); sMsg += "\n";
	file.WriteString(sMsg);
}


////////////////////////////////////////////////////////////////////////////////////
// Description:
//	Calculate the Gaussian kernels for computing the first, second derivatives for image data.
//	On Aug. 23, 2007, it is modified based on the code of Zigang Wang, which was based on paper
//	"1992_Using partial derivatives of 3D images to extract typical surface features". In the 
//	following code, the points with the f_values under LOW_LIMIT are cut off.
//
// Parameters:
//	alpha ------ the input alpha value.
//	number ----- how many points there are in the kernel.
//	table ------ returns the kernels for computing 0th, 1st, 2nd-derivatives.
////////////////////////////////////////////////////////////////////////////////////
void calGaussian(Float alpha, int number[3], Float table[3][100])
{
	static int Max_Width = 100;
//	static Float LOW_LIMIT = 0.01;// why 0.00005? it's an empirical value.

	Float expNA = exp(-alpha), expNA2 = exp(-2*alpha), expNA3 = exp(-3*alpha), expNA4 = exp(-4*alpha);
	Float c0 = (1 - expNA) * (1 - expNA) / (1 + 2*expNA*alpha - expNA2);
	Float c1 = -(1 - expNA) * (1 - expNA) * (1 - expNA) / (2*alpha*alpha*expNA*(1 + expNA));
	Float c2 = -2*(1 - expNA) * (1 - expNA) * (1 - expNA) * (1 - expNA) /
				(1 + 2*expNA - 2*expNA3 - expNA4);
	Float c3 = (1 - expNA2) / (2*alpha*expNA); 

	Float	t1;
	int		num1;

// for zero derivative
	i  = -Max_Width;
	t1 = c0 * (1 + alpha*abs(i)) * exp(-alpha*abs(i));
	// get rid of those points with f0 under gPara.GAUSSIAN_LOW_LIMIT
	while (fabs(t1) < gPara.GAUSSIAN_LOW_LIMIT && i++ <= 0)
     	t1 = c0 * (1 + alpha * abs(i)) * exp(-alpha * abs(i)); 

	if(i > 0) number[0] = 0;
	else {
     	num1 = 0;
//		for(k = -i; k >= i; k--)
		for(k = i; k <= -i; k++)
		{
			t1 = c0 * (1 + alpha * abs(k)) * exp(-alpha * abs(k));
			table[0][num1++] = t1;
		}
     	number[0] = num1;
	}

// for first derivative
	i = -Max_Width;
	t1 = -c1 * i * alpha * alpha * exp(-alpha * abs(i));
	// get rid of those points with f1 under gPara.GAUSSIAN_LOW_LIMIT
	while(fabs(t1) < gPara.GAUSSIAN_LOW_LIMIT && i++ <= 0)  
		t1 = -c1 * i * alpha * alpha * exp(-alpha * abs(i)); 
	if (i > 0) number[1] = 0;
	else {
     	num1 = 0;
//		for(k = -i; k >= i; k--)
		for(k = i; k <= -i; k++)
		{
	    	t1 = -c1 * k * alpha * alpha * exp(-alpha * abs(k)); 
			table[1][num1++] = t1;
		}
     	number[1] = num1;
	}

// for second derivative

	i = -Max_Width;
	t1 = c2 * (1 - c3 * alpha * abs(i)) * exp(-alpha * abs(i)); 
	// get rid of those points with f2 under gPara.GAUSSIAN_LOW_LIMIT
	while(fabs(t1) < gPara.GAUSSIAN_LOW_LIMIT && i++ <= 0)  
     	t1 = c2 * (1 - c3 * alpha * abs(i)) * exp(-alpha * abs(i)); 
	if (i > 0) number[2] = 0;
	else {
     	num1 = 0;
//		for(k = -i; k >= i; k--)
		for(k = i; k <= -i; k++)
		{
        	t1 = c2 * (1 - c3 * alpha * abs(k)) * exp(-alpha * abs(k)); 
			table[2][num1++] = t1;
		}
     	number[2] = num1;
	}
}


////////////////////////////////////////////////////////////////////////////////////
// Description:
//	Calculate the Gaussian kernels, which is returned in a 1D float array, for computing 
//	the first, second derivatives for image data.
//	On Sept. 19, 2007, it is modified based on the above code. In the following code, 
//	a fixed window size is used to limit the 0th, 1st and 2nd order derivative.
//
// Add normalizing for each kernel, Oct. 11, 2007.
//
// Parameters:
//	alpha ------ the input alpha value.
//	len -------- the length of the kernel, should be a prime number and len/2 is the radius of the window.
//	zero ------- returns the kernels in 1D array for the 0-order derivatives. If it is NULL, do not return.
//	fst, snd --- returns the kernels in 1D array for the 1st, 2nd-derivatives in x,y,z-directions respectively. If they are NULL, do not return.
////////////////////////////////////////////////////////////////////////////////////
void calGaussian(Float alpha, int len, Float* zero, Float* fst[3], Float* snd[6])
{
	Float* temp[3], t[3]={0};
	temp[0] = new Float[len]; // save the t-values of 0th-order derivatives
	temp[1] = new Float[len]; // save the t-values of 1st-order derivatives
	temp[2] = new Float[len]; // save the t-values of 2nd-order derivatives

	Float expNA = exp(-alpha), expNA2 = exp(-2*alpha), expNA3 = exp(-3*alpha), expNA4 = exp(-4*alpha);
	Float c0 = (1 - expNA) * (1 - expNA) / (1 + 2*expNA*alpha - expNA2);
	Float c1 = -(1 - expNA) * (1 - expNA) * (1 - expNA) / (2*alpha*alpha*expNA*(1 + expNA));
	Float c2 = -2*(1 - expNA) * (1 - expNA) * (1 - expNA) * (1 - expNA) /
				(1 + 2*expNA - 2*expNA3 - expNA4);
	Float c3 = (1 - expNA2) / (2*alpha*expNA); 

	int rad = len/2, ii, jj, kk;
	for (ii=-rad; ii<=rad; ii++)
	{
		temp[0][ii+rad] =  c0 * (1 + alpha * abs(ii)) * exp(-alpha * abs(ii));
		temp[1][ii+rad] = -c1 * ii * alpha * alpha * exp(-alpha * abs(ii));		 
		temp[2][ii+rad] =  c2 * (1 - c3 * alpha * abs(ii)) * exp(-alpha * abs(ii));
	}
	for (ii=-rad; ii<=rad; ii++)
	{
		if (fabs(temp[0][ii+rad]) < gPara.GAUSSIAN_LOW_LIMIT) temp[0][ii+rad] = 0.;
		if (fabs(temp[1][ii+rad]) < gPara.GAUSSIAN_LOW_LIMIT) temp[1][ii+rad] = 0.;
		if (fabs(temp[2][ii+rad]) < gPara.GAUSSIAN_LOW_LIMIT) temp[2][ii+rad] = 0.;
	}

	// return back the 0-order coefficients
	if (zero != NULL)
		memcpy(zero, temp[0], len*sizeof(Float));
	
	int ig = 0;
//	TRACE("\n");
	for (kk=-rad; kk<=rad; kk++){ for (jj=-rad; jj<=rad; jj++){ for (ii=-rad; ii<=rad; ii++){
//	for (k=0; k<len; k++){ for (j=0; j<len; j++){ for (i=0; i<len; i++){
		ig = (kk+rad)*len*len + len*(jj+rad) + ii+rad;
		if (fst != NULL)
		{
			fst[0][ig] = temp[1][ii+rad]*temp[0][jj+rad]*temp[0][kk+rad]; // I_x = (f1(x)f0(y)f0(z))*I(x,y,z)
			fst[1][ig] = temp[0][ii+rad]*temp[1][jj+rad]*temp[0][kk+rad]; // I_y = (f0(x)f1(y)f0(z))*I(x,y,z)
			fst[2][ig] = temp[0][ii+rad]*temp[0][jj+rad]*temp[1][kk+rad]; // I_z = (f0(x)f0(y)f1(z))*I(x,y,z)
		}
		if (snd != NULL)
		{
			snd[0][ig] = temp[2][ii+rad]*temp[0][jj+rad]*temp[0][kk+rad]; // Ixx = (f2(x)f0(y)f0(z))*I(x,y,z)
			snd[1][ig] = temp[0][ii+rad]*temp[2][jj+rad]*temp[0][kk+rad]; // Iyy = (f0(x)f2(y)f0(z))*I(x,y,z)
			snd[2][ig] = temp[0][ii+rad]*temp[0][jj+rad]*temp[2][kk+rad]; // Izz = (f0(x)f0(y)f2(z))*I(x,y,z)
			snd[3][ig] = temp[1][ii+rad]*temp[1][jj+rad]*temp[0][kk+rad]; // Ixy = (f1(x)f1(y)f0(z))*I(x,y,z)
			snd[4][ig] = temp[1][ii+rad]*temp[0][jj+rad]*temp[1][kk+rad]; // Ixz = (f1(x)f0(y)f1(z))*I(x,y,z)
			snd[5][ig] = temp[0][ii+rad]*temp[1][jj+rad]*temp[1][kk+rad]; // Iyz = (f0(x)f1(y)f1(z))*I(x,y,z)
//			TRACE(" %f %f %f %f %f %f\n", snd[0][ig], snd[1][ig], snd[2][ig],
//				snd[3][ig], snd[4][ig], snd[5][ig]);
		}
	}}}
//	TRACE("\n");

	delete [] temp[0];
	delete [] temp[1];
	delete [] temp[2];
}

