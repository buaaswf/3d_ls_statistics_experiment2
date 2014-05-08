#pragma once
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#include <vector>
using std::vector;
#include <fstream>
using namespace std;

#include <algorithm>

//template <class T>
//void endswap(T *objp)
//{
//  unsigned char *memp = reinterpret_cast<unsigned char*>(objp);
//  std::reverse(memp, memp + sizeof(T));
//}



/////////////////////////////////////////////////////////////
// message types
#define	MSG_FILE_LOADED				(WM_USER+MSG_NUM_UI)
#define MSG_MOUSEMOVE_VIEW_XYZ		(WM_USER+MSG_NUM_UI+1)
#define MSG_CANDIDATE_SEL			(WM_USER+MSG_NUM_UI+2)

/////////////////////////////////////////////////////////////
// data types
#define RAW_DATA_TYPE	short
#define SEG_DATA_TYPE	BYTE
//#define Float			double
#define Float			float

#define THRESHOLD_MUCOSA_MIN		0
#define THRESHOLD_MUCOSA_MAX		200

#define BYTE_BITS					8

typedef enum VXL_TYPE
{
	VXL_UNKNOWN			= 0,	// unknown type
	VXL_MUCOSA			= 1,	// mucosa layer (voxels valuated 1~199)
	VXL_AIR				= 2,	// pure air (voxels valuated 200)
	VXL_TISSUE			= 4,	// pure tissue (voxels valuated 0)
	VXL_TISSUE_BORDER	= 8,	// tissue border immediately adjacent to mucosa (part of pure tissue, voxels valuated 0)
	VXL_AIR_BORDER		= 16,	// air border immediately adjacent to mucosa (part of pure tissue, voxels valuated 200)
	VXL_0_LEVEL			= 32,	// voxels on the 0-level set
	VXL_INNER			= 64,	// voxels inside the 0-level set
	VXL_FLAG			= 128	// voxels flag, used to flag voxels
}* PVXL_TYPE;
#define VXL_EDM		(VXL_AIR|VXL_MUCOSA) // |VXL_TISSUE_BORDER
typedef CMap<int, int, int, int> IndexMap;

/////////////////////////////////////////////////////////////
// tool functions
#define D_MAX			(double(1.0e+18))
#define F_MAX			(float(1.0e+10))

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
void calGaussian(Float alpha, int number[3], Float table[3][100]);
// another version
void calGaussian(Float alpha, int len, Float* zero, Float* fst[3], Float* snd[6]);


/////////////////////////////////////////////////////////////
// indices
extern int i, j, k;
extern int NX, NY, NZ, dj, dk, size;
extern float FX, FY, FZ, RZ;

extern Float djInv, dkInv;
#define	FOR_THREE		for (k=0; k<NZ; k++) { for (j=0; j<NY; j++) { for (i=0; i<NX; i++) {
#define FOR_EACH		for (i=0; i<size; i++)
#define END_THREE		}}}
#define	FOR_YX			for (j=0; j<NY; j++) { for (i=0; i<NX; i++) {
#define FOR_ZY			for (k=0; k<NZ; k++) { for (j=0; j<NY; j++) {
#define FOR_ZX			for (k=0; k<NZ; k++) { for (i=0; i<NX; i++) {
#define END_TWO			}}
#define GI(i,j,k)		((i) + dj*(j) + dk*(k))
#define inX(i)			((i)>=0 && (i)<NX)
#define inY(j)			((j)>=0 && (j)<NY)
#define inZ(k)			((k)>=0 && (k)<NZ)
#define inDomain(i,j,k)	(inX((i)) && inY((j)) && inZ((k)))

inline void GIJK(int n, int& i, int& j, int& k)
{	
	k = int(n * dkInv);		j = Mod(n,dk,dkInv);
	i = Mod(j,dj,djInv);	j = int(j * djInv); 
}

// index for the 26 neighbors of (x,y,z)
#define ZZZ (GI(x, y, z)) // it's me

#define	MZZ (GI(x-1, y, z))
#define	PZZ (GI(x+1, y, z))
#define	ZMZ (GI(x, y-1, z))
#define	ZPZ (GI(x, y+1, z))
#define	ZZM (GI(x, y, z-1))
#define	ZZP (GI(x, y, z+1))

#define	ZMP (GI(x, y-1, z+1))
#define	ZMM (GI(x, y-1, z-1))
#define	ZPP (GI(x, y+1, z+1))
#define	ZPM (GI(x, y+1, z-1))
#define	MZP (GI(x-1, y, z+1))
#define	MZM (GI(x-1, y, z-1))
#define	PZP (GI(x+1, y, z+1))
#define	PZM (GI(x+1, y, z-1))
#define	PMZ (GI(x+1, y-1, z))
#define	MMZ (GI(x-1, y-1, z))
#define	PPZ (GI(x+1, y+1, z))
#define	MPZ (GI(x-1, y+1, z))

#define	MMM (GI(x-1, y-1, z-1))
#define	MMP (GI(x-1, y-1, z+1))
#define	MPM (GI(x-1, y+1, z-1))
#define	MPP (GI(x-1, y+1, z+1))
#define	PMM (GI(x+1, y-1, z-1))
#define	PMP (GI(x+1, y-1, z+1))
#define	PPM (GI(x+1, y+1, z-1))
#define	PPP (GI(x+1, y+1, z+1))

/////////////////////////////////////////////////////////////
// header for raw file
class RAW_HEADER
{
public:
	RAW_HEADER()	{nX = nY = nZ = 0; fX = fY = fZ = 0.;}
	~RAW_HEADER()	{}
	void	clear() {nX = nY = nZ = 0; fX = fY = fZ = 0.;}

	int		nX, nY, nZ;			// resolutions in x,y,z-direction
	float	fX, fY, fZ;			// actual lengths (in mm) for each voxel in x,y,z-direction

public:
	void	load(CFile& file);
};

/////////////////////////////////////////////////////////////
// header for seg file
class SEG_HEADER
{
public:
	SEG_HEADER()	{nX = nY = nZ = 0; f = 0; }
	~SEG_HEADER()	{}
	void	clear() {nX = nY = nZ = 0; f = 0; }
	int		nX, nY, nZ;		// resolutions in x,y,z-direction
	float	f;

public:
	void	load(CFile& file);
	void	write(CFile& file);
	
};

class GEO_HEADER
{
public:
	GEO_HEADER()
	{
		bOkGrad=bOkCurv=bOkLICC=bOkSICN=FALSE;num=0;
		maxGrad = maxCurv = maxCurvedness = -F_MAX;
		minGrad = minCurv = minCurvedness =  F_MAX;
	}
	~GEO_HEADER(){}

	void	clear() 
	{
		bOkGrad=bOkCurv=bOkLICC=bOkSICN=FALSE;num=0;
		maxGrad = maxCurv = maxSI = maxCurvedness = -F_MAX;
		minGrad = minCurv = minSI = minCurvedness =  F_MAX;
	}

	void	load(CFile& file);
	void	write(CFile& file);
	void    write_new(ofstream &out );
public:
	BOOL	bOkGrad;
	BOOL	bOkCurv;
	BOOL	bOkLICC;
	BOOL	bOkSICN;
	int		num;
	Float	maxGrad, minGrad, maxCurv, minCurv, 
			minSI, maxSI, minCurvedness, maxCurvedness;
};

#define MAX_NUM_TP	10
class PLP_HEADER
{
public:
	PLP_HEADER()
	{
		nIPC = nTP = 0;
	}
	~PLP_HEADER(){}

	void	load(CFile& file);
	void	write(CFile& file);

public:
	int		nIPC;
	int		nTP;				// number of TPs
	int		dtcIPC[MAX_NUM_TP];	// detection result
};

typedef enum GEO_TYPE
{
	GRAD_X			= 0,	// x-component of gradient
	GRAD_Y			= 1,	// y-component of gradient
	GRAD_Z			= 2,	// z-component of gradient
	GRAD_XYZ		= 3,	// norm of gradient
	FST_CURV		= 4,	// k1
	SND_CURV		= 5,	// k2
	FST_LIC_CURV	= 6,	// k1' after smoothing along gradient direction in the mucosa layer (CRM paper)
	SND_LIC_CURV	= 7,	// k2' after smoothing along gradient direction in the mucosa layer (CRM paper)
	SI_ORG			= 8,	// shape index from k1 and k2
	CV_ORG			= 9,	// curvedness from k1 and k2
	SI_SMOOTH		= 10,	// shape index from k1' and k2'
	CV_SMOOTH		= 11,	// curvedness from k1' and k2'
	CAN_ID			= 12,	// the IPC this voxel belongs to
	PHI				= 13,	// the phi value on this voxel for level set (CRM paper)
	GEO_NUM			= 14	// 
};


/////////////////////////////////////////////////////////////
// geometric and other information on each voxel
class MucosaVoxel
{
public:
	MucosaVoxel() {memset(m_data, 0, GEO_NUM*sizeof(Float)); m_idCan = -1; /*adapt = NULL;*/}
	~MucosaVoxel(){/*clearAdapt();*/}

	Float&	operator [] (GEO_TYPE t)		{return m_data[t];}
	Float	operator [] (GEO_TYPE t) const	{return m_data[t];}
	Float*	data()							{return m_data;	}
	int		idCan()							{return m_idCan;}

public:
	Float	m_data[GEO_NUM];// the geometry information
	int		m_idCan;		// indicates to which this voxel belongs
};


typedef enum FEAT_TYPE
{
	// si and cv are computed on seeds.
	MEAN_SI		= 0,
	VAR_SI		= 1,
	SKEW_SI		= 2,
	KURT_SI		= 3,
	ENTROPY_SI	= 4,
	MEAN_CV		= 5,
	VAR_CV		= 6,
	SKEW_CV		= 7,
	KURT_CV		= 8,
	ENTROPY_CV	= 9,
	
	// ct are computed within the whole IPC on original un-cleansed image
	MEAN_CT		= 10,
	VAR_CT		= 11,
	SKEW_CT		= 12,
	KURT_CT		= 13,
	ENTROPY_CT	= 14,
	
	// features relative to volume
	VOL_SEEDS	= 15,
	VOL_SURF	= 16,
	VOL_RATIO	= 17,
	SG_RATIO	= 18,
	SG_NUM		= 19,
	
	// morphology feature
	AXIS_RATIO	= 20,

	// features in axial gray image
	AXI_HR		= 21,	// highlighting ratio
	AXI_DL		= 22,	// disk-likeness 
	
	// features in coronal gray image
	COR_LR		= 23,	// lightness ratio
	COR_GCX		= 24,	// normalized cx of gray patch 
	COR_GCY		= 25,	// normalized cy of gray patch
	COR_BCX		= 26,	// normalized cx of bright patch 
	COR_BCY		= 27,	// normalized cy of bright patch

	// features in sagittal gray image
	SAG_LR		= 28,	// lightness ratio
	SAG_GCX		= 29,	// normalized cx of gray patch
	SAG_GCY		= 30,	// normalized cy of gray patch
	SAG_BCX		= 31,	// normalized cx of bright patch
	SAG_BCY		= 32,	// normalized cy of bright patch

	// feature in axial color image
	CLR_AH		= 33,
	CLR_AS		= 34,
	CLR_AI		= 35,

	// for debugging purpose
	TEMP_00		= 36,
	TEMP_01		= 37,
	TEMP_02		= 38,
	TEMP_03		= 39,
	TEMP_04		= 40,
	TEMP_05		= 41,
	TEMP_06		= 42,
	TEMP_07		= 43,
	TEMP_08		= 44,
	TEMP_09		= 45,
	TEMP_10		= 46,
	TEMP_11		= 47,
	TEMP_12		= 48,
	TEMP_13		= 49,
	TEMP_14		= 50,
	TEMP_15		= 51,
	TEMP_16		= 52,
	TEMP_17		= 53,
	TEMP_18		= 54,
	TEMP_19		= 55,
	TEMP_20		= 56,
	TEMP_21		= 57,
	TEMP_22		= 58,
	TEMP_23		= 59,
	TEMP_24		= 60,
	TEMP_25		= 61,
	TEMP_26		= 62,
	TEMP_27		= 63,
	TEMP_28		= 64,
	TEMP_29		= 65,
	TEMP_30		= 66,
	TEMP_31		= 67,
	TEMP_32		= 68,
	TEMP_33		= 69,
	TEMP_34		= 70,
	TEMP_35		= 71,
	TEMP_36		= 72,
	TEMP_37		= 73,
	TEMP_38		= 74,
	TEMP_39		= 75,
	TEMP_40		= 76,
	TEMP_41		= 77,
	TEMP_42		= 78,
	TEMP_43		= 79,
	TEMP_44		= 80,

	FEAT_NUM	= 81

};


// extensions for outputting the features into ascii files, so that matlab can be used to plot
const CString FEAT_NAME = "features";
const CString FEAT_EXT[] = 
{
	".sim",		// mean(SI)
	".siv",		// var(SI)
	".sis",		// skew(SI)
	".sik",		// Kurt(SI)
	".sie",		// Entropy(SI)
	".cvm",		// mean(CV)
	".cvv",		// var(CV)
	".cvs",		// Skew(CV)
	".cvk",		// Kurt(CV)
	".cve",		// Entropy(CV)
	".ctm",		// mean(CT)
	".ctv",		// var(CT)
	".cts",		// Skew(CT)
	".ctk",		// Kurt(CT)
	".cte",		// Entropy(CT)

	".vsd",		// volume of seeds
	".vsf",		// volume of surface
	".vsgr",	// volume ratio of seed over grow
	".sgr",		// number ratio of seed over grow
	".sgn",		// number of seeds and grows voxels

	".aro",		// the axis ratio

	".ahr",		// highlighting ratio
	".adl",		// disk-likeness
	".clr",		// lightness ratio in coronal image
	".cgx",		// cx of gray patch in coronal image
	".cgy",		// cy of gray patch in coronal image
	".cbx",		// cx of bright patch in coronal image
	".cby",		// cy of bright patch in coronal image
	".slr",		// lightness ratio in sagittal image
	".sgx",		// cx of gray patch in sagittal image
	".sgy",		// cy of gray patch in sagittal image
	".sbx",		// cx of bright patch in sagittal image
	".sby",		// cy of bright patch in sagittal image

	".cah",		// h-component in color core of axial image
	".cas",		// s-component in color core of axial image
	".cai",		// i-component in color core of axial image

	".tmp00",
	".tmp01",
	".tmp02",
	".tmp03",
	".tmp04",
	".tmp05",
	".tmp06",
	".tmp07",
	".tmp08",
	".tmp09",
	".tmp10",
	".tmp11",
	".tmp12",
	".tmp13",
	".tmp14",
	".tmp15",
	".tmp16",
	".tmp17",
	".tmp18",
	".tmp19",
	".tmp20",
	".tmp21",
	".tmp22",
	".tmp23",
	".tmp24",
	".tmp25",
	".tmp26",
	".tmp27",
	".tmp28",
	".tmp29",
	".tmp30",
	".tmp31",
	".tmp32",
	".tmp33",
	".tmp34",
	".tmp35",
	".tmp36",
	".tmp37",
	".tmp38",
	".tmp39",
	".tmp40",
	".tmp41",
	".tmp42",
	".tmp43",
	".tmp44",

	//
	".tf"		// FP(0)/TP(1)
};

// Designed for gray stuff
#define	GRY_SEED	0x1
#define	GRY_GROW	0x2
#define GRY_BGD		0x4	// used to represent the "background" 
// Designed for colorful stuff
#define CLR_ROI		GRY_SEED
// Pixel types from the 3D voxels
#define PLP_SEED	0x8
#define PLP_GROW	0x10
#define PLP_EXT		0x20
#define PLP_OUT		0x40

/////////////////////////////////////////////////////////////
// Class candidate
class Candidate
{
public:
	Candidate() { ID = genID(); num = 0; TPID = -1;
		ZeroMemory(feat, FEAT_NUM*sizeof(float));
		ZeroMemory(noml, 3*sizeof(float)); m_bCADTP = FALSE; }
	~Candidate(){}

	static bool		RED(Float H, Float S, Float I);
	static bool		WHITE(Float H, Float S, Float I);
	static bool		GREEN(Float H, Float S, Float I);
public:
	inline int		id()						{return ID;			}
	inline void		id(int d)					{ID = d;			}
	inline void		cid(int d)					{CID= d;			}
	inline int		cid()						{return CID;		}
	inline int		getNum()					{return num;		}
	inline int		getVxl(int n)				{return vxls[n];	}
	inline void		getVxls(vector<int>& rs)	{rs.insert(rs.end(), vxls.begin(), vxls.end());}
	inline int		getCenter()					{return center;		}
	inline float&	feature(FEAT_TYPE t)		{return feat[t];	}
	inline float	feature(FEAT_TYPE t) const	{return feat[t];	}
	inline bool		isTruth()					{return TPID>0?true:false;}
	inline void		setTPID(int id=-1)			{TPID = id;			}
	inline int		getTPID()					{return TPID;		}
	inline void		setNoml(float* nrl)			{memcpy(noml, nrl, 3*sizeof(float)); }
	inline void		getNoml(float* nrl)			{memcpy(nrl, noml, 3*sizeof(float)); }
	inline float*	getLFBox()					{return lfBox;		}
	inline int*		getOFBox()					{return ofBox;		}

	bool			isInShowBox(int x, int y, int z);
	float			volEstimation(const IndexMap& mpMV, float* mv);
	void			save(CFile& file);
	void			load(CFile& file);
	void			saveCAD(CStdioFile& file, float zf);
	void			setVxlCanID(const IndexMap& mpMV, float* MV);
	void			setupLocalFrame();
	void			genBiopsyFeats(RAW_DATA_TYPE* cleanRaw, RAW_DATA_TYPE* dirtyRaw, const CString& path);

	typedef vector<int>::iterator Iterator;
	typedef vector<int>::const_iterator cIterator;
	cIterator begin() const		{return vxls.begin();}
	cIterator end()   const		{return vxls.end(); }

	inline void	grow(int index, BYTE t);
	inline void	merge(Candidate* other);
	inline void setPlpVxls()	{ plpVxls.insert(plpVxls.end(), vxls.begin(), vxls.end()); }
	inline void	calcCenter();
	inline CMatrix4X4* getLocalFrame(){return &mtx;}

	static void clearID(int d=0){canID = d;}
	static int	getMaxID()		{return canID;}
	static void writeFeatursMat(Candidate* c);
	
	static CStdioFile	featFiles[FEAT_NUM+2];
	static bool			isFeatFilesOpen;
	static void			openFeatFiles(const CString& sPath);
	static void			closeFeatFiles();

	BOOL	m_bCADTP;
protected:
	void	genBiopsyImgs(RAW_DATA_TYPE* cleanRaw, RAW_DATA_TYPE* dirtyRaw, 
		Grid2dF& gA, Grid2dF& gC, Grid2dF& gS, 
		Grid2d4v& cA, Grid2d4v& cC, Grid2d4v& cS,
		Grid2dU& tA, Grid2dU& tC, Grid2dU& tS);
	void	proGrayImg(Grid2dF& gA, Grid2dF& gC, Grid2dF& gS, 
		Grid2dU& tA, Grid2dU& tC, Grid2dU& tS);
	void	proClrImg(Grid2d4v& cA, Grid2d4v& cC, Grid2d4v& cS,
		Grid2dU& edgeA, Grid2dU& edgeC, Grid2dU& edgeS);

	//
	void	genGrayFeats(Grid2dU& tA, Grid2dU& tC, Grid2dU& tS,
		Grid2dF& gA, Grid2dF& gC, Grid2dF& gS);
	void	genCCASeeds(Grid2dF& g, Grid2dU& t, unsigned char s);
	void	genClrFeats(Grid2d4v& c, Grid2dU& t, unsigned char or);
protected:
	static int	canID;
	static int	genID() {return ++canID;}

private:
	int			ID;
	int			CID;	// indicates which candidate = "case id + 0/1(Prone/Supine) + candidate ID", like 30370001, means wi_3037_P, candidate #01.
	vector<int> vxls;	// voxel indices array
	int			num;	// number of voxels
	vector<int> plpVxls;// voxels must be belong to the polyp
	vector<BYTE> type;

	// the center and the normal vector -> the origin and z-axis of the local frame
	int			center;	// index of the central voxel
	float		noml[3];// normal vector

	CMatrix4X4	mtx, mV; // mtx: mapping points to the local frame; mV: mapping vector to the local frame
	int			ofBox[6];// the bounding box in the original frame
	float		lfBox[6];// the bounding box in the local frame

	int			TPID;	 // the #ID in polyp database
	float		feat[FEAT_NUM];
};

////////////////////////////////////////////////////////////////////////////////////
// grow by one voxel
////////////////////////////////////////////////////////////////////////////////////
inline void Candidate::grow(int index, BYTE t)
{
	vxls.push_back(index);
	type.push_back(t);
	num ++;
}

////////////////////////////////////////////////////////////////////////////////////
// merge the candidate of "other" into this candidate
////////////////////////////////////////////////////////////////////////////////////
inline void Candidate::merge(Candidate* other)
{
	vxls.insert(vxls.end(), other->vxls.begin(), other->vxls.end());
	type.insert(type.end(), other->type.begin(), other->type.end());
	num += other->num;
}

inline void	Candidate::calcCenter()
{
	int x = 0, y = 0, z = 0, x0, y0, z0;
	vector<int>::iterator it;
	for (it=vxls.begin(); it<vxls.end(); it++)
	{
		GIJK(*it, x0, y0, z0);
		x += x0; y += y0; z += z0;
	}
	x /= num; y /= num; z /= num;
	center = GI(x, y, z);
}

/////////////////////////////////////////////////////////////
// geometric and other information on each voxel
class GlobalPara
{
public:
	GlobalPara()
	{
		maxSeedSI = 0.09f;
		minSeedSI = 0.f;

		maxSeedCN = 0.19f;
		minSeedCN = 0.05f;

		maxGrowSI = 0.25f;
		minGrowSI = 0.f;

		maxGrowCN = 0.5f;
		minGrowCN = 0.02f;

		mergeDist = 15.f;	// 10 mm
		minVol	  = 5;
		maxVol	  = 2000;

		ALPHA_0TH = 1.5f;
		ALPHA_1ST = 1.f;
		ALPHA_2ND = 0.8f;

		GAUSSIAN_SIZE_1ST = 19;
		GAUSSIAN_SIZE_2ND = 27;
		GAUSSIAN_LOW_LIMIT= 0.00005f;
		GAUSSIAN_MAX_LEN  = 27;

		FMM_LIMIT_P = 2.5f;
		FMM_LIMIT_N = 23.5f;

		ifUseAdapt	= true;
		ifUseLIC	= true;
	}

public:
	// used for clustering candidates
	Float	maxSeedCN, minSeedCN, maxGrowCN, minGrowCN;
	Float	maxSeedSI, minSeedSI, maxGrowSI, minGrowSI;

	Float	mergeDist;
	int		minVol, maxVol;

	Float	ALPHA_0TH, ALPHA_1ST, ALPHA_2ND;
	int		GAUSSIAN_SIZE_1ST, GAUSSIAN_SIZE_2ND, GAUSSIAN_MAX_LEN;

	Float	FMM_LIMIT_P, FMM_LIMIT_N;	// used for fast marching method
	Float	GAUSSIAN_LOW_LIMIT;

	bool	ifUseAdapt, ifUseLIC;

};

extern GlobalPara gPara;

inline bool isSeed(float si, float cv)
{
	if (si>gPara.minSeedSI&&si<gPara.maxSeedSI &&
		cv>gPara.minSeedCN&&cv<gPara.maxSeedCN )
		return true;
	else return false;
}
inline bool isGrow(float si, float cv)
{
	if (isSeed(si,cv)) return false;
	else
	{
		if (si>gPara.minGrowSI&&si<gPara.maxGrowSI &&
			cv>gPara.minGrowCN&&cv<gPara.maxGrowCN )
			return true;
		else return false;
	}
}
inline bool	isSG(float si, float cv)
{
	if (si>gPara.minGrowSI&&si<gPara.maxGrowSI &&
		cv>gPara.minGrowCN&&cv<gPara.maxGrowCN )
		return true;
	else return false;
}

#endif
