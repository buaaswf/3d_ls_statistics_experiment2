#pragma once
#include "vol_math_RawImage.h"
#include "vol_math_Raw3D_Independt.h"
#include <vector>
class Point
{
	
public :
	long long xindex;
	long long yindex;
	long long zindex;
	PIXTYPE value;
public:
	Point (long long x,long long y, long long z,PIXTYPE val)
	{
		this->xindex = x;
		this->yindex = y;
		this->zindex = z;
		this->value = val;

	}
	
};
class Fast_ls
{
public:
	vector<Point> narrow;
	Raw initial;//the data to be segmented.
	Fast_ls(void);
	~Fast_ls(void);
	Raw minimal_surface(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction );
	void array_surface(Raw *src);
	void outerwall(Raw &src,Raw &phi,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void outerwallauto(Raw &src,Raw &phi,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void couple(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void NeumannBoundCond( Raw &phi );
	void Narrowband (Raw &phi);
	Raw ImageFSqrt( Raw &phi_x, Raw &phi_y,Raw &phi_z );
	vector<Point>ImageFSqrt( vector<Point> phi_x, vector<Point> phi_y,vector<Point> phi_z );
	void initialg(Raw &raw);
	Raw minimal_surface(Raw &src,vector <Point> phi,Raw &g,double lambda,double mu,
		double alfa,float epsilon,int timestep,int iter,char *potentialFunction );
	vector<Point> findNarrow(Raw *src);
};

