#pragma once
#include "vol_math_LevelSet.h"
#include "CImg.h" 
#include "vol_math_RawImage.h"
#include <vector>
class ThreeDim_LevelSet
{
public:
	ThreeDim_LevelSet(void);
	~ThreeDim_LevelSet(void);
	Raw minimal_surface(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction );

	void array_surface(Raw *src);
	void outerwall(Raw &src,Raw &phi,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void outerwallauto(Raw &src,Raw &phi,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void couple(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction);
	void NeumannBoundCond( Raw &phi );
	void Narrowband (Raw &phi);
	Raw ImageFSqrt( Raw &phi_x, Raw &phi_y,Raw &phi_z );
	void initialg(Raw &raw);
};

