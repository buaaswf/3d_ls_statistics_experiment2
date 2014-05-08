#include "ThreeDim_LevelSet.h"
#include <iostream>
#define pi 3.141592653
#include<math.h>
using namespace cimg_library;
using namespace std;
//PIXTYPE *p;
//PIXTYPE	*q;
//int i,j,k;
Raw del2( Raw &phi ) 
{
	int m=phi.getXsize();
	int n=phi.getYsize();
	int l=phi.getZsize();
	Raw ret2(m,n,l);
	for (int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			for (int k=0;k<l;k++)
			{
				if (i+1< m && j+1< n &&k+1< l && i-1>=0 && j-1 >= 0 && k-1>= 0)
				{
					PIXTYPE value = (1.0/6.0)*(phi.get(i+1, j, k) + phi.get(i-1, j, k) + phi.get(i, j-1, k) + phi.get(i, j+1, k) + phi.get(i,j,k+1)+phi.get(i,j,k-1)- 6*(phi.get(i,j,k)));
					ret2.put(i, j,k, (PIXTYPE)value);
				}
				else 
				{
					ret2.put(i, j,k,0);
				}
			}

			
		}
	}

	return ret2.set_shared(true);
}

Raw gradientxgc( Raw &g ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	Raw ret(g);
	int i,j,k,temp1,temp2;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			for ( k=0;k < l;k++)
			{
				if(i>0)
					temp1=i-1;
				else
					temp1=0;
				if (i<n-1)
					temp2=i+1;
				else 
					temp2=n-1;
				ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
				//if (ret.get(i,j,k)!=0)
				//{
				//	cout<<"i="<<i<<",j="<<j<<",k="<<k<<ret.get(i,j,k)<<endl;
				//}
			}			
		}
	}
	return ret.set_shared(true);
}

Raw gradientygc( Raw & g ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	Raw ret(n, m,l);
	int i,j,k,temp1,temp2;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			for ( k = 0;k < l;k++)
			{
				if(j>0)
					temp1=j-1;
				else
					temp1=0;
				if (j<n-1)
					temp2=j+1;
				else 
					temp2=m-1;
				ret.put(i,j,k,0.5*(g.get(i,temp2,k)-g.get(i,temp1,k)));
			}			
		}
	}
	return ret.set_shared(true);
	
}
Raw  gradientzgc( Raw &g ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	Raw ret(n, m,l);
	int i,j,k,temp1,temp2;
	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			for ( k = 0;k < l;k++)
			{
				if(k>0)
					temp1=k-1;
				else
					temp1=0;
				if (k<l-1)
					temp2=k+1;
				else 
					temp2=l-1;
				ret.put(i,j,k,0.5*(g.get(i,j,temp2)-g.get(i,j,temp1)));
			}			
		}
	}
	return ret.set_shared(true);
}
Raw cos(Raw &x)
{
	Raw ret(x);
	for (int i=0;i<x.getXsize();i++)
	{
		for (int j=0;j< x.getYsize();j++)
		{
			for (int k=0;k< x.getZsize();k++)
			{
				ret.put(i,j,k,cos((double)x.get(i,j,k)));
			}
		}
	}
	return ret.set_shared(true);
}
Raw div(Raw &x,Raw &y,Raw &z)
{
	return (gradientxgc(x) += gradientygc(y) += gradientzgc(z)).set_shared(true);
}
void ThreeDim_LevelSet::array_surface(Raw *src)
{
	for (int i=0;i<src->getZsize();i++)
	{

		
	}
}
Raw regFunction(Raw &s,double m,double n)
{
	int x = s.getXsize(),
		y = s.getYsize(),
		z = s.getZsize();
	Raw ss(x,y,z);
	PIXTYPE val=0;
	for (int i=0;i<x;i++)
	{
		for(int j=0;j<y;j++)
		{
			for (int k= 0;k<z;k++)
			{
				val=s.get(i,j,k);
				if(val>=m && val<=n )
				{
					//ss->put(i,j,255);//unsigned char version
					ss.put(i,j,k,1);
				}
				else if(s.get(i,j,k)==m||s.get(i,j,k)==n) 
				{
					ss.put(i,j,k,0);
				}
				else 
				{
					ss.put(i,j,k,0);
				}
			}
			
		}
	}

	return ss.set_shared(true);
}
Raw  Dirac( Raw & x, double sigma ) 
{
	Raw ret(x);
	PIXTYPE temp=((1.0/2.0)/sigma);
	//double temp2=(cos((2/sigma)*pi)+1)*temp;
	Raw f= (cos((ret/sigma)*pi)+1)*temp;
	Raw b = regFunction(ret, -sigma, sigma);
	//IShowImg(b);
	ret = f*b;
	//IShowImg(ret);

	return ret.set_shared(true);
}
Raw sin(Raw & src)
{
	for (int i = 0; i < src.getXsize(); i++)
	{
		for (int j = 0; j < src.getYsize(); j++)
		{
			for (int k = 0; k < src.getZsize(); k++)
			{
				src.put( i, j, k, sin(src.get( i, j ,k)));
			}
		}
	}
	return src;
}
Raw distReg_p2( Raw  &phi ) 
{
	ThreeDim_LevelSet *ls;
		int m=phi.getXsize();
		int n=phi.getYsize();
		Raw phi_x(gradientxgc(phi));
		Raw phi_y(gradientygc(phi));
		Raw phi_z(gradientzgc(phi));
		Raw s(ls->ImageFSqrt(phi_x,phi_y,phi_z));
		Raw &a=regFunction(s,0,1);
		Raw &b=regFunction(s,1,10);//need to be changed.
		Raw &ps=a*(sin((s*2.0*pi)))/(2.0*pi)+b*(s-1);
		Raw &dps=((regFunction(ps,0,0)*(-1)+1)*ps+regFunction(ps,0,0))/((regFunction(s,0,0)*(-1)+1)*s+regFunction(s,0,0));
		Raw &div2=div(dps*phi_x-phi_x,dps*phi_y-phi_y,dps*phi_z-phi_z);
		Raw &f=div2-del2(phi)*4.0;
		return f;

}
void ThreeDim_LevelSet::initialg(Raw &g)
{
	Raw gx=gradientxgc(g);
	Raw gy=gradientygc(g);
	Raw gz=gradientzgc(g);
	gx*=gx;
	gy*=gy;
	gz*=gz;
	gx+= gy += gz;
	gx += 1.0;
	g = 1.0/(gx);
}
Raw ThreeDim_LevelSet::minimal_surface(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction )
{
	int m=g.getXsize();
	int n=g.getYsize();
	int l=g.getZsize();
	Raw vx=gradientxgc(g);
	Raw vy=gradientygc(g);
	Raw vz=gradientzgc(g);
	for(int i=0;i<iter;i++)
	{
		float smallNumber=1e-10;

		NeumannBoundCond(phi);
		Raw phi_x = gradientxgc(phi);
		Raw phi_y = gradientygc(phi);
		Raw phi_z = gradientzgc(phi);
		Raw s = ImageFSqrt(phi_x, phi_y, phi_z) += smallNumber;

		phi_x /= s;
		phi_y /= s;
		phi_z /= s;

		s.~Raw();

		Raw curvature = div(phi_x, phi_y, phi_z);
	
		char *p1="single_well";
		Raw distRegTerm;
		if (0 == strcmp(potentialFunction, p1))
		{
			/*
			compute distance regularization term in equation (13) 
			with the single-well potential p1.
			*/
			distRegTerm = (del2(phi) *= 6.0) -= curvature;
		} else if (0 == strcmp(potentialFunction, "double_well")) {
			distRegTerm = distReg_p2(phi);  // compute the distance regularization term in eqaution (13) with the double-well potential p2.
		} else {
			cout << "EEROR" << endl;
		}

		phi_x *= vx;
		phi_y *= vy;
		phi_z *= vz;


		phi_x += phi_y += phi_z += g*curvature;
		phi_x *= lambda;
		phi_x += g*alfa;
		phi_x *= Dirac(phi,epsilon);

		phi += (distRegTerm)*mu*((double)timestep);
		phi += phi_x;
		cout<<"iterator i="<<i<<endl;
	}	
	return phi;
}
Raw outwallpull(Raw &src)
{
	Raw pull(src);
	for (int i=0;i<src.getXsize();i++)
	{
		for (int j=0;j<src.getYsize();j++)
		{
			for (int k=0;k<src.getZsize();k++)
			{
				if (i<3&&j<3&&k<3)
				{
					pull.put(i,j,k,100);
				} 
				else
				{
					pull.put(i,j,k,0);
				}
				

			}
		}
	}



	return pull;
}



void ThreeDim_LevelSet::outerwall(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
{

	this->initialg(g);
	Raw pull=outwallpull(phi);
	this->initialg(pull);
	for (int i=0; i < 5; i++)
	{
		phi = this->minimal_surface(phi,pull+g,lambda,mu,-10,epsilon,timestep,8,potentialFunction)+\
			this->minimal_surface(phi,g,lambda,mu,-alfa,epsilon,timestep,1,potentialFunction);

		cout << "outer wall iter = " << i <<endl;
	}

}

void ThreeDim_LevelSet::outerwallauto(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
{

	this->initialg(g);
	for (int i=0; i < 30; i++)
	{
		phi=this->minimal_surface(phi,g,lambda,mu,alfa,epsilon,timestep,1,potentialFunction);
		cout << "outer wall auto stop iter = " << i <<endl;
	}

}


void ThreeDim_LevelSet::couple(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
{

}
ThreeDim_LevelSet::ThreeDim_LevelSet(void)
{



}


ThreeDim_LevelSet::~ThreeDim_LevelSet(void)
{
}

void ThreeDim_LevelSet::NeumannBoundCond( Raw &img )
{
	int nrow=img.getXsize();
	int ncol=img.getYsize();
	int ndepth=img.getZsize();
	int i=0,j=0,k=0;
	//the eight point SDF
	img.put(0,0,0,img.get(3,3,3));
	img.put(nrow-1,0,0,img.get(nrow-3,3,3));
	img.put(0,ncol-1,0,img.get(3,ncol-3,3));
	img.put(0,0,ndepth-1,img.get(3,3,ndepth-3));
	img.put(nrow-1,ncol-1,0,img.get(nrow-3,ncol-3,3));
	img.put(nrow-1,0,ndepth-1,img.get(nrow-3,3,ndepth-3));
	img.put(0,ncol-1,ndepth-1,img.get(3,ncol-3,ndepth-3));
	img.put(nrow-1,ncol-1,ndepth-1,img.get(nrow-3,ncol-3,ndepth-3));
	//first and the last column SDF
	for(i=2;i<nrow-2;i++)
	{
		img.put(i,0,0,img.get(i,3,3));
		img.put(i,ncol-1,0,img.get(i,ncol-3,3));
		img.put(i,ncol-1,ndepth-1,img.get(i,ncol-3,ndepth-3));
		img.put(i,0,ndepth-1,img.get(i,3,ndepth-3));
	}
	//first and last row SDF
	for(j=2;j<ncol-2;j++)
	{
		img.put(0,j,0,img.get(3,j,3));
		img.put(0,j,ndepth-1,img.get(3,j,ndepth-3));
		img.put(nrow-1,j,0,img.get(nrow-3,j,3));
		img.put(nrow-1,j,ndepth-1,img.get(nrow-3,j,ndepth-3));
	}
	//first and last depth SDF
	for(k=2;k<ndepth-2;k++)
	{
		img.put(0,0,k,img.get(3,3,k));
		img.put(nrow-1,0,k,img.get(nrow-3,3,k));
		img.put(0,ncol-1,k,img.get(3,ncol-3,k));
		img.put(nrow-1,ncol-1,k,img.get(nrow-3,ncol-3,k));
	}
	//front and back surface SDF
	for(i=2;i<nrow-2;i++)
	{
		for (k=2;k<ndepth-2;k++)
		{
			img.put(i,0,k,img.get(i,3,k));
			img.put(i,ncol-1,k,img.get(i,ncol-3,k));

		}

	}
	//up and below surface SDF
	for(i=2;i<nrow-2;i++)
	{
		for (j=2;j<ncol-2;j++)
		{
			img.put(i,j,0,img.get(i,j,3));
			img.put(i,j,ndepth-1,img.get(i,j,ndepth-3));

		}

	}
	//left and right surface SDF
	for(j=2;j<ncol-2;j++)
	{
		for (k=2;k<ndepth-2;k++)
		{
			img.put(0,j,k,img.get(3,j,k));
			img.put(nrow-1,j,k,img.get(nrow-3,j,k));

		}

	}

	//for (i=0;i<nrow;i++)
	//{
	//	for(j=0;j<ncol;j++)
	//	{
	//		img.put(i,j,k,img.get(i,j,k));
	//	}
	//}
}

Raw ThreeDim_LevelSet::ImageFSqrt( Raw &phi_x, Raw &phi_y,Raw &phi_z )
{
	//return sqrt(phi_x*phi_x+phi_y*phi_y+phi_z*phi_z);
	int m=phi_x.getXsize();
	int n=phi_x.getYsize();
	int l=phi_x.getZsize();
	Raw ret(m,n,l);
	for (int i=0;i< m; i++)
	{
		for ( int j = 0; j< n; j++)
		{
			for ( int k = 0; k< l; k++)
			{
				ret.put(i,j,k,sqrt((double)phi_x.get(i,j,k)*phi_x.get(i,j,k)+phi_y.get(i,j,k)*phi_y.get(i,j,k)+phi_z.get(i,j,k)*phi_z.get(i,j,k)));
			}
		}
	}

	return ret.set_shared(true);
}

void Narrowband (Raw &phi)
{

	
}