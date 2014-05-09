#include "Fast_ls.h"
#define pi 3.141592653


Fast_ls::Fast_ls(void)
{
}


Fast_ls::~Fast_ls(void)
{
}
vector<Point> del2( Raw &phi,vector<Point> narrow ) 
{
	int m=phi.getXsize();
	int n=phi.getYsize();
	int l=phi.getZsize();
	//Raw ret2(m,n,l);
	vector<Point>ret2;
	//for (int i=0;i<m;i++)
	//{
	//	for(int j=0;j<n;j++)
	//	{
	//		for (int k=0;k<l;k++)
	//		{
	//			if (i+1< m && j+1< n &&k+1< l && i-1>=0 && j-1 >= 0 && k-1>= 0)
	//			{
	//				PIXTYPE value = (1.0/6.0)*(phi.get(i+1, j, k) + phi.get(i-1, j, k) + phi.get(i, j-1, k) + phi.get(i, j+1, k) + phi.get(i,j,k+1)+phi.get(i,j,k-1)- 6*(phi.get(i,j,k)));
	//				ret2.put(i, j,k, (PIXTYPE)value);
	//			}
	//			else 
	//			{
	//				ret2.put(i, j,k,0);
	//			}
	//		}

	//		
	//	}
	//}

	for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	{
		int i = it->xindex;
		int j = it->yindex;
		int k =  it->zindex;
		if (i+1< m && j+1< n &&k+1< l && i-1>=0 && j-1 >= 0 && k-1>= 0)
		{
			PIXTYPE value = (1.0/6.0)*(phi.get(i+1, j, k) + phi.get(i-1, j, k) + phi.get(i, j-1, k) + phi.get(i, j+1, k) + phi.get(i,j,k+1)+phi.get(i,j,k-1)- 6*(phi.get(i,j,k)));
			ret2.push_back(Point(i, j,k, (PIXTYPE)value));
		}
		else 
		{
			ret2.push_back(Point(i, j,k,0));
		}
	}


	return ret2;
}
//real in use  Wenfeng SONG 20140508
vector<Point> gradientxgc( Raw &g ,vector <Point> narrow) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	//Raw ret(g);
	int i,j,k,temp1,temp2;
	vector<Point> ret;
	//for(i=0;i<n;i++)
	//{
	//	for(j=0;j<m;j++)
	//	{
	//		for ( k=0;k < l;k++)
	//		{
	//			if(i>0)
	//				temp1=i-1;
	//			else
	//				temp1=0;
	//			if (i<n-1)
	//				temp2=i+1;
	//			else 
	//				temp2=n-1;
	//			ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
	//			//if (ret.get(i,j,k)!=0)
	//			//{
	//			//	cout<<"i="<<i<<",j="<<j<<",k="<<k<<ret.get(i,j,k)<<endl;
	//			//}
	//		}			
	//	}
	//}
	//
	for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	{
		int i=it->xindex;
		if(i>0)
			temp1=i-1;
		else
			temp1=0;
		if (i<n-1)
			temp2=i+1;
		else 
			temp2=n-1;
		ret.push_back(Point(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0));

	}
	return ret;
}

Raw& gradientxgc( Raw &g ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	Raw ret(g);
	int i,j,k,temp1,temp2;
	//vector<Point> ret;
	//for(i=0;i<n;i++)
	//{
	//	for(j=0;j<m;j++)
	//	{
	//		for ( k=0;k < l;k++)
	//		{
	//			if(i>0)
	//				temp1=i-1;
	//			else
	//				temp1=0;
	//			if (i<n-1)
	//				temp2=i+1;
	//			else 
	//				temp2=n-1;
	//			ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
	//			//if (ret.get(i,j,k)!=0)
	//			//{
	//			//	cout<<"i="<<i<<",j="<<j<<",k="<<k<<ret.get(i,j,k)<<endl;
	//			//}
	//		}			
	//	}
	//}
	//
	//for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	//{
	//	int i=it->xindex;
	//	if(i>0)
	//		temp1=i-1;
	//	else
	//		temp1=0;
	//	if (i<n-1)
	//		temp2=i+1;
	//	else 
	//		temp2=n-1;
	//	ret.push_back(Point(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0));

	//}
	for(int i=0; i<n; ++i)
	{
		for(int j = 0; j < m; ++j )
		{
			for(int k = 0; k < l; ++k)
			{
				ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
			}
		}
	}

	return ret;
}
//real in use swf 20140508
vector<Point> gradientygc( Raw & g,vector<Point> narrow ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	//Raw ret(n, m,l);
	vector<Point> ret;
	int i,j,k,temp1,temp2;
	//for(i=0;i<n;i++)
	//{
	//	for(j=0;j<m;j++)
	//	{
	//		for ( k = 0;k < l;k++)
	//		{
	//			if(j>0)
	//				temp1=j-1;
	//			else
	//				temp1=0;
	//			if (j<n-1)
	//				temp2=j+1;
	//			else 
	//				temp2=m-1;
	//			ret.put(i,j,k,0.5*(g.get(i,temp2,k)-g.get(i,temp1,k)));
	//		}			
	//	}
	//}
	
	for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	{
		int j=it->yindex;
		if(j>0)
			temp1=j-1;
		else
			temp1=0;
		if (j<n-1)
			temp2=j+1;
		else 
			temp2=m-1;
		ret.push_back(Point(i,j,k,0.5*(g.get(i,temp2,k)-g.get(i,temp1,k))));
	}
	return ret;
	
}
Raw& gradientygc( Raw &g ) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	Raw ret(g);
	int i,j,k,temp1,temp2;
	//vector<Point> ret;
	//for(i=0;i<n;i++)
	//{
	//	for(j=0;j<m;j++)
	//	{
	//		for ( k=0;k < l;k++)
	//		{
	//			if(i>0)
	//				temp1=i-1;
	//			else
	//				temp1=0;
	//			if (i<n-1)
	//				temp2=i+1;
	//			else 
	//				temp2=n-1;
	//			ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
	//			//if (ret.get(i,j,k)!=0)
	//			//{
	//			//	cout<<"i="<<i<<",j="<<j<<",k="<<k<<ret.get(i,j,k)<<endl;
	//			//}
	//		}			
	//	}
	//}
	//
	//for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	//{
	//	int i=it->xindex;
	//	if(i>0)
	//		temp1=i-1;
	//	else
	//		temp1=0;
	//	if (i<n-1)
	//		temp2=i+1;
	//	else 
	//		temp2=n-1;
	//	ret.push_back(Point(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0));

	//}
	for(int i=0; i<n; ++i)
	{
		for(int j = 0; j < m; ++j )
		{
			for(int k = 0; k < l; ++k)
			{
				ret.put(i,j,k,(g.get(temp2,j,k)-g.get(temp1,j,k))/2.0);
			}
		}
	}

	return ret;
}
// real in use swf 20140508
vector<Point> gradientzgc( Raw &g ,vector <Point> narrow) 
{
	int n=g.getXsize();
	int m=g.getYsize();
	int l=g.getZsize();
	//Raw ret(n, m,l);
	vector<Point> ret;
	int i,j,k,temp1,temp2;

	for (vector<Point>::iterator it = narrow.begin() ; it != narrow.end(); ++it)
	{
		int i=it->zindex;
		if(k>0)
			temp1=k-1;
		else
			temp1=0;
		if (k<l-1)
			temp2=k+1;
		else 
			temp2=l-1;
		ret.push_back(Point(i,j,k,0.5*(g.get(i,j,temp2)-g.get(i,j,temp1))));

	}
	return ret;
}
vector<Point> cos(Raw &x,vector<Point> narrow)
{
	vector <Point> ret;

	for (vector<Point>::iterator it=narrow.begin(); it !=narrow.end(); ++iter)
	{
		PIXTYPE val = cos(it->value);
		ret.push_back(Point(it->xindex,it->yindex,it->zindex,val));
	}
	return ret;
}
vector <Point> div(Raw &x,Raw &y,Raw &z)
{
	return (gradientxgc(x) += gradientygc(y) += gradientzgc(z)).set_shared(true);
}
void Fast_ls::array_surface(Raw *src)
{
	for (int i=0;i<src->getZsize();i++)
	{

		
	}
}
vector<Point> regFunction(vector<Point> narrow,double m,double n)
{

	//Raw ss(x,y,z);
	vector<Point> ss;
	PIXTYPE val=0;
	for (vector<Point>::iterator it =narrow.begin(); it !=narrow.end();++it)
	{
		val=it->value;
		if(val>=m && val<=n )
		{
			//ss->put(i,j,255);//unsigned char version
			ss.push_back(Point(i,j,k,1));
		}
		else if(it->value==m||it->value==n) 
		{
			ss.push_back(Point(i,j,k,0));
		}
		else 
		{
			ss.push_back(Point(i,j,k,0));
		}

	}
	

	return ss;
}
vector<Point>  Dirac( vector<Point> narrow, double sigma ) 
{
	//Raw ret(x);
	vector<Point> ret;
	PIXTYPE temp=((1.0/2.0)/sigma);
	//double temp2=(cos((2/sigma)*pi)+1)*temp;
	vector<Point> f= (cos((ret/sigma)*pi)+1)*temp;
	vector<Point> b = regFunction(ret, -sigma, sigma);
	//IShowImg(b);
	ret = f*b;
	//IShowImg(ret);

	return ret;
}
vector<Point> sin(vector<Point> src)
{
	vector <Point> ret;

	for (vector<Point>::iterator it=narrow.begin(); it !=narrow.end(); ++iter)
	{
		PIXTYPE val = sin(it->value);
		ret.push_back(Point(it->xindex,it->yindex,it->zindex,val));
	}
	return ret;
}
Raw distReg_p2( Raw  &phi ) 
{

		//int m=phi.getXsize();
		//int n=phi.getYsize();
		//Raw phi_x(gradientxgc(phi));
		//Raw phi_y(gradientygc(phi));
		//Raw phi_z(gradientzgc(phi));
		//Raw s(ImageFSqrt(phi_x,phi_y,phi_z));
		//Raw &a=regFunction(s,0,1);
		//Raw &b=regFunction(s,1,10);//need to be changed.
		//Raw &ps=a*(sin((s*2.0*pi)))/(2.0*pi)+b*(s-1);
		//Raw &dps=((regFunction(ps,0,0)*(-1)+1)*ps+regFunction(ps,0,0))/((regFunction(s,0,0)*(-1)+1)*s+regFunction(s,0,0));
		//Raw &div2=div(dps*phi_x-phi_x,dps*phi_y-phi_y,dps*phi_z-phi_z);
		//Raw &f=div2-del2(phi)*4.0;
		//return f;

}
void Fast_ls::initialg(Raw &g)
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
Raw Fast_ls::minimal_surface(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction )
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



void Fast_ls::outerwall(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
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

void Fast_ls::outerwallauto(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
{

	this->initialg(g);
	for (int i=0; i < 30; i++)
	{
		phi=this->minimal_surface(phi,g,lambda,mu,alfa,epsilon,timestep,1,potentialFunction);
		cout << "outer wall auto stop iter = " << i <<endl;
	}

}


void Fast_ls::couple(Raw &phi,Raw &g,double lambda,double mu,double alfa,float epsilon,int timestep,int iter,char *potentialFunction)
{

}




void Fast_ls::NeumannBoundCond( Raw &img )
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

vector<Point> Fast_ls::ImageFSqrt( vector<Point> phi_x, vector<Point> phi_y,vector<Point> phi_z )
{
	vector<Point> ret;
	vector<Point>::iterator it1=phi_y.begin();  
	vector<Point>::iterator it2=phi_z.begin(); 
	for (vector<Point>::iterator it=phi_x.begin(); it!=phi_x.end(); ++it)
	{
		++it1;
		++it2;
		ret.push_back(Point(it->xindex,it->yindex,it->zindex,sqrt((double)it->value*it->value+it1->value*it1->value+it2->value*it2->value)));
	}
	return ret;
}
// find the zero crossings
vector<Point> Fast_ls::findNarrow(Raw *src)
{
	vector<Point> narrowband;
	vector<Point> zero;
	for (int i = 0; i < src->getXsize(); i++ )
	{
		for (int j=0; j <src->getYsize(); j++)
		{
			for ( int k = 0; k<src->getZsize(); k++ )
			{
				PIXTYPE val =src->get(i,j,k);
				PIXTYPE val_x1=src->get(i-1,j,k);
				PIXTYPE val_x2=src->get(i+1,j,k);
				PIXTYPE val_y1=src->get(i,j-1,k);
				PIXTYPE val_y2=src->get(i,j+1,k);
				PIXTYPE val_z1=src->get(i,j,k-1);
				PIXTYPE val_z2=src->get(i,j,k+1);

				int r=1;
				//
				if (val_x1*val_x2<0||val_y1*val_y2<0|| val_z1*val_z2<0)
				{
					zero.push_back(Point(i,j,k,val));
					//check the (2*r+1)*(2*r+1)neighbor hood
					for (int l = -r+i; l < r+i; l++)
					{
						for (int m = -r+j; m < r+j; m++)
						{
							for (int n = -r+k; n < r+k; n++)
							{
								PIXTYPE neighbour= src->get(l,m,n);
								narrowband.push_back(Point(l,m,n,neighbour));
							}
						}
					}
				}

			}
		}
	}
	return narrowband;
}
//src is the full domain of the phi first time is initial
Raw Fast_ls::minimal_surface(Raw &phi,vector <Point> narrow,Raw &g,double lambda,double mu,
	double alfa,float epsilon,int timestep,int iter,char *potentialFunction )
{

		this->narrow= findNarrow(&phi);
	int m=g.getXsize();
	int n=g.getYsize();
	int l=g.getZsize();
	vector<Point> vx=gradientxgc(g,narrow);
	vector<Point> vy=gradientygc(g,narrow);
	vector<Point> vz=gradientzgc(g,narrow);
	for(int i=0;i<iter;i++)
	{
		float smallNumber=1e-10;
		
		NeumannBoundCond(phi);
		vector<Point> phi_x = gradientxgc(phi,narrow);
		vector<Point> phi_y = gradientygc(phi,narrow);
		vector<Point> phi_z = gradientzgc(phi,narrow);
		vector<Point> s = ImageFSqrt(phi_x, phi_y, phi_z) += smallNumber;

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