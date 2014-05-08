#include "vol_math_Raw3D_Independt.h"
#include "vol_math_RawImage.h"
#include<iostream>
#include <queue>
#include <vector>
#include <queue>
#include"statistics.h"
using namespace std;
/*
1\delete the tree 
2\compute the length of the center line
3\divide the region 
*/
struct Point
{
	long long x;
	long long y;
	long long z;
	PIXTYPE value;
	
	Point(long long x,long long y,long long z,PIXTYPE value)
	{
		this->x=x;
		this->y=y;
		this->z =z;
		this->value=value;
	}
	Point ()
	{

	}
};
PIXTYPE  biglittleedianv2(PIXTYPE  val)
{

	PIXTYPE *q=&val;
	unsigned char * p= (unsigned char*) q;
	std::swap(p[0],p[3]);
	std::swap(p[1],p[2]);
	return val;
}
//vector <Point> SeedGrow(Point seed)
//{
//
//}
//deserve for graph storage structure
typedef struct edge
{
	int adjvex;
	int weight;
	struct edge *next;

}ELink;
typedef struct ver
{
	Point vertex;
	ELink *link;
}VLink  ;
class DivideRegion
{

public:
	queue<Point> q;
	Raw *skeletondata;
	vector<Point> center;
	long long  count;
	DivideRegion(queue<Point> q,Raw *s,vector<Point> center)
	{
		this->q=q;
		this->skeletondata=new Raw(*s);
		this->center=center;
		count=0;
	}
	void checkadj(Point firstseed,int ff)
	{
		
		for (int l = -1; l <=1 ;l++)
		{
			for (int m =-1; m<= 1; m++)
			{
				for(int n=-1; n<= 1;n++)
				{
					int posx = firstseed.x+l;
					int posy = firstseed.y+m;
					int posz = firstseed.z+n;
					if (posx<skeletondata->getXsize()&&posy<skeletondata->getYsize()&&posz<skeletondata->getZsize()
						&&skeletondata->get(posx,posy,posz)!=0&&skeletondata->get(posx,posy,posz)!=ff)
					{
						//firstseed.value=0;
						Point pos(posx,posy,posz,ff);
						skeletondata->put(posx,posy,posz,ff);
						center.push_back(pos);
						this->count++;
						q.push(pos);
					}


				}
			}
		}
	}
	//delete the noise points whose neighbor is less than one 
	int deletenoise(Point firstseed,int size,int threshold)
	{
		int neighbor=0;
		for (int l = -size; l <=size ;l++)
		{
			for (int m =-size; m<= size; m++)
			{
				for(int n=-size; n<= size;n++)
				{
					int posx = firstseed.x+l;
					int posy = firstseed.y+m;
					int posz = firstseed.z+n;
					if (posx<skeletondata->getXsize()&&posy<skeletondata->getYsize()&&posz<skeletondata->getZsize()
						&&skeletondata->get(posx,posy,posz)!=0)
					{

						neighbor++;
						//firstseed.value=0;
						/*Point pos(posx,posy,posz,ff);
						skeletondata->put(posx,posy,posz,ff);
						center.push_back(pos);
						this->count++;
						q.push(pos);*/
					}
					

				}
			}
		}
		if (neighbor < threshold)
		{
			skeletondata->put(firstseed.x,firstseed.y,firstseed.z,0);
			return 1;
		}
		return 2;
	}
	void wirtecenterline(int x,int y ,int z)
	{
		PIXTYPE * data= new PIXTYPE[x*y*z];
		memset(data,0,sizeof(PIXTYPE)*x*y*z);
		Raw *centerdata=new Raw(x,y,z,data);

		for (std::vector<Point>::iterator it = center.begin() ; it != center.end(); ++it)
		{
			centerdata->put(it->x,it->y,it->z,255);
		}
		RawImage test;
		test.writeImagecolon(*centerdata);
	}
	void recoverskeletonata()
	{
		for (int i=0;i< center.size();i++)
		{
			skeletondata->put(center[i].x,center[i].y,center[i].z,255);
		}

	}
	void markneighbor(Point p)
	{
		int i=p.x,j=p.y,k=p.z;
			int ff=200;
		for (int l=i-1 ; l <= i+1; l++)
		{
			for (int m = j-1; m <= j+1; m++)
			{
				for (int n = k-1; n <= k+1; n++)
				{

					if (skeletondata->get(l,m,n)==3)
					{
						//ssscount++;
						//Point point(l,m,n,skeletondata->get(l,m,n));
						skeletondata->put(l,m,n,ff);
					
					}
				}
			}
		}

	}

	//void checkadjminusone(Point firstseed,int ff)
	//{
	//	
	//	for (int l = -1; l <=1 ;l++)
	//	{
	//		for (int m =-1; m<= 1; m++)
	//		{
	//			for(int n=-1; n<= 1;n++)
	//			{
	//				int posx = firstseed.x+l;
	//				int posy = firstseed.y+m;
	//				int posz = firstseed.z+n;
	//				if (posx<skeletondata->getXsize()&&posy<skeletondata->getYsize()&&posz<skeletondata->getZsize()
	//					&&skeletondata->get(posx,posy,posz)!=ff)
	//				{
	//					//firstseed.value=0;
	//					Point pos(posx,posy,posz,ff);
	//					skeletondata->put(posx,posy,posz,ff);
	//					center.push_back(pos);
	//					this->count++;
	//					q.push(pos);
	//				}


	//			}
	//		}
	//	}
	//}

	//the BFS for store the skeleton in order add a queue store different vertexes
	void  putskletoninorder()
	{
		//int i=0;
		//vector <Point> center;
		int i,j,k;
		Point firstseed;
		//queue<Point> q;
		//int count=0;
		bool flag =true;
		int ff=0;
		for ( i =0;i< skeletondata->getXsize()&&flag ==true ;i++)
		{
			for ( j=0;j<skeletondata->getYsize()&&flag ==true;j++)
			{
				for ( k=0;k<skeletondata->getZsize()&&flag ==true;k++)
				{
					
					if (skeletondata->get(i,j,k)!=0)
					{
						count++;
						firstseed.x=(long long)i;
						firstseed.y=(long long)j;
						firstseed.z=(long long)k;
						firstseed.value=0;
						flag =false;
						
					}
				}
			}
		}
		center.push_back(firstseed);
		checkadj(firstseed,ff);
		while(!q.empty())
		{
			Point pos=q.front();
			skeletondata->put(pos.x,pos.y,pos.z,0);
			checkadj(pos,ff);
			q.pop();
		}
		// swf change if to while
		while (center.size()< count-100)
		{
			bool flag =true;
			for ( i =0;i< skeletondata->getXsize()&&flag ==true ;i++)
			{
				for ( j=0;j<skeletondata->getYsize()&&flag ==true;j++)
				{
					for ( k=0;k<skeletondata->getZsize()&&flag ==true;k++)
					{

						if (skeletondata->get(i,j,k)!=0)
						{
							count++;
							firstseed.x=(long long)i;
							firstseed.y=(long long)j;
							firstseed.z=(long long)k;
							firstseed.value=0;
							flag =false;

						}
					}
				}
			}
			
			center.push_back(firstseed);
			checkadj(firstseed,ff);
			while(!q.empty())
			{
				Point pos=q.front();
				skeletondata->put(pos.x,pos.y,pos.z,0);
				checkadj(pos,ff);
				q.pop();
			}
		}

	}
	// given the first seed
	void  putskletoninorderpublic(Point firstseed)
	{
		int i,j,k;
		int ff=0;
		center.clear();
		center.push_back(firstseed);
		checkadj(firstseed,ff);
		while(!q.empty())
		{
			Point pos=q.front();
			skeletondata->put(pos.x,pos.y,pos.z,ff);
			checkadj(pos,ff);
			q.pop();
		}
		//swf change if to while
		while (center.size()< count/20)
		{
			bool flag =true;
			for ( i =0;i< skeletondata->getXsize()&&flag ==true ;i++)
			{
				for ( j=0;j<skeletondata->getYsize()&&flag ==true;j++)
				{
					for ( k=0;k<skeletondata->getZsize()&&flag ==true;k++)
					{

						if (skeletondata->get(i,j,k)!=0)
						{

							firstseed.x=(long long)i;
							firstseed.y=(long long)j;
							firstseed.z=(long long)k;
							firstseed.value=ff;
							flag =false;

						}
					}
				}
			}

			center.push_back(firstseed);
			checkadj(firstseed,ff);
			while(!q.empty())
			{
				Point pos=q.front();
				skeletondata->put(pos.x,pos.y,pos.z,ff);
				checkadj(pos,ff);
				q.pop();
			}
		}

	}
	// solution 1,count the little tree
	void  putskletoninorders1()
	{
		//int i=0;
		//vector <Point> center;
		int i,j,k;
		int ff=-1;
		Point firstseed;
		//queue<Point> q;
		//int count=0;
		bool flag =true;
		for (int iter =0;iter< 2;iter++)
		{
	
				for ( i =0;i< skeletondata->getXsize() ;i++)
				{
					for ( j=0;j<skeletondata->getYsize();j++)
					{
						for ( k=0;k<skeletondata->getZsize();k++)
						{

							if (skeletondata->get(i,j,k)>0)
							{
						
								Point *pos=new Point(i,j,k,skeletondata->get(i,j,k));
								int ifdelete=deletenoise(*pos,1,2);
								if (ifdelete==2)
								{
									count++;
								}
						
						
							}//end k
						}//end j
					}
				}//end i
		}
		for ( i =0;i< skeletondata->getXsize()&&flag ==true ;i++)
		{
			for ( j=0;j<skeletondata->getYsize()&&flag ==true;j++)
			{
				for ( k=0;k<skeletondata->getZsize()&&flag ==true;k++)
				{

					if (skeletondata->get(i,j,k)>0)
					{
						
						firstseed.x=(long long)i;
						firstseed.y=(long long)j;
						firstseed.z=(long long)k;
						firstseed.value=ff;
						flag =false;

					}
				}
			}
		}
		center.push_back(firstseed);
		skeletondata->put(firstseed.x,firstseed.y,firstseed.z,0);
		checkadj(firstseed,ff);
		while(!q.empty())
		{
			Point pos=q.front();
			skeletondata->put(pos.x,pos.y,pos.z,ff);
			checkadj(pos,ff);
			q.pop();
		}
		//swf 20140418 delete the while ,just for find the end point
		//// swf change if to while
		while (center.size()< count/100)
		{
			bool flag =true;
			for ( i =0;i< skeletondata->getXsize()&&flag ==true ;i++)
			{
				for ( j=0;j<skeletondata->getYsize()&&flag ==true;j++)
				{
					for ( k=0;k<skeletondata->getZsize()&&flag ==true;k++)
					{

						if (skeletondata->get(i,j,k) >0)
						{
						
							firstseed.x=(long long)i;
							firstseed.y=(long long)j;
							firstseed.z=(long long)k;
							firstseed.value=ff;
							flag =false;

						}
					}
				}
			}
			center.push_back(firstseed);
			checkadj(firstseed,ff);
			while(!q.empty())
			{
				Point pos=q.front();
				skeletondata->put(pos.x,pos.y,pos.z,ff);
				checkadj(pos,ff);
				q.pop();
			}
		}//end while
		putskletoninorderpublic(center.front());
		for (int i=0;i< center.size();i++)
		{
			skeletondata->put(center[i].x,center[i].y,center[i].z,200);
		}
	}
	//use cos tell which line to go
	// input:skeleton data ,first start point
	// output:center
	double pointcos(Point a,Point b)
	{
		int ff=1;
		//Point ab(a.x-b.x,a.y-b.y,(a.z-b.z),ff);
		double ret=(double)(a.x*b.x+a.y*b.y+a.z*b.z)/(double)(sqrt((double)a.x*a.x+a.y*a.y+a.z*a.z)*sqrt((double)(b.x*b.x+b.y*b.y+b.z*b.z)));
		return ret;
	}
	Point makevector(Point a,Point b)
	{
		int ff = 3;
		Point ret(a.x-b.x,a.y-b.y,a.z-b.z,ff);
		return ret;
	}
	void putneighborzero()
	{

	}
	Point* coutadjcos(Point pos,Point pre)
	{
		int ff=3;
	//	Point pre =center.back();
		Point centerdirection=makevector(pre,pos);
		int count =0;
		Point *next=new Point();
		int i=pos.x;
		int j=pos.y;
		int k= pos.z;
		
		vector<Point> adj;
		double max=0;
		for (int l=i-1 ; l <= i+1; l++)
		{
			for (int m = j-1; m <= j+1; m++)
			{
				for (int n = k-1; n <= k+1; n++)
				{
					//if (l==52&&m==358&&n==221)
					//{
					//	cout <<skeletondata->get(l,m,n)<<endl;
					//}
					//if (l==52&&m==356&&n==223)
					//{
					//	cout <<skeletondata->get(l,m,n)<<endl;
					//}
					bool f1= skeletondata->get(l,m,n)==3;
					bool f2= !(l==pos.x&&m==pos.y&&n==pos.z);
					if (f1 &&f2 )
					{
						Point nextadj(l,m,n,3);
						skeletondata->put(l,m,n,1);//cos
						nextadj.x=l;
						nextadj.y=m;
						nextadj.z=n;
						Point treedir=makevector(pos,nextadj);
						double val=abs(pointcos(centerdirection,treedir));
						if (max <= val)
						{
							max = val;
							next->x=l;
							next->y=m;
							next->z=n;
							next->value=ff;
							
						}
						
					
					}
				}
			}
		}
		return next;
	}
	vector<Point> coutadj(int i,int j,int k)
	{
		int count =0;
		int ff = 3;
		vector<Point> adj;
		for (int l=i-1 ; l <= i+1; l++)
		{
			for (int m = j-1; m <= j+1; m++)
			{
				for (int n = k-1; n <= k+1; n++)
				{
					
					if (skeletondata->get(l,m,n)==200)
					{
						count++;
						Point point(l,m,n,skeletondata->get(l,m,n));
						skeletondata->put(l,m,n,ff);
						adj.push_back(point);
					}
				}
			}
		}
		return adj;
	}
	int  countall(Point pos,int count)
	{
		bool flag=true;
		int threshold=100;
		while (flag)
		{
			int ff=1;
			vector<Point>  adj=coutadj(pos.x,pos.y,pos.z);
			for (int i = 0; i < adj.size() && count <= threshold; i++)
			{
				countall(adj[i],count+=1);
			}
			flag =false;
		}
		return count;
	}
	// look up which point adj is white without change flag just record
	vector<Point> justcountadj(int i,int j,int k)
	{
		int count =0;
		vector<Point> adj;
		
		for (int l=i-1 ; l <= i+1; l++)
		{
			for (int m = j-1; m <= j+1; m++)
			{
				for (int n = k-1; n <= k+1; n++)
				{

					if (skeletondata->get(l,m,n)!=0)
					{
						count++;
						Point point(l,m,n,skeletondata->get(l,m,n));
						adj.push_back(point);
					}
				}
			}
		}
		return adj;
	}
	//void putskeletoninorders2(Point start,Point end);
	void putskeletoninorders2(Point start ,Point end)
	{
		vector<Point> adj;
		vector<Point> centerbak=center;
		int bakindex=0;
		center.clear();
		int total=0;
		//for (int i = 1; i < skeletondata->getXsize()-1; i++)
		//{
		//	for ( int j = 1; j< skeletondata->getYsize()-1 ; j++ )
		//	{
		//		for(int k=1 ; k <skeletondata->getZsize()-1 ; k++)
		//		{
		Point pos=start;
		bool flag =true;
		while (flag)
		{
		

					
					adj= justcountadj(pos.x,pos.y,pos.z);
					if ((adj.size())>=3)
					{
						Point maxpoint;
						int max=0;
						for (int ll=0;  ll< adj.size(); ++ll)
						{
							if (adj[ll].x==pos.x&&adj[ll].y==pos.y&&adj[ll].z!=pos.z)
							{

							}
							else
							{
									pos =adj[ll];
									int count=0;
									int allregion=countall(pos,count);
									if (allregion > max)
									{
										max =allregion;
										maxpoint=pos;
										pos=maxpoint;
									}
							}
						}
						center.push_back(pos);
						if (pos.x==end.x&&pos.y==end.y&&pos.z==end.z)
						{
							flag=false;
							break;
						}
						
					}
					else if(adj.size()==2) 
					{
						pos=adj[0];
						center.push_back(pos);
						if (pos.x==end.x&&pos.y==end.y&&pos.z==end.z)
						{
							flag=false;
							break;
						}
					}
					else if ((total=center.size())<count/20)
					{
						pos=centerbak[++bakindex];
					}
					else flag=false;

						
				
			}
		//		}//end k
		//	}//end j
		//}//end i
		
	}
	//void putskeletoninorders2(Point start,Point end);
	void putskeletoninorders3(Point start ,Point end)
	{
		vector<Point> newcenter;
		vector<Point> adj;
		bool flag =true;
		Point pos=start;
		while (flag)
		{
			if (newcenter.size()==358)
			{
				cout <<"10"<<endl;
			}
			vector<Point> val=coutadj(pos.x,pos.y,pos.z);//flag==2;
			if (val.size()==1)
			{
				newcenter.push_back(val.back());
				skeletondata->put(val.back().x,val.back().y,val.back().z,0);
				if (newcenter.back().x == end.x && newcenter.back().y==end.y&& newcenter.back().z==end.z)
				{
					flag=false;
					break;
				}
				pos=val.back();
			}
			else if(val.size()!=0)
				{
					
					if (pos.x==start.x&&pos.y==start.y&&pos.z==start.z)
					{
						//newcenter.push_back(pos);
						//pos=center[2];
						//newcenter.push_back(pos);
						//skeletondata->put(center[1].x,center[1].y,center[1].z,0);
						//skeletondata->put(pos.x,pos.y,pos.z,0);
						int size=val.size();
						newcenter.push_back(pos);
						pos=*(coutadjcos(pos,center[size]));
					}
					else 
					{
						int size=newcenter.size();
						if (pos.x==newcenter.back().x&&pos.y==newcenter.back().y&&pos.z==newcenter.back().z)
						{
							pos=*(coutadjcos(pos,newcenter[size-4]));//change 2 to 4
						}
						else 
						{
							newcenter.push_back(pos);
							if (newcenter.size()>10)
							{
								pos=*(coutadjcos(pos,newcenter[size-6]));
							}
							else pos=*(coutadjcos(pos,newcenter[size-1]));

						}
						
						if (newcenter.back().x == end.x && newcenter.back().y==end.y&& newcenter.back().z==end.z)
						{
							flag=false;
							break;
						}
					}
					
			}
			// no points neighbour
			else 
			{
				// point is enough
				if (newcenter.size()> center.size()*5/6 )
				{
					flag = false;
				}
				else 
				{
					for (int i=0;i<center.size();i++)
					{
						if (center[i].x==pos.x&&center[i].y==pos.y&&center[i].z==pos.z)
						{
							pos=center[i+1];
							newcenter.push_back(pos);
							skeletondata->put(pos.x,pos.y,pos.z,0);
							break;
						}
					}
				}
			}

		}
		center.clear();
		center=newcenter;
	}
	void putskeletoninorders3s2(Point start ,Point end)
	{
		vector<Point> newcenter;
		vector<Point> adj;
		bool flag =true;
		Point pos=start;
		vector <Point> nodeseperateflag;
		Point nodechoosed;
		while (flag)
		{
			if (newcenter.size()==358)
			{
				cout <<"10"<<endl;
			}
			vector<Point> val=coutadj(pos.x,pos.y,pos.z);//flag==2;
			if (val.size()==1)
			{
				newcenter.push_back(val.back());
				skeletondata->put(val.back().x,val.back().y,val.back().z,0);
				if (newcenter.back().x == end.x && newcenter.back().y==end.y&& newcenter.back().z==end.z)
				{
					flag=false;
					break;
				}
				pos=val.back();
			}
			else if(val.size()!=0)
			{
				
				nodeseperateflag=val;
				if (pos.x==start.x&&pos.y==start.y&&pos.z==start.z)
				{
					//newcenter.push_back(pos);
					//pos=center[2];
					//newcenter.push_back(pos);
					//skeletondata->put(center[1].x,center[1].y,center[1].z,0);
					//skeletondata->put(pos.x,pos.y,pos.z,0);
					int size=val.size();
					newcenter.push_back(pos);
					pos=*(coutadjcos(pos,center[size]));
				}
				else 
				{
					int size=newcenter.size();
					if (pos.x==newcenter.back().x&&pos.y==newcenter.back().y&&pos.z==newcenter.back().z)
					{
						if (size>5)
						{
							pos=*(coutadjcos(pos,newcenter[size-4]));//change 2 to 4
						} 
						else
						{
							pos=*(coutadjcos(pos,newcenter[size-2]));//change 2 to 4
						}
						
					}
					else 
					{
						newcenter.push_back(pos);
						if (newcenter.size()>10)
						{
							pos=*(coutadjcos(pos,newcenter[size-6]));
						}
						else pos=*(coutadjcos(pos,newcenter[size-1]));

					}

					if (newcenter.back().x == end.x && newcenter.back().y==end.y&& newcenter.back().z==end.z)
					{
						flag=false;
						break;
					}
				}

				nodechoosed=pos;
			}
			// no points neighbour
			else 
			{
				// point is enough
				if (newcenter.size()> center.size()*5/6 )
				{
					flag = false;
				}
				else 
				{
					bool fflag=true;
					for (int i=0;i<nodeseperateflag.size();i++)
					{
						
						if (!(nodeseperateflag[i].x==nodechoosed.x&&nodeseperateflag[i].y==nodechoosed.y&&nodeseperateflag[i].z==nodechoosed.z))
						{
							//the error tree come out from center
							for ( int i=newcenter.size()-1; i>0 && flag; i--)
							{
								flag=skeletondata->get(newcenter[i].x,newcenter[i].y,newcenter[i].z)==0;
								
								skeletondata->put(newcenter[i].x,newcenter[i].y,newcenter[i].z,200);
								
								if (newcenter[i].x=nodechoosed.x&&newcenter[i].y==nodechoosed.y&&newcenter[i].z==nodechoosed.z)
								{
									//newcenter.pop_back();
									//skeletondata->put(newcenter[i].x,newcenter[i].y,newcenter[i].z,255);

									newcenter.pop_back();
									break;
								}
								else
									newcenter.pop_back();
							}
							if (!flag)//
							{
								pos = nodeseperateflag[i];
								newcenter.push_back(pos);
								nodechoosed=pos;
								markneighbor(pos);
								skeletondata->put(pos.x,pos.y,pos.z,0);
								break;

							}
							else 
							{
								for (int i=0;i<center.size();i++)
								{
									if (center[i].x==pos.x&&center[i].y==pos.y&&center[i].z==pos.z)
									{
										pos=center[i+1];
										newcenter.push_back(pos);
										skeletondata->put(pos.x,pos.y,pos.z,0);
										break;
									}
								}
							}
						}//if
					}//for
				}
			}

		}
		center.clear();
		center=newcenter;
	}
	//compute the HU statistics
	void DivideRegionv1(Raw *colonskeleton,Raw *data)
	{
		//skeletondata=colonskeleton;
		int count=0;
		vector<Point> skeleton;

		for (int i = 0; i < colonskeleton->getXsize(); i++)
		{
			for (int j=0 ; j < colonskeleton->getYsize() ; j++)
			{
				for (int k=0; k<colonskeleton->getZsize();k++)
				{
					//big --little edian
					PIXTYPE val =biglittleedianv2(colonskeleton->get(i,j,k));
					if (val!= 0)
					{
						Point point(i,j,k,255);
						skeleton.push_back(point);
						count++;
					}
				}
			}


		}
		short histgram[4][3000]={0};
		for (int i = 1; i < skeleton.size(); i++)
		{
			int box=10;
			if (i <= count/4 &&i >=1)
			{
				Point pos=skeleton[i];

				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

				
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[0][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}

			}
			else if (i<= 2*count/4)
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[1][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else if (i<= 3*count/4)
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[2][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[3][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			}



		}
		for (int i=0;i<4;i++)
		{
			//for (int j = 0; j < 255; j++)
			//{
			ofstream os("histgram.txt",ios::app);
			if (os)
			{
				for (int j=0; j <1000; j++)
				{
					os << histgram[i][j] <<" ";

				}
				os << endl;
			} 
			/*FILE *p=fopen("histgram.txt","a+");
			fwrite((void *)histgram[i],sizeof(short),3000,p);*/
			//}
		}

	}
	void DivideRegionv2(Raw *colonskeleton,Raw *data)
	{
		//skeletondata=colonskeleton;
		int count=center.size();
		//vector<Point> skeleton;
		//for (int i = 0; i < colonskeleton->getXsize(); i++)
		//{
		//	for (int j=0 ; j < colonskeleton->getYsize() ; j++)
		//	{
		//		for (int k=0; k<colonskeleton->getZsize();k++)
		//		{
		//			//big --little edian
		//			PIXTYPE val =biglittleedianv2(colonskeleton->get(i,j,k));
		//			if (val!= 0)
		//			{
		//				Point point(i,j,k,0);
		//				center.push_back(point);
		//			
		//				//count++;
		//			}
		//		}
		//	}
		//	

		//}
		short histgram[4][3000]={0};
		for (int i = 1; i < center.size(); i++)
		{
			int box=10;
			if (i <= count/4 &&i >=1)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[0][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}

			}
			else if (i<= 2*count/4)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

				
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[1][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else if (i<= 3*count/4)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

			
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[2][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
			}
			} 
			else
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{

							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[3][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			}



		}
		for (int i=0;i<4;i++)
		{
			//for (int j = 0; j < 255; j++)
			//{
			ofstream os("histgram.txt",ios::app);
			if (os)
			{
				for (int j=0; j <1000; j++)
				{
					os << histgram[i][j] <<" ";

				}
				os << endl;
			} 
			/*FILE *p=fopen("histgram.txt","a+");
			fwrite((void *)histgram[i],sizeof(short),3000,p);*/
			//}
		}

	}
	void ComputeLength()
	{

	}
	void  HistHU(Raw *thickness,Raw *origion,Raw *hu)
	{

		for (int i=0;i<thickness->size();i++)
		{
			if (thickness->getXYZ(i)!=0)
			{
				hu->putXYZ(i,origion->getXYZ(i));
			}
		}

	}
//compute the divide region thickness,base on the origion data order 
	void DivideRegionthickness(Raw *colonskeleton,Raw *data)
	{
		int count=0;

		vector<Point> skeleton;
		for (int i = 0; i < colonskeleton->getXsize(); i++)
		{
			for (int j = 0; j < colonskeleton->getYsize() ; j++)
			{
				for (int k = 0; k < colonskeleton->getZsize(); k++)
				{
					if (colonskeleton->get(i,j,k)!= 0)
					{
						Point point(i,j,k,255);
						skeleton.push_back(point);
						count++;
					}
				}
			}


		}
		int histgram[4][20]={0};
		for (int i = 1; i < skeleton.size(); i++)
		{
			int box=10;
			if (i <= count/4)
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

				
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							PIXTYPE val =(data->get(l,m,n)) ;
							if (val !=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[0][(int)val]++;
								}
							}
						}
					}
				}
				}

			}
			else if (i<= 2*count/4)
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[1][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else if (i<= 3*count/4)
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[2][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else
			{
				Point pos=skeleton[i];
				Point prepos=skeleton[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[3][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			}



		}
		for (int i = 0; i < 4; i++)
		{
			//for (int j = 0; j < 255; j++)
			//{
			ofstream os("thicknesshistgram.txt",ios::app);
			for (int j=0;j<20;j++)
			{
				os<<histgram[i][j]<<" ";
			}
			os <<endl;
			/*FILE *p=fopen("thicknesshistgram.txt","a+");
			fwrite((void *)histgram[i],sizeof(int),20,p);*/
			//}
			//1\thickness 4
			// 2
			//\thickness new3 = num3-2*num4
			// \thickness new2 = num2-2*num4-2*new3
			// \thickneww new1= num1-2*new2-new3

		}

	}
	//store the skeleton in a line by the putskeleton in order function score the result 
	// in skeletondata which is the member of DivideRegion
	void DivideRegionthicknessv2(Raw *colonskeleton,Raw *data)
	{
		//skeletondata=colonskeleton;
		int count=center.size();
		int histgram[4][20]={0};
		int box=10;
		for (int i = 1; i < center.size(); i++)
		{

			if (i <= count/4)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							PIXTYPE val = data->get(l,m,n);
							if (val !=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[0][(int)val]++;
								}
							}
						}
					}
				}
				}
			}
			else if (i<= 2*count/4)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

				
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[1][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else if (i<= 3*count/4)
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{

				
				for (int l=pos.x-box; l< pos.x+box ;l++)
				{
					for (int m= pos.y-box; m< pos.y+box ; m++)
					{
						for (int n= pos.z-box ; n< pos.z+box; n++)
						{
							if (data->get(l,m,n)!=0)
							{
								bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
								if (flag)
								{
									histgram[2][(int)data->get(l,m,n)]++;
								}
							}
						}
					}
				}
				}
			} 
			else
			{
				Point pos=center[i];
				Point prepos=center[i-1];
				if (pos.x-box > 0 && pos.x + box <skeletondata->getXsize()&& 
					pos.y-box > 0 && pos.y + box <skeletondata->getYsize()&&
					pos.z-box > 0 && pos.z + box <skeletondata->getZsize())
				{
					for (int l=pos.x-box; l< pos.x+box ;l++)
					{
						for (int m= pos.y-box; m< pos.y+box ; m++)
						{
							for (int n= pos.z-box ; n< pos.z+box; n++)
							{
								if (data->get(l,m,n)!=0)
								{
									bool flag=(l-pos.x)*(pos.x-prepos.x)+(m-pos.y)*(pos.y-prepos.y)*(pos.z-prepos.z)==0;
									if (flag)
									{
										histgram[3][(int)data->get(l,m,n)]++;
									}
								}
							}
						}
					}
				}

			}



		}
		for (int i = 0; i < 4; i++)
		{
			//for (int j = 0; j < 255; j++)
			//{
			ofstream os("thicknesshistgram.txt",ios::app);
			for (int j=0;j<20;j++)
			{
				os<<histgram[i][j]<<" ";
			}
			os <<endl;
			/*FILE *p=fopen("thicknesshistgram.txt","a+");
			fwrite((void *)histgram[i],sizeof(int),20,p);*/
			//}
			//1\thickness 4
			// 2
			//\thickness new3 = num3-2*num4
			// \thickness new2 = num2-2*num4-2*new3
			// \thickneww new1= num1-2*new2-new3

		}

	}
};