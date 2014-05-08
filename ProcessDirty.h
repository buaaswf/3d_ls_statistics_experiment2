#include "vol_math_Raw3D_Independt.h"
#include "vol_math_RawImage.h"
#include<iostream>
#include <queue>
#include <vector>
#include <queue>
void processDirty(Raw *dirty)
{
	for (int i=0;i<dirty->getXsize(); i++)
	{
		for (int j=2; j <dirty->getYsize()-1; j++)
		{
			for (int k=100 ; k<dirty->getZsize()-62; k++)
			{
				PIXTYPE val=dirty->get(i,j,k);
				PIXTYPE val2= dirty->get(i,j-2,k);
				PIXTYPE val3= dirty->get(i,j-1,k);
				PIXTYPE val4 = dirty->get(i,j+1,k);
				if ( (val> 300 && val2< -600 ) )
				{
					
					dirty->put(i,j,k,val2);
				}
				else if ((val> 300 && val3< -600 ))
				{
					dirty->put(i,j,k,val3);
				}
				else if (val4-val3>1000&& val3<-500)
				{
					dirty->put(i,j,k,val3);
				}
				if (dirty->get(i,j,k)<-600)
				{
					dirty->put(i,j,k,-937);
				}
			}//..end for k
		}//..end for j
	}//..end for i

	for (int i=1;i<dirty->getXsize()-1; i++)
	{
		for (int j=2; j <dirty->getYsize()-1; j++)
		{
			for (int k=100 ; k<dirty->getZsize()-62; k++)
			{
				PIXTYPE val=dirty->get(i,j,k);
				PIXTYPE val2= dirty->get(i+1,j,k);
				PIXTYPE val3= dirty->get(i-1,j,k);
				PIXTYPE val4=dirty->get(i,j+1,k);
				PIXTYPE val5=dirty->get(i,j-1,k);
				PIXTYPE val6=dirty->get(i,j,k+1);
				PIXTYPE val7=dirty->get(i,j,k-1);
				if ( (val>500&&val3<-500 && val2<-500) )
				{

					dirty->put(i,j,k,-937);
				}
				else if ((val>500&& val4<-600) ||(val>500&& val5<-600))
				{
					dirty->put(i,j,k,-937);
				}
				else if ((val>500&& val6<-600) ||(val>500&& val7<-600)||(val>300&&abs(val6-val7)>1000)&&val7 <-500)
				{
					dirty->put(i,j,k,-937);
				}
				else if (val<-900 && val2>900)
				{
					dirty->put(i+1,j,k,-937);
				}
				else if (val<-400 && val2>900)
				{
					dirty->put(i,j,k,-937);
					dirty->put(i+1,j,k,-937);
				}
				//else if (val2>800&& val3>800&& val4>800&&val5>800)
				//{
				//	dirty->put(i,j,k,-937);
				//}

			}//..end for k
		}//..end for j
	}//..end for i
	for (int i=dirty->getXsize()-1;i>2; i--)
	{
		for (int j=2; j <dirty->getYsize()-1; j++)
		{
			for (int k=10 ; k<dirty->getZsize()-62; k++)
			{
				PIXTYPE val=dirty->get(i,j,k);
				PIXTYPE val2= dirty->get(i+1,j,k);
				PIXTYPE val3= dirty->get(i-1,j,k);
				if (val3>800&& val<-300)
				{
					dirty->put(i-1,j,k,-937);
				}

			}//..end for k
		}//..end for j
	}//..end for i
}

void test()
{
	int l = 0, m = 0, n = 0;
	RawImage test;
	short*orgiondata=test.readStream("F:\\data\\dirty\\bone\\WIS_4859_2_P.raw",&l,&m,&n);
	PIXTYPE *indata=new PIXTYPE[l*m*n];
	for (int i=0;i< l*m*n ;i++)
	{
		indata[i]=(PIXTYPE)orgiondata[i];
	}
	delete [] orgiondata;
	Raw * dirty =new Raw(l,m,n,indata);
	processDirty(dirty);
	test.writenormal(*dirty,"F:\\data\\dirty\\test4859_2.raw");
}