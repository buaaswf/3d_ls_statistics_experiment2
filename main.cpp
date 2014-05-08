#include "vol_math_LevelSet.h"
#include "statistics.h"
#include "test.h"
#include "ProcessDirty.h"
//#include <iostream>
//#include <crtdbg.h> 
//#include "CImg.h" 
//#include "ThreeDim_LevelSet.h"
//#include "Filter.h"
//#include "WipeNioisePde.h"
//#include "string.h"
////================
//#include <fstream>
//#include <vector>
//#include <string>
//#ifdef _WIN32
//#include <Windows.h>
//#include <strsafe.h>
//#else
//#include <dirent.h>
//#endif
//#define output "D:\\sdfdata\\" 
//#define input1  "D:\\sdfdata\\pvaluethickness\\" //"K:\\sdf\\volume\\clean\\clean\\ep\\""E:\\volume\\cleantestdata\\test2\\"
//#define input2  "K:\\sdf\\volume\\clean\\clean\\ep\\test3\\" 
//using namespace cimg_library;
//using namespace std;


void testcolon(int argc,string dir)
{
	
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=10;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	for (int i = 0; i < l*m*n; i++)
	{
		inputo[i]=(float) indata[i];		
	}

	Raw *input=new Raw(l,m,n,inputo);
	Filter *f=new Filter();
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//20140405 delete because of the existance of 
	ls->initialg(*input);
	for (int i=0; i<input->getXsize(); i++)
	{
		for (int j=0; j<input->getYsize(); j++)
		{
			for (int k=0; k<input->getZsize(); k++)
			{
				if (input->get(i,j,k)>=1)
				{
					initial->put(i,j,k,-2);
				}
				else initial->put(i,j,k,2);

			}
		}

	}
	*initial=ls->minimal_surface(*initial,*input,5.0,0.1,-3,1.5,1,iter_outer,pt);
	char *outname1="inner5-8_2.raw";
	char outdir[200]=output;

	strcat(outdir,dirbody);
	strcat(outdir,outname1);
	//test.readImage2(initial->getdata(),outdir,l*m*n);
	test.writeImageName(*initial,outdir);
	//Raw temp(*initial);
	ls->outerwall(*initial,*input,5.0,0.1,-3,1.5,1,10,pt);
	//*initial -=temp;
	char *outname2="outer5-8_2_20140405.raw";
	char outdir2[200]=output;
	strcat(outdir2,dirbody);
	strcat(outdir2,outname2);
	test.writeImageName(*initial,outdir2);
	evaluate(dir,l,m,n);
}
void testcolontest()
{
	char *pt="single_well";
	int l = 512,m = 512,n = 144, l1=0,l2=0,iter_outer = 10;
	RawImage test;
	char dirbody[100];

	short * indata;///=test.readMRI("K:\\sdf\\MRI\\SE1512_512_144.raw",&l,&m,&n);//F:\\PA1\\ST1\\SE1\\  //K:\\sdf\\MRI\\
	
	indata =new short[l*m*n];
	for (int k =0; k < n; k++ )
	{
		for ( int j = 0; j < m; j++)
		{
			for ( int i = 0; i < l; i++)
			{
				//if ( i < 300 && i > 100 && j < 300 && j > 200 && k > 30 && k < 60 )
			
				if ( (i >= 190 && i <= 230 && j >= 200 && j <= 270 && k > 35 && k < 55))
				{
					indata[i + j*l + k*m*l] = 1;
				} 
				else
				{
					indata[i + j*l + k*m*l] = 0;
				}
			}
		}
		
	}
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	short min = 1000,max = -100;
	for (int i = 0; i < l*m*n; i++)
	{
		float * p= (float *)(indata+i);
		//unsigned char * bp= (unsigned char *)p;
		//std:swap(bp[0],bp[3]);
		//std::swap(bp[1],bp[2]);
		min < indata[i] ? min=min:min=indata[i];
		max > indata[i] ? max=max:max=indata[i];
		
	/*	if ( indata[i] >= 864 && indata[i] <= 1063 )
		{
			inputo[i] = 100;
		} 
		else
		{
			inputo[i] = (short )0;
		}
		*/
		inputo[i]=(short) indata[i];		
	}

	cout <<min << max <<endl;

	Raw *input=new Raw(l,m,n,inputo);

	//Filter *f=new Filter();
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//ls->initialg(*input);
	for (int i=0; i<input->getXsize(); i++)
	{
		for (int j=0; j<input->getYsize(); j++)
		{
			for (int k=0; k<input->getZsize(); k++)
			{
				//if (input->get(i,j,k) >= 1)
				//{
				//	initial->put(i,j,k,-2);
				//}
				//else 
					//if ((i >= 172 && i <= 352 && j >= 164 && j <= 376 && z>19 && z <))
				if ((i >= 196 && i <= 220 && j >= 202 && j <= 267 && k > 40 && k < 50))
				{
					initial->put(i, j, k, -2);
				} 
				else
				{
					initial->put(i, j, k, 2);
				}


			}
		}

	}
	*initial=ls->minimal_surface(*initial,*input,5.0,0.1,-3,1.5,1,iter_outer,pt);
	//if you available this, don,t
	//forget to change the next line to initial
	test.writeMRI(*initial,"K:\\sdf\\MRI\\SE1512_512_144.testdata.raw");//F:\\PA1\\ST1\\SE1

}
void testsesmic()
{
	char *pt="single_well";
	int l = 201,m = 201,n = 851, l1=0,l2=0,iter_outer = 10;
	RawImage test;
	char dirbody[100];
	unsigned char * indata=new unsigned char[l*m*n];
	test.readImage(indata,"K:\\sdf\\geo\\Probe_fault_Amp.probe .raw",l*m*n);//F:\\PA1\\ST1\\SE1\\  //K:\\sdf\\MRI\\
	//
	//indata =new short[l*m*n];
	//for (int k =0; k < n; k++ )
	//{
	//	for ( int j = 0; j < m; j++)
	//	{
	//		for ( int i = 0; i < l; i++)
	//		{
	//			//if ( i < 300 && i > 100 && j < 300 && j > 200 && k > 30 && k < 60 )
	//		
	//			if ( (i >= 190 && i <= 230 && j >= 200 && j <= 270 && k > 35 && k < 55))
	//			{
	//				indata[i + j*l + k*m*l] = 1;
	//			} 
	//			else
	//			{
	//				indata[i + j*l + k*m*l] = 0;
	//			}
	//		}
	//	}
	//	
	//}
	Raw *initial=new Raw(l,m,n);
	float *inputo=new float[l*m*n];
	short min = 1000,max = -100;
	for (int i = 0; i < l*m*n; i++)
	{
		//change the big --little
		//float * p= (float *)(indata+i);
		//unsigned char * bp= (unsigned char *)p;
		//std:swap(bp[0],bp[3]);
		//std::swap(bp[1],bp[2]);
		min < indata[i] ? min=min:min=indata[i];
		max > indata[i] ? max=max:max=indata[i];
		//cal the max and min data
	/*	if ( indata[i] >= 864 && indata[i] <= 1063 )
		{
			inputo[i] = 100;
		} 
		else
		{
			inputo[i] = (short )0;
		}
		*/
		inputo[i]=(float) indata[i];		
	}

	cout <<min << max <<endl;

	Raw *input=new Raw(l,m,n,inputo);

	//Filter *f=new Filter();
	//input=f->guass3DFilter(input,3);
	RawImage *write=new RawImage();
	ThreeDim_LevelSet *ls=new ThreeDim_LevelSet();
	//ls->initialg(*input);
	for (int i=0; i<input->getXsize(); i++)
	{
		for (int j=0; j<input->getYsize(); j++)
		{
			for (int k=0; k<input->getZsize(); k++)
			{
				//if (input->get(i,j,k) >= 1)
				//{
				//	initial->put(i,j,k,-2);
				//}
				//else 
					//if ((i >= 172 && i <= 352 && j >= 164 && j <= 376 && z>19 && z <))
				if ((i >= 196 && i <= 220 && j >= 202 && j <= 267 && k > 40 && k < 50))
				{
					initial->put(i, j, k, -2);
				} 
				else
				{
					initial->put(i, j, k, 2);
				}


			}
		}

	}
	*initial=ls->minimal_surface(*initial,*input,5.0,0.1,-3,1.5,1,iter_outer,pt);
	//if you available this, don,t
	//forget to change the next line to initial
	test.writeMRI(*initial,"K:\\sdf\\geo\\data.raw");//F:\\PA1\\ST1\\SE1

}
void testhistgram()
{
	HUandThickness();
}
void rate(int argc,string dir)
{
	char *pt="single_well";
	int l=0,m=0,n=0,l1=0,l2=0,iter_outer=10;
	RawImage test;
	char dirhead[200]=input2;  //K:\\sdf\\volume\\clean\\clean\\ep\\

	char dirbody[100];
	strcpy(dirbody,dir.c_str());
	cout <<"dirbody"<<dirbody<<endl;
	strcat(dirhead,dirbody);
	cout << "dirhead" << dirhead<< endl;
	short * indata=test.readStream(dirhead,&l,&m,&n);
	delete [] indata;
	
}
int main(int argc,char **argv)
{
	string dir2(input2);
	vector<string> files2;
	GetFileNameFromDir(dir2,files2);
	vector<string>::iterator iterFile2;
	for ( iterFile2 = files2.begin(); iterFile2 != files2.end(); iterFile2++ )
	{

		
		iterFile2->assign(iterFile2->substr(dir2.size()+1));
		cout<<*iterFile2 <<endl;
		//ddcircle(*iterFile);
		//testcolon(argc,*iterFile2);
		//float2uchar(512,512,700,*iterFile2);
		//testsesmic();
		//thincknessstdv2(*iterFile2);
		//roc(*iterFile2);
		//rate(argc,*iterFile2);
	}
	testhistgram();
	//test();//delete dirty success
	//cout<<endl;
	//roc3();
	//threshold();
	//rocway2();
	//testcolontest();
	//testsesmic();
	system("pause");
	return 0;
}
