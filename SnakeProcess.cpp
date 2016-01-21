   // SnakeProcess.cpp: implementation of the SnakeProcess class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "mammocad.h"
#include "SnakeProcess.h"
#include "ROI.h"
#include "imageprocess.h"
#include <fstream.h>
#include <list>
#include "Matrix.h"
#include <afxtempl.h>
#include <VECTOR>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
extern void gradient(USHORT*f ,	// image buffer of growth region defined by active contour
			    int rows, int cols);
SnakeProcess::SnakeProcess()
{

}

SnakeProcess::~SnakeProcess()
{

}
 // Snake.cpp : Defines the class behaviors for the application.
//



/////////////////////////////////////////////////////////////////////////////
// The one and only CSnakeApp object

LONG glHeight,glWidth;
bool glFlag;
struct DIRECTION
{
	bool direction[8];
	DIRECTION()
	{
		for(int i=0;i<8;i++)
			direction[i]=FALSE;
	}
};

struct Ppoint 
{
	int x,y;
	Ppoint(int a=-1,int b=-1)
	{
		x=a;
		y=b;
	}
	
	bool operator==(Ppoint tempPoint)
	{
		return (x==tempPoint.x && y == tempPoint.y);
	}
	Ppoint& operator=(Ppoint& tempPoint)
	{
		x=tempPoint.x;
		y=tempPoint.y;
		return *this;
	}
};
/* NOTE:
 * 
 * You must have "mex" command working in your matlab.  In matlab, 
 * type "mex gvfc.c" to compile the code. See usage for details of 
 * calling this function.
 *
 * Replace GVF(...) with GVFC in the snakedeform.m. The speed is 
 * significantly faster than the Matlab version.
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"
#define PG_MAX(a,b)  (a>b? a: b)
#define PG_MIN(a,b)  (a<b? a: b)
#define MAX_12_BIT 4095

typedef unsigned char  PGbyte;
typedef          char  PGchar;
typedef unsigned short PGushort;
typedef          short PGshort; 
typedef unsigned int   PGuint;
typedef          int   PGint;
typedef float          PGfloat;
typedef double         PGdouble;
typedef void           PGvoid;
typedef unsigned short USHORT;


/* Prints error message: modified to work in mex environment */
void MinQuare(USHORT *Para_a, USHORT *Para_b, USHORT *Para_c,USHORT *Para_d,USHORT *pDataSrc)
{
	int i,j,k,l;
	int nRows,nCols;
	nRows = 125;
	nCols = 125;
	double* pData_a = new double[nRows*nCols];
	double* pData_b = new double[nRows*nCols];
	double* pData_c = new double[nRows*nCols];
	double* pData_d = new double[nRows*nCols];

	int sizeWnd=5;
	for(k=sizeWnd/2;k<nRows-sizeWnd/2;k++)
	{
		for(l=sizeWnd/2;l<nCols-sizeWnd/2;l++)
		{
			double a,b,c;	
			double sumII,sumJJ,sumIJ,sumI,sumJ,sum1;
			sumII=sumJJ=sumIJ=sumI=sumJ=sum1=0;
			double sumIPij,sumJPij,sumPij;
			sumIPij=sumJPij=sumPij=0;
			for(i=k-sizeWnd/2;i<k+sizeWnd/2+1;i++)
			{
				for(j=l-sizeWnd/2;j<l+sizeWnd/2+1;j++)
				{
					sumII += (i-k)*(i-k);
					sumJJ += (j-l)*(j-l);
					sumIJ += (i-k)*(j-l);
					sumI  += (i-k);
					sumJ  += (j-l);
					sum1  += 1;

					sumIPij += (double)(i-k) * (*(pDataSrc+125*i+j));
					sumJPij += (double)(j-l) *(*(pDataSrc+125*i+j));
					sumPij  += (*(pDataSrc+125*i+j));
				}
			}
			double value[9]={
				sumII,sumIJ,sumI,
				sumIJ,sumJJ,sumJ,
				sumI ,sumJ ,sum1   
			};
			CMatrix matrix(3,3,value);
			BOOL flag = matrix.InvertSsgj(); 
			if(flag)
			{
				double valueTmp[3]={
					sumIPij,sumJPij,sumPij
				};
				CMatrix matrixTmp(3,1,valueTmp);
				CMatrix para(3,1);
				para = matrix*matrixTmp;
				a=para.GetElement(0,0);
				b=para.GetElement(1,0);
				c=para.GetElement(2,0);
			}
			pData_a[k*nCols+l]=a;
			pData_b[k*nCols+l]=b;
			pData_c[k*nCols+l]=c;
			double sin,cos,r;
			r=sqrt((k-nRows/2.0)*(k-nRows/2.0)+(l-nCols/2.0)*(l-nCols/2.0));
			sin = (nRows/2.0-k)/r;
			cos = (l-nCols/2.0)/r;
			pData_d[k*nCols+l]=a*sin-b*cos;
		}//End for(l=2;l<nCols-2;l++)
	}//End for(k=2;k<nRows-2;k++)
	//进行标定
	float min_a,min_b,min_c,min_d,max_a,max_b,max_c,max_d;
	min_a=max_a=pData_a[sizeWnd/2*nCols+sizeWnd/2];
	min_b=max_b=pData_b[sizeWnd/2*nCols+sizeWnd/2];
	min_c=max_c=pData_c[sizeWnd/2*nCols+sizeWnd/2];
	min_d=max_d=pData_d[sizeWnd/2*nCols+sizeWnd/2];
	
	for(k=sizeWnd/2;k<nRows-sizeWnd/2;k++)
	{
		for(l=sizeWnd/2;l<nCols-sizeWnd/2;l++)
		{
			//Para_a
			if(pData_a[k*nCols+l]<min_a)
			{
				min_a=pData_a[k*nCols+l];
			}
			else if(pData_a[k*nCols+l]>max_a)
			{
				max_a=pData_a[k*nCols+l];
			}
            //Para_b
			if(pData_b[k*nCols+l]<min_b)
			{
				min_b=pData_b[k*nCols+l];
			}
			else if(pData_b[k*nCols+l]>max_b)
			{
				max_b=pData_b[k*nCols+l];
			}
            //Para_c
			if(pData_c[k*nCols+l]<min_c)
			{
				min_c=pData_c[k*nCols+l];
			}
			else if(pData_c[k*nCols+l]>max_c)
			{
				max_c=pData_c[k*nCols+l];
			}
			//Para_d
			if(pData_d[k*nCols+l]<min_d)
			{
				min_d=pData_d[k*nCols+l];
			}
			else if(pData_d[k*nCols+l]>max_d)
			{
				max_d=pData_d[k*nCols+l];
			}
		}//End for(l=2;l<nCols-2;l++)
	}//End for(k=2;k<nRows-2;k++)
	for(k=sizeWnd/2;k<nRows-sizeWnd/2;k++)
	{
		for(l=sizeWnd/2;l<nCols-sizeWnd/2;l++)
		{
			//Para_a  
            *(Para_a+125*k+l) =  4095*(pData_a[k*nCols+l]-min_a)/(max_a-min_a);
            //Para_b
			*(Para_b+125*k+l) =  4095*(pData_b[k*nCols+l]-min_b)/(max_b-min_b);
            //Para_c
			*(Para_c+125*k+l) =  4095*(pData_c[k*nCols+l]-min_c)/(max_c-min_c);
            //Para_d	

			*(Para_d+125*k+l) =  4095*(pData_d[k*nCols+l]-min_d)/(max_d-min_d);
		}//End for(l=2;l<nCols-2;l++)
	}//End for(k=2;k<nRows-2;k++)
}
PGvoid pgError(char error_text[])//错误消息处理
{
   char buf[255];
   
   sprintf(buf, "GVFC.MEX: %s", error_text);
   //mexErrMsgTxt(buf);
   
   return;
}


/* Allocates a matrix of doubles with range [nrl..nrh][ncl..nch] */
PGdouble **pgDmatrix(int nrl, int nrh, int ncl, int nch)//分配二维数组，数组下标可以任意开始A[3..10][5..15]
{
   int j;
   long bufsize,bufptr;
   PGdouble **m;
   nrh++;
   nch++;
   bufsize = (nrh-nrl+1)*sizeof(PGdouble*)
      + (nrh-nrl+1)*(nch-ncl+1)*sizeof(PGdouble);//节省空间，sizeof(PGdouble*)=4
   
   m=(PGdouble **) malloc(bufsize);
   if (!m) pgError("allocation failure 1 in pgDmatrix()");
   m -= nrl;//如果此处不这样写的话，下面m[j]就不是那个位置了!
   
   bufptr = ((long) (m)) + (nrh-nrl+1)*sizeof(PGdouble*);//相当于到数据区的偏移量
   for(j=nrl;j<=nrh;j++) {
      m[j] = ((PGdouble *) (bufptr+(j-nrl)*(nch-ncl+1)*sizeof(PGdouble)));//此时相当于m+j
      m[j] -= ncl;//此时将其地址向前移动,目的是调用时m[nrl][ncl]相当于m[0][0]
   }
   
   return m;
}


/* Frees a matrix allocated by dmatrix */
PGvoid pgFreeDmatrix(PGdouble **m, int nrl, int nrh, int ncl, int nch)//二维数组指针释放
{
   free((char*) (m+nrl));
   return;
}

double *dvector(long nl, long nh)
{
	double *v;
	v=(double *)malloc((size_t)((nl*nh)*sizeof(double)));
	if (!v)
	{
		AfxMessageBox("内存分配失败！");
		//exit(1);
	}
	return (v);
}


int *ivector(long nl, long nh)//申请int类型向量空间
{
	int *v;

	v = (int *)malloc((size_t) ((nh-nl+1+1)*sizeof(int)));
	if( !v )
		exit(1);

	return v-nl+1;//?
}

USHORT *svector(long nl, long nh)//申请USHORT类型向量空间//?
{
	USHORT *v;

	v = (USHORT *)malloc((size_t) ((nl*nh)*sizeof(USHORT)));
	if( !v )
		exit(1);

	return v;

}

void free_dvector(double *v, long nl, long nh)//释放double类型向量空间
{
	free((char*) (v));
}

void free_ivector(int *v, long nl, long nh)//释放int类型向量空间
{
	free((char*) (v+nl-1));
}
void free_svector(USHORT *v)//释放double类型向量空间//?
{
	free(v);
}

void copy_svector(USHORT* img_buf, USHORT* EdgeMap_buf, int rows,int cols)
{
	memcpy(EdgeMap_buf, img_buf, rows*cols*sizeof(USHORT));

}
//输入：图像的高和宽：int YN, int XN  受限的梯度图像：double *f　固定参数：mu＝0.1 ITER=80 ，mu越大GVF场越模糊；
//输出：U场和V场：double *ou, double *ov ；
//功能：创建GVF场 ；
void GVFC(int YN, int XN, double *f, double *ou, double *ov, 
	  double mu, int ITER)
{
   double mag2, temp, tempx, tempy, fmax, fmin;
   double *d_base, *d_line, *d_buf1, *d_buf2, *d_buf3, *d_buf4;
   int count, x, y, XN_1, XN_2, YN_1, YN_2;
   
   PGdouble **fx, **fy, **u, **v, **Lu, **Lv, **g, **c1, **c2, **b;
   
   /* define constants and create row-major double arrays */
   XN_1 = XN - 1;
   XN_2 = XN - 2;
   YN_1 = YN - 1;
   YN_2 = YN - 2;
   fx = pgDmatrix(0,YN_1,0,XN_1);
   fy = pgDmatrix(0,YN_1,0,XN_1);
   u  = pgDmatrix(0,YN_1,0,XN_1);
   v  = pgDmatrix(0,YN_1,0,XN_1);
   Lu = pgDmatrix(0,YN_1,0,XN_1);
   Lv = pgDmatrix(0,YN_1,0,XN_1);
   g  = pgDmatrix(0,YN_1,0,XN_1);
   c1 = pgDmatrix(0,YN_1,0,XN_1);
   c2 = pgDmatrix(0,YN_1,0,XN_1);
   b  = pgDmatrix(0,YN_1,0,XN_1);
   
   /************** I: Normalize the edge map to [0,1] **************/
   fmax = -4095.0;
   fmin = 4095;
   count = XN * YN;
   x = 0;
   while( count -- )
   {
      fmax = PG_MAX(fmax,f[x]);
      fmin = PG_MIN(fmin,f[x]);
	  x ++;
   }

   //if (fmax==fmin) 
     // mexErrMsgTxt("Edge map is a constant image.");      

   count = XN * YN;
   x = 0;
   while( count -- )
   {
      f[x] = (f[x] - fmin) / (fmax - fmin);
	  x ++;
   }//Normalize the intensity values of the edge map so that their values fall between 0 and 1

   /**************** II: Compute edge map gradient *****************/
   /* I.1: Neumann boundary condition: 
    *      zero normal derivative at boundary 
    //处理边界梯度
   /* Deal with corners */
   fx[0][0] = fy[0][0] = fx[0][XN_1] = fy[0][XN_1] = 0;
   fx[YN_1][XN_1] = fy[YN_1][XN_1] = fx[YN_1][0] = fy[YN_1][0] = 0;

   /* Deal with left and right column */
   for (y=1; y<YN_1; y++) 
   {
      fx[y][0] = fx[y][XN_1] = 0;
      fy[y][0] = 0.5 * (f[y+1] - f[y-1]);//中心差分求梯度
      fy[y][XN_1] = 0.5 * (f[y+1 + XN_1*YN] - f[y-1 + XN_1*YN]);//中心差分求梯度
   }

   /* Deal with top and bottom row */
   for (x=1; x<XN_1; x++) 
   {
      fy[0][x] = fy[YN_1][x] = 0;
      fx[0][x] = 0.5 * (f[(x+1)*YN] - f[(x-1)*YN]);//中心差分求梯度
      fx[YN_1][x] = 0.5 * (f[YN_1 + (x+1)*YN] - f[YN_1 + (x-1)*YN]);//中心差分求梯度
   }
   
   /* I.2: Compute interior derivative using central difference *///处理图像内部梯度
   d_base = f+XN+1;
   for (y=1; y<=YN_2; y++)
   {
	   d_line = d_base;
       for (x=1; x<=XN_2; x++)
	   {
		   /* NOTE: f is stored in column major */
		   d_buf1 = d_line - 1;
		   d_buf2 = d_line + 1;
		   d_buf3 = d_line - XN;
		   d_buf4 = d_line + XN;
		   fx[y][x] = 0.5 * (*d_buf2 - *d_buf1); //中心差分求梯度
	       fy[y][x] = 0.5 * (*d_buf4 - *d_buf3); //中心差分求梯度
		   d_line ++;
       }
	   d_base += XN;
   }


    
   /******* III: Compute parameters and initializing arrays **********/
   temp = -1.0/(mu*mu);//?无用
   for (y=0; y<=YN_1; y++)
   {
      for (x=0; x<=XN_1; x++) 
	  {
		  tempx = fx[y][x];//x方向梯度
	      tempy = fy[y][x];//y方向梯度
	      /* initial GVF vector */
	      u[y][x] = tempx;
	      v[y][x] = tempy;
	      /* gradient magnitude square */
	      mag2 = tempx*tempx + tempy*tempy; //总梯度的模的平方
	 
	      g[y][x] = mu;
	      b[y][x] = mag2;

	      c1[y][x] = b[y][x] * tempx;//总梯度的模的平方*x方向梯度
	      c2[y][x] = b[y][x] * tempy;//总梯度的模的平方*y方向梯度
      }
   }
   
   /* free memory of fx and fy */
   pgFreeDmatrix(fx,0,YN_1,0,XN_1);
   pgFreeDmatrix(fy,0,YN_1,0,XN_1);

   /************* Solve GVF = (u,v) iteratively ***************/
   for (count=1; count<=ITER; count++) //80次
   {
      /* IV: Compute Laplace operator using Neuman condition */
      /* IV.1: Deal with corners */
      Lu[0][0] = (u[0][1] + u[1][0])*0.5 - u[0][0]; 
      Lv[0][0] = (v[0][1] + v[1][0])*0.5 - v[0][0];
      Lu[0][XN_1] = (u[0][XN_2] + u[1][XN_1])*0.5 - u[0][XN_1];
      Lv[0][XN_1] = (v[0][XN_2] + v[1][XN_1])*0.5 - v[0][XN_1];
      Lu[YN_1][0] = (u[YN_1][1] + u[YN_2][0])*0.5 - u[YN_1][0];
      Lv[YN_1][0] = (v[YN_1][1] + v[YN_2][0])*0.5 - v[YN_1][0];
      Lu[YN_1][XN_1] = (u[YN_1][XN_2] + u[YN_2][XN_1])*0.5 - u[YN_1][XN_1];
      Lv[YN_1][XN_1] = (v[YN_1][XN_2] + v[YN_2][XN_1])*0.5 - v[YN_1][XN_1];
      
      /* IV.2: Deal with left and right columns */
      for (y=1; y<=YN_2; y++) 
	  {
		  Lu[y][0] = (2*u[y][1] + u[y-1][0] + u[y+1][0])*0.25 - u[y][0];
	      Lv[y][0] = (2*v[y][1] + v[y-1][0] + v[y+1][0])*0.25 - v[y][0];
		  Lu[y][XN_1] = (2*u[y][XN_2] + u[y-1][XN_1] 
				+ u[y+1][XN_1])*0.25 - u[y][XN_1];
		  Lv[y][XN_1] = (2*v[y][XN_2] + v[y-1][XN_1] 
				+ v[y+1][XN_1])*0.25 - v[y][XN_1];
      }
      
      /* IV.3: Deal with top and bottom rows */
      for (x=1; x<=XN_2; x++) 
	  {
		  Lu[0][x] = (2*u[1][x] + u[0][x-1] + u[0][x+1])*0.25 - u[0][x];
		  Lv[0][x] = (2*v[1][x] + v[0][x-1] + v[0][x+1])*0.25 - v[0][x];
		  Lu[YN_1][x] = (2*u[YN_2][x] + u[YN_1][x-1] 
				+ u[YN_1][x+1])*0.25 - u[YN_1][x];
		  Lv[YN_1][x] = (2*v[YN_2][x] + v[YN_1][x-1] 
				+ v[YN_1][x+1])*0.25 - v[YN_1][x];
      }
      
      /* IV.4: Compute interior */
      for (y=1; y<=YN_2; y++)
	  {
		  for (x=1; x<=XN_2; x++) 
		  {
			  Lu[y][x] = (u[y][x-1] + u[y][x+1] 
					+ u[y-1][x] + u[y+1][x])*0.25 - u[y][x];
			  Lv[y][x] = (v[y][x-1] + v[y][x+1]
					+ v[y-1][x] + v[y+1][x])*0.25 - v[y][x];
		  }
	  }
      
      /******** V: Update GVF ************/
      for (y=0; y<=YN_1; y++)
	  {
		  for (x=0; x<=XN_1; x++) 
		  {
			  u[y][x] = (1- b[y][x])*u[y][x] + g[y][x]*Lu[y][x]*4 + c1[y][x];
			  v[y][x] = (1- b[y][x])*v[y][x] + g[y][x]*Lv[y][x]*4 + c2[y][x];
		  }//t=1时间步长为1
	  }
      
   } 
   /* copy u,v to the output in column major order */
   d_buf1 = ou;
   d_buf2 = ov;
   for (y=0; y<=YN_1; y++)
   {
      for (x=0; x<=XN_1; x++)
	  {
		*d_buf1 ++ = u[y][x];
		*d_buf2 ++ = v[y][x];
      }
   }

   /* free all the array memory */
   pgFreeDmatrix(u,0,YN_1,0,XN_1);
   pgFreeDmatrix(v,0,YN_1,0,XN_1);
   pgFreeDmatrix(Lu,0,YN_1,0,XN_1);
   pgFreeDmatrix(Lv,0,YN_1,0,XN_1);
   pgFreeDmatrix(g,0,YN_1,0,XN_1);
   pgFreeDmatrix(c1,0,YN_1,0,XN_1);
   pgFreeDmatrix(c2,0,YN_1,0,XN_1);
   pgFreeDmatrix(b,0,YN_1,0,XN_1);

   return;
}


// functions copied from Numerical Recipes in C
//输出错误信息到error_record.txt
#define TINY  1.0e-20

void nrerror(char* chmsg)
{
	FILE *fpt;

	if((fpt=fopen("error_record.txt", "w"))==NULL)
		exit( 1 );

	fprintf(fpt, "Error Message: %s...\n", chmsg);
	exit( 1 );
	fclose(fpt);

	return;
}
//inverse the matrix
void ludcmp(double **a, int n, int *indx, double *d)//和ludksb共同实现逆矩阵
{
	int    i, imax, j , k;
	double big, dum, sum, temp;
	double *vv;

	vv = dvector(1, n+1);

	*d = 1.0;

	for(i=1; i<=n; i++)
	{
		big = 0.0;
		for(j=1; j<=n; j++)
			if((temp=fabs(a[i][j])) > big)
				big = temp;

		if(big == 0.0)
			nrerror(" Singular matrix in routine ludcmp! ");

		vv[i] = 1.0 / big;
	}

	for(j=1; j<=n; j++)
	{
		for(i=1; i<j; i++)
		{
			sum = a[i][j];
			for(k=1; k<i; k++)
				sum -= (a[i][k] * a[k][j]);
			a[i][j] = sum;
		}

		big = 0.0;
		for(i=j; i<=n; i++)
		{
			sum = a[i][j];
			for(k=1; k<j; k++)
				sum -= (a[i][k] * a[k][j]);
			a[i][j] = sum;

			if((dum = vv[i] * fabs(sum)) >= big)
			{
				big = dum;
				imax = i;
			}//求最大值位置
		}

		if(j != imax)
		{
			for(k=1; k<=n; k++)
			{
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}

		indx[j] = imax;
		if(a[j][j] == 0.0) 
			a[j][j] = TINY;

		if(j != n)
		{
			dum = 1.0 / a[j][j];
			for(i=j+1; i<=n; i++)
				a[i][j] *= dum;
		}
	}

	free_dvector(vv,1,n+1);

	return;
}

//和ludcmp共同实现逆矩阵
void lubksb(double **a, int n, int *indx, double b[])
{
	int    i, ii, ip, j;
	double sum;

	ii = 0;

	for(i=1; i<=n; i++)
	{
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if( ii )
		{
			for(j=ii; j<=i-1; j++)
				sum -= (a[i][j] * b[j]);
		}
		else if( sum )
		{
			ii = i;
		}

		b[i] = sum;
	}

	for(i=n; i>=1; i--)
	{
		sum = b[i];
		for(j=i+1; j<=n; j++)
			sum -= (a[i][j] * b[j]);

		b[i] = sum / a[i][i];
	}

	return;
}

//输入：轮廓点x、y坐标数组， 轮廓点数目 array_length；
//输出：轮廓点梯度图 diff；
//功能：求轮廓点序列的差分图；
void ComputeDifference(int array_length, 
					   double* diff, 
					   double* x, double* y
					  )//计算每个点的梯度－梯度图
{
	int    i, N;
	double dx, dy;
	double x_shift[24000], y_shift[24000];

	// shift array position by one bit.

	N = array_length - 1;

	for(i=0; i<N; i++)
	{
		x_shift[i] = x[i+1];
		y_shift[i] = y[i+1];
	}
	x_shift[N] = x[0];
	y_shift[N] = y[0];

	// compute difference

	for(i=0; i<=N; i++)
	{
		dx = x_shift[i] - x[i];//利用前向差分
		dy = y_shift[i] - y[i];//利用前向差分
		diff[i] = fabs(dx) + fabs(dy);  //差分的绝对值之和；
	}

	return;
}
//将double类型数组拷贝到double类型数组中
void CopyDDArray(int N, double* array1, double* array2)
{
	int i;

	for(i=0; i<N; i++)
		array2[i] = array1[i];

	return;
}
//将double类型数组拷贝到int类型数组中
void CopyDIArray(int N, double* array1, int* array2)
{
	int i;

	for(i=0; i<N; i++)
		array2[i] = (int)(array1[i]);

	return;
}
//输入：轮廓点数目：N　　轮廓点梯度图：diff　　标志：index　阈值：threshold　；
//输出：IDX ；
//功能：根据index和threshold将diff数组(梯度图)二值化后存于IDX中
void DefineIDX(int N, int* IDX, double* diff, int index, double threshold)
{
	int i;

	if(index > 0)    // 如果index > 0 ，那么IDX中存的是插值索引 ，此时：threshold=max_dist
	{
		for(i=0; i<N; i++)
		{
			if(diff[i] > threshold)
				IDX[i] = 1;
			else
				IDX[i] = 0;
		}
	}
	else             //   如果index <= 0 ，那么IDX中存的是删除多余点索引 ，此时：threshold=min_dist
	{
		for(i=0; i<N; i++)
		{
			if(diff[i] < threshold)
				IDX[i] = 0;
			else
				IDX[i] = 1;
		}
	}

	return;
}

//寻找diff_array中最大值，并返回；
double FindMaxValue(int N, double* diff_array)
{
	int    i;
	double max_value;

	max_value = 0.0;
	for(i=0; i<N; i++)
	{
		if(diff_array[i] > max_value)
			max_value = diff_array[i];
	}

	return( max_value );
}
//输入：插值标识数组：IDX   轮廓点x或者y方向数组：array_i　　轮廓点的数目：array_length　；
//输出：插值后的轮廓点x或y方向数组：array_d　；
//功能：线性插值，用相邻元素取均值进行插值；
int interp1(int* IDX, double* array_i, int array_length,
			 double* array_d
			)
{
	int  i, k;
	int  new_length;
	k = array_length - 1;
	new_length = 0;
	for(i=0; i<k; i++)
	{
		array_d[new_length ++] = array_i[i];
   		if(IDX[i] == 1)  //如果标识数组的元素值IDX[i] == 1 ，那么就要进行插值；
		{
			array_d[new_length++] = (array_i[i] + array_i[i+1]) *0.5;

		}
	}
	array_d[new_length++] = array_i[k];
	if(IDX[k] == 1)  //处理最后一个点；
	{
		array_d[new_length++] = (array_i[k] + array_i[0]) *0.5;
	}
	return( new_length );
}
//功能：实现双线性插值；
void interp2(double* EFvec_buf, double* array_x, double* array_y, double* vfx, 
			 int array_length, int rows, int cols) 
{
	int      i, j, k, m;
	double   x1, x2;
	double   y1, y2, y3, y4, t, u;
	double   a1, a2, a3, a4;
	double   *d_buf;
	PGdouble **fx;

	fx = pgDmatrix(1, rows, 1, cols);

	d_buf = EFvec_buf;
	for(j=1; j<=rows; j++)
	{
		for(i=1; i<=cols; i++)
		{
			fx[j][i] = *d_buf ++;
		}
	}

	m = array_length;
	if(m==1)
	{
		pgFreeDmatrix(fx, 1, rows, 1, cols);
		return ;
	}

	for(i=0; i<m; i++)
	{
		x1 = array_y[i];
		x2 = array_x[i];
		j = (int)(x1);
		k = (int)(x2);
		if(j<1||j>rows-1||k<1||k>cols-1)
		{
			pgFreeDmatrix(fx, 1, rows, 1, cols);
			return ;
		}

		y1 = fx[j][k];
		y2 = fx[j+1][k];
		y3 = fx[j+1][k+1];
		y4 = fx[j][k+1];

		t = x1 - (double)(j);
		u = x2 - (double)(k);
 //双线性插值
		a1 = (1.0 - t) * (1.0 - u) * y1;
		a2 = t * (1.0 - u) * y2;
		a3 = t * u * y3;
		a4 = (1.0 - t) * u * y4;
		vfx[i] = a1 + a2 + a3 + a4;
	}

	pgFreeDmatrix(fx, 1, rows, 1, cols);

	return;
}

//输入：轮廓点x、y坐标数组 contour_x、contour_y , 轮廓点的数目 contour_length ，
//      最大距离和最小距离max_dist 和min_dist ；
//输出：插值之后的轮廓x、y坐标数组 contour_x、contour_y
//返回：轮廓点的数目；
//功能：snake演变之后插值；
int SnakeInterp(double* contour_x, double* contour_y, int contour_length,
				 double max_dist, double min_dist
				)
{
	int    i, k, m, N;
	int    Interp_length;
	double max_d;
	double *diff;
	int    *IDX;
	double *xi, *yi;
	double *xo, *yo;

	// convert to double format
	N = contour_length;
	m = 30 * N;
    //申请空间；
	diff = dvector(1, m);
	xo   = dvector(1, m);
	xi   = dvector(1, m);
	yo   = dvector(1, m);
	yi   = dvector(1, m);
	IDX  = ivector(1, m);

	for(i=0; i<N; i++)
	{
		xo[i] = (double)(contour_x[i]);
		yo[i] = (double)(contour_y[i]);
	}

	ComputeDifference(N, diff, xo, yo);  //求轮廓点序列的差分图，并将其存于数组diff中；

	DefineIDX(N, IDX, diff, 0, min_dist);

	// remove the points which distance to neighbor points is shorter than min_dist

	k = 0;
	for(i=0; i<N; i++)
	{
		if(IDX[i] == 1)
		{
			xo[k] = xo[i];
			yo[k] = yo[i];
			k ++;
		}	
	}

	// shift array position by one bit.

	N = k;

	ComputeDifference(N, diff, xo, yo);
 
	DefineIDX(N, IDX, diff, 1, max_dist);//找出所有大于max_dist的点用IDX图标识

	// call function of interp1

	Interp_length = interp1(IDX, xo, N, xi);//大于max_dist的点要进行插值
	Interp_length = interp1(IDX, yo, N, yi);

	N = Interp_length;

	CopyDDArray(N, xi, xo);
	CopyDDArray(N, yi, yo);

	ComputeDifference(N, diff, xo, yo);

	max_d = FindMaxValue(N, diff);

	k = 4;
	while( max_d > max_dist )  //当轮廓点序列中的最大差分值比max_dist大时，需要继续插值；
	{
		if(Interp_length>3000)
		{
			max_d=0;
			break;
		}
		DefineIDX(N, IDX, diff, 1, max_dist);  //定义插值标识；

		Interp_length = interp1(IDX, xo, N, xi);   //x方向线性插值两次；
		Interp_length = interp1(IDX, yo, N, yi);   //y方向线性插值两次；

		N = Interp_length;

		CopyDDArray(N, xi, xo);
		CopyDDArray(N, yi, yo);

		ComputeDifference(N, diff, xo, yo);  //计算差分图；
	
		max_d = FindMaxValue(N, diff);      //找差分图中最大值；
		k ++;
	}

	CopyDDArray(N, xo, contour_x);
	CopyDDArray(N, yo, contour_y);
	contour_x[N] = 0;
	contour_y[N] = 0;
	contour_length = N;
    //释放空间；
	free_dvector(diff,1,m);
	free_dvector(xo,1,m);
	free_dvector(xi,1,m);
	free_dvector(yo,1,m);
	free_dvector(yi,1,m);
	free_ivector(IDX,1,m);

	return( contour_length );
}
//输出：array_x和array_y，其他参数均为输入；
//功能：snake形变
void SnakeDeform(double* array_x, double* array_y, int array_length,  //轮廓点x、y数组以及轮廓的长度；
				double alpha_single, double beta_single,   //固定snake形变参数；
				double gamma_single, double kappa_single,   //固定snake形变参数
				double* OutputU_buf, double* OutputV_buf,   //U、V场
				int rows, int cols, int ITER   //图像的行和列，以及迭代的次数；
			   )//snake形变
{
	int    i, j, k, m, N;
	int    pixel_count;
	int    *indx;
	double dv;
	double *d_buf1, *d_buf2;
	double *alpha, *beta, *gamma, *kappa;
	double *alpham1, *alphap1;
	double *betam1, *betap1;
	double *a, *b, *c, *d, *e;
	double *Avec, *Gvec;
	double *col;
	double *vfx, *vfy;
	double *tmpX, *tmpY;
	PGdouble **Amatrix;
	PGdouble **InvAmat;
 
	N = array_length;
	m = 4 * N;
    
	indx = ivector(1, m);

	col   = dvector(1, m);
	alpha = dvector(1, m);
	beta  = dvector(1, m);
	gamma = dvector(1, m);
	kappa = dvector(1, m);

	alpham1 = dvector(1, m);
	alphap1 = dvector(1, m);
	betam1  = dvector(1, m);
	betap1  = dvector(1, m);

	a = dvector(1, m);
	b = dvector(1, m);
	c = dvector(1, m);
	d = dvector(1, m);
	e = dvector(1, m);

	tmpX = dvector(1, m);
	tmpY = dvector(1, m);

	vfx = dvector(1, m);
	vfy = dvector(1, m);

	Avec = dvector(N, N);
	Gvec = dvector(N, N);

	Amatrix = pgDmatrix(1, N, 1, N);
	InvAmat = pgDmatrix(1, N, 1, N);///

	for(i=0; i<N; i++)
	{
		alpha[i] = alpha_single;
		beta[i]  = beta_single;
	}

	// produce five diagnal vectors计算五对角阵A

	for(i=1; i<N; i++)
	{
		alpham1[i-1] = alpha[i];
	}
	alpham1[N-1] = alpha[0];

	alphap1[0] = alpha[N-1];
	for(i=1; i<N; i++)
	{
		alphap1[i] = alpha[i-1];
	}

	for(i=1; i<N; i++)
	{
		betam1[i-1] = beta[i];
	}
	betam1[N-1] = beta[0];

	betap1[0] = beta[N-1];
	for(i=1; i<N; i++)
	{
		betap1[i] = beta[i-1];
	}

	for(i=0; i<N; i++)
		a[i] = betam1[i];//beta[i+1]=beta

	for(i=0; i<N; i++)
		b[i] = -alpha[i] - 2.0 * beta[i] - 2.0 * betam1[i];//-alpha[i] - 2.0 *( beta[i] + beta[i+1])=-alpha - 4.0 *beta 

	for(i=0; i<N; i++)
		c[i] = alpha[i] + alphap1[i] + betam1[i] + 4.0 * beta[i] + betap1[i];//alpha[i] + alpha[i-1] + beta[i+1] + 4.0 * beta[i] + beta[i-1]=2.0*alpha  +6.0 * beta

	for(i=0; i<N; i++)
		d[i] = -alphap1[i] - 2.0 * beta[i] - 2.0 * betap1[i];//-alpha[i-1] - 2.0 * (beta[i] + beta[i-1])=-alpha - 4.0 * beta

	for(i=0; i<N; i++)
		e[i] = betap1[i];//beta[i-1]=beta

	// generate the parameters matrix "Amatrix"求得N*N阶五对角矩阵，注意矩阵元素边界处理
    //注意这种编程方式，他的时间复杂度为O(n)，而一次二重循环做的话时间复杂度为O(n*n)，但代码简洁
	pixel_count = N * N;
	d_buf1 = Avec;

	while( pixel_count -- )
		*(d_buf1 ++) = 0.0;

	m = N + 1;

	// A = diag(a(1:N-2), -2) + diag(a(N-1:N), N-2);

	k = N - 2;
	d_buf1 = Avec + 2 * N;
	for(i=0; i<k; i++)
	{
		*d_buf1 = a[i];
		d_buf1 += m;
	}

	d_buf1 = Avec + k;
	for(i=0; i<2; i++)
	{
		*d_buf1 = a[k+i];    //???
		d_buf1 += m;
	}

	// A = A + diag(b(1:N-1), -1) + diag(b(N), N-1);

	k = N - 1;
	d_buf1 = Avec + N;
	for(i=0; i<k; i++)
	{
		*d_buf1 += b[i];
		d_buf1 += m;
	}

	d_buf1 = Avec + k;
	*d_buf1 += b[k];

	// A = A + diag(c)

	d_buf1 = Avec;
	for(i=0; i<N; i++)
	{
		*d_buf1 += c[i];
		d_buf1 += m;
	}

	// A = A + diag(d(1:N-1), 1) + diag(d(N), -(N-1));

	k = N - 1;
	d_buf1 = Avec + 1;
	for(i=0; i<k; i++)
	{
		*d_buf1 += d[i];
		d_buf1 += m;
	}

	d_buf1 = Avec + k * N;
	*d_buf1 += d[k];

	// A = A + diag(e(1:N-2), 2) + diag(e(N-1:N), -(N-2));

	k = N - 2;
	d_buf1 = Avec + 2;
	for(i=0; i<k; i++)
	{
		*d_buf1 += e[i];
		d_buf1 += m;
	}

	d_buf1 = Avec + k * N;
	for(i=0; i<2; i++)
	{
		*d_buf1 += e[k+i];
		d_buf1 += m;
	}

	// invAI = inv(A + gamma * diag(ones(1,N)));求(A+gamma*I)的逆

	pixel_count = N * N;
	d_buf2 = Gvec;

	while( pixel_count -- )
		*d_buf2 ++ = 0.0;

	m = N + 1;
	d_buf2 = Gvec;
	for(i=0; i<N; i++)
	{
		*d_buf2 = gamma_single;
		d_buf2 += m;
	}

	d_buf1 = Avec;
	d_buf2 = Gvec;
	pixel_count = N * N;

	while( pixel_count -- )
	{
		*d_buf1 += *d_buf2;
		d_buf1 ++;
		d_buf2 ++;
	}

		// create a matrix

	d_buf1 = Avec;
	for(j=1; j<=N; j++)
	{
		for(i=1; i<=N; i++)
		{
			Amatrix[j][i] = *d_buf1 ++;
		}
	}

	// inverse a matrix

	ludcmp(Amatrix, N, indx, &dv);
	for(j=1; j<=N; j++)
	{
		for(i=1; i<=N; i++)
			col[i] = 0.0;
		col[j] = 1.0;

		lubksb(Amatrix, N, indx, col);

		for(i=1; i<=N; i++)
			InvAmat[i][j] = col[i];
	}//end of 求(A+gamma*I)的逆

	for(k=0; k<ITER; k++)
	{
		interp2(OutputU_buf, array_x, array_y, vfx, N, rows, cols);//将外力场x方向插值
		interp2(OutputV_buf, array_x, array_y, vfy, N, rows, cols);//将外力场y方向插值

		// deform snake
		
		for(i=0; i<N; i++)
		{
			tmpX[i] = gamma_single * array_x[i] + kappa_single * vfx[i];
			tmpY[i] = gamma_single * array_y[i] + kappa_single * vfy[i];
		}

		for(j=1; j<=N; j++)
		{
			dv = 0.0;
			for(i=1; i<=N; i++)
			{
				dv += (InvAmat[j][i] * tmpX[i-1]);
			}
			array_x[j-1] = dv;
		}

		for(j=1; j<=N; j++)
		{
			dv = 0.0;
			for(i=1; i<=N; i++)
			{
				dv += (InvAmat[j][i] * tmpY[i-1]);
			}
			array_y[j-1] = dv;
		}
	}
    //释放空间；
	free_ivector(indx,1,m);
	free_dvector(col,1,m);
	free_dvector(alpha,1,m);
	free_dvector(beta,1,m);
	free_dvector(gamma,1,m);
	free_dvector(kappa,1,m);

	free_dvector(alpham1,1,m);
	free_dvector(alphap1,1,m);
	free_dvector(betap1,1,m);
	free_dvector(betam1,1,m);

	free_dvector(a,1,m);
	free_dvector(b,1,m);
	free_dvector(c,1,m);
	free_dvector(d,1,m);
	free_dvector(e,1,m);
	free_dvector(vfx,1,m);
	free_dvector(vfy,1,m);
	free_dvector(tmpX,1,m);
	free_dvector(tmpY,1,m);

	free_dvector(Avec,N,N);
	free_dvector(Gvec,N,N);
   
	pgFreeDmatrix(Amatrix,1,N,1,N);
	pgFreeDmatrix(InvAmat,1,N,1,N);

	return;
}

// ----------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------

// Related computing subroutines for active contour computing
/*
 ************************************************************************************
 *																					*
 *		This function is used to record contour array (x, y) in order				*
 *																					*
 ************************************************************************************
 */

int record_contour_array(USHORT* contour_buf, // image buffer contains contour results generateed by "Define_RegionContour"  
						 int rows, int cols,  // rows and columes of image buffer
						 int* x_array, int* y_array		// two (x, y) arrays contain contour pixels
						)//将曲线上的点按行走路线有序存储
{
	USHORT *r_buf1;
	int    i, j, k;
	int    ib, jb;
	int    dd, min_dd, min_i;
	int	   contour_count;
	int    srow, scol;

	srow = rows; 
	scol = cols;

	r_buf1 = contour_buf;
	k = 0;
	for(j=0; j<srow; j++)
	{
		for(i=0; i<scol; i++)
		{
			if(*r_buf1 == MAX_12_BIT)
			{
				x_array[k] = i;     //初始化轮廓数组；
				y_array[k] = j;
				k ++;
			}
			r_buf1 ++;
		}
	}
	contour_count = k;

	// reorganize contour arrays

	for(k=0; k<contour_count-1; k++)
	{
		ib = x_array[k];
		jb = y_array[k];
		min_dd = 4000;
		for(i=k; i<contour_count; i++)
		{
			dd = (ib - x_array[i]) * (ib - x_array[i]) + (jb - y_array[i]) * (jb - y_array[i]);
			if(dd < min_dd && dd > 0)
			{
				min_dd = dd;
				min_i = i;
			}
		}
		i = x_array[k+1];
		j = y_array[k+1];
		x_array[k+1] = x_array[min_i];
		y_array[k+1] = y_array[min_i];
		x_array[min_i] = i;
		y_array[min_i] = j;
	}

	return( contour_count);
}  // record_contour_array()

/*
 ************************************************************************************
 *																					*
 *		This function is used to compute active contour				*
 *																					*
 ************************************************************************************
 */
int Computing_ActiveContour(USHORT* img_buf,	  // image buffer of original image
							USHORT* growth_buf,  // image buffer recorded the original growth region (binary image) 
							int rows, int cols,  // rows and columns of image buffer
							int orig_contour_count,   // count of contour pixels
							int* x_array, int* y_array  // x, y array of contour pixels
						   )//计算GVF场并计算活动轮廓
{
	USHORT *r_buf1, *r_buf2;
	USHORT *EdgeMap_buf;
	int    i, j, k, n;
	int    ITER;
	int    pixel_count;
	int    contour_count;
	double mu;        //固定参数 mu=0.1;
	double u, v, mag;
	double max_dist, min_dist;    //snakeInterp参数，相邻两点最大的距离max_dist和最小距离min_dist；
	double alpha, beta, gamma, kappa;    //固定sake形变参数；
	double *FMap_buf;
	double *OutputU_buf, *OutputV_buf;   
	double *px_buf, *py_buf;
	double *contour_x, *contour_y;
	double *f_buf1, *f_buf2, *f_buf3, *f_buf4;
    //申请空间；
	contour_count = orig_contour_count;
	EdgeMap_buf = svector(rows, cols);
	FMap_buf    = dvector(rows, cols);
	OutputU_buf = dvector(rows, cols);
	OutputV_buf = dvector(rows, cols);
	px_buf      = dvector(rows, cols);
	py_buf      = dvector(rows, cols);
	contour_x   = dvector(rows*cols/3, 1);//??
	contour_y   = dvector(rows*cols/3, 1);

	copy_svector(img_buf, EdgeMap_buf, rows, cols);


	// read initial contour data array

	for(k=0; k<contour_count; k++)
	{
		contour_x[k] = (double)(x_array[k]);
		contour_y[k] = (double)(y_array[k]);
	}
	contour_x[contour_count] = 0.0;
	contour_y[contour_count] = 0.0;

	//计算EdgeMap_buf的梯度图置入EdgeMap_buf
	USHORT  * Para_a, * Para_b, * Para_c, * OrigImg;    //求最小二乘法梯度图用到的变量；
	Para_a = svector(rows, cols);
	Para_b = svector(rows, cols);
	Para_c = svector(rows, cols);
	OrigImg = svector(rows, cols);
	copy_svector(EdgeMap_buf, OrigImg, rows, cols);	
	MinQuare(Para_a,Para_b,Para_c,EdgeMap_buf,OrigImg);   //最小二乘法梯度；
 
	//测试
			IplImage* temp = cvCreateImage(cvSize(125,125), 8, 1);
			int stepsPL= temp->widthStep;
			BYTE * pData=(BYTE*)temp->imageData;
			for(i=0;i<125;i++)
				for(j=0;j<125;j++)
				{
					int jyx;
					jyx=(*(EdgeMap_buf+125*i+j))/16;
					if(jyx>255)
					pData[i*stepsPL+j]=255;
					else
					pData[i*stepsPL+j]=jyx;
				}
/*
        cvNamedWindow("jyx",1);
        cvShowImage("jyx",temp);      //显示梯度图；
        AfxMessageBox("nihao",MB_OK);*/

	
   // gradient(EdgeMap_buf,rows,cols);   //差分梯度；
 				
    // remove pixels inside growth region		
	pixel_count = rows * cols;
	r_buf1 = growth_buf;
	r_buf2 = EdgeMap_buf;
	while( pixel_count-- )
	{
        if(*r_buf1 == MAX_12_BIT)
			*r_buf2 = 0;
			
		r_buf1 ++;
		r_buf2 ++;
	}
	pixel_count = rows * cols;
	f_buf1 = FMap_buf;
	r_buf1 = EdgeMap_buf;
	while( pixel_count -- )
	{
		*f_buf1 ++ = (double)(*r_buf1 ++); 
	}

	// compute the GVF of the edge map

	mu = 0.1;
	ITER = 80;  

	GVFC(rows, cols, FMap_buf, OutputU_buf, OutputV_buf, mu, ITER);	
	// normalize the GVF external force.
    
	f_buf1 = OutputU_buf;
	f_buf2 = OutputV_buf;

	f_buf3 = px_buf;
	f_buf4 = py_buf;

	pixel_count = rows * cols;
    double minva,maxva;
    minva=1.0;
	maxva=-1.0;
	while( pixel_count -- )
	{
		u = *f_buf1;
		v = *f_buf2;

		mag = sqrt(u * u + v * v) + 1e-10;

		*f_buf3 = u / mag;
		*f_buf4 = v / mag;
        minva=(*f_buf3)<minva?(*f_buf3):minva;
		maxva=(*f_buf3)>maxva?(*f_buf3):maxva;
		f_buf1 ++;
		f_buf2 ++;
		f_buf3 ++;
		f_buf4 ++;
	}

	//测试；
    for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			pData[i*stepsPL+j]=(BYTE)(((double)255*((*(px_buf+125*i+j))-minva)/(maxva-minva))+0.5);
		//	pData[i*stepsPL+j]=(BYTE)((double)255*(*(px_buf+125*i+j)));
		}

/*
        cvNamedWindow("jyx",1);
        cvShowImage("jyx",temp);
        AfxMessageBox("nihao",MB_OK);*/
    
	// snake deformation

	max_dist = 2.0;  
	min_dist = 0.5;
	k = SnakeInterp(contour_x, contour_y, contour_count, max_dist, min_dist);   //插值

	contour_count = k;

	// iteration for snake deformation and interpretation.
    
	ITER = 5;
	alpha = 0.05;
	beta  = 0.0;
	gamma = 1.0;
	kappa = 0.6;

	int snake_iter =20;   //snack_iter =9；迭代次数多些没关系；

	for(k=1; k<snake_iter; k++)//插值与否没有效果上的改变，只是效果不同
	{
		if(contour_count>1&&contour_count<3000)  
		{
			SnakeDeform(contour_x, contour_y, contour_count, alpha, beta, gamma, kappa, 
				        px_buf, py_buf, rows, cols, ITER);

		    j = SnakeInterp(contour_x, contour_y, contour_count, max_dist, min_dist);
		    contour_count = j;
		}
	}
	max_dist=1.0;
	min_dist=0.5;
	if(contour_count<3000)
	{
        j = SnakeInterp(contour_x, contour_y, contour_count, max_dist, min_dist);
	    contour_count = j;
	}
	// record the detected contour cordinate (x, y) into x_array and y_array
	pixel_count = 0;
	for(k=0; k<contour_count; k++)
	{
		i = (int)(0.5 + contour_x[k]);
		j = (int)(0.5 + contour_y[k]);
		if(i<0||j<0||i>124||j>124)
		{
			continue;
		}
		n = *(growth_buf+j*cols+i);//用到ActiveContourRegion函数返回的growth_buf值
		if(n == 0)
		{
			x_array[pixel_count] = (int)(0.5 + contour_x[k]);
			y_array[pixel_count] = (int)(0.5 + contour_y[k]);
			pixel_count ++;
		}
	}
	contour_count = pixel_count;
	x_array[contour_count] = 0;
	y_array[contour_count] = 0;
    //释放空间；
	free_svector(EdgeMap_buf);
	free_dvector(FMap_buf,rows,cols);
	free_dvector(OutputU_buf,rows,cols);
	free_dvector(OutputV_buf,rows,cols);
	free_dvector(px_buf,rows,cols);
	free_dvector(py_buf,rows,cols);
	free_dvector(contour_x,rows*cols/3, 1);
	free_dvector(contour_y,rows*cols/3, 1);
	free_svector(Para_a);
	free_svector(Para_b);
	free_svector(Para_c);
	free_svector(OrigImg);
	cvReleaseImage(&temp);

	return( contour_count );
}  // End of Computing_ActiveContour()

//输入：rows，scol；contour_count；x_array，y_array
//输出：growth_buf
//功能：计算包含肿块区域的矩形，该矩形的功能实际上是限制梯度图和GVF场的范围；
void ActiveContourRegion(USHORT* growth_buf,	// image buffer of growth region defined by active contour
						 int rows, int cols,    // rows and columns of image buffer
						 int contour_count,		// number of pixels in boundary contour
						 int* x_array, int* y_array  // contour (x, y) arrays
						)//将growth_buf赋值
{
	USHORT *LocalGrowth_buf;
	USHORT *LocalMorph_buf;
	USHORT *r_buf1;

	int    i, j, k, pixel_count;
	int    contour_pixels;
	int    cent_x, cent_y;
	int    max_x, max_y, min_x, min_y;
	int    scol, srow;
	double *contour_x, *contour_y;
    contour_x   = dvector(rows*cols/3, 1);//??
	contour_y   = dvector(rows*cols/3, 1);
	for(k=0; k<contour_count; k++)
	{
		contour_x[k] = (double)(x_array[k]);
		contour_y[k] = (double)(y_array[k]);
	}
	contour_x[contour_count] = 0.0;
	contour_y[contour_count] = 0.0;
	scol = cols;
	srow = rows;

	LocalGrowth_buf = svector(srow, scol);
	LocalMorph_buf  = svector(srow, scol);

	pixel_count = srow * scol;
	r_buf1 = LocalGrowth_buf;
	while( pixel_count -- )
	{
		*r_buf1 ++ = 0;
	}
	for(k=0; k<contour_count; k++)
	{
		i = contour_x[k];
		j = contour_y[k];
		*(LocalGrowth_buf+j*scol+i) = MAX_12_BIT;
	}

	// Find contour center
	contour_pixels = 0;
	cent_x = 0;
	cent_y = 0;
	min_x = 4000;
	min_y = 4000;
	max_x = 0;
	max_y = 0;
	r_buf1 = LocalGrowth_buf;
	for(j=0; j<srow-2; j++)
	{
		for(i=0; i<scol; i++)
		{
			if(*r_buf1 == MAX_12_BIT)
			{
				cent_x += i;
				cent_y += j;
				contour_pixels ++;

				if(i < min_x)
					min_x = i;
				if(i > max_x)
					max_x = i;
				if(j < min_y)
					min_y = j;
				if(j > max_y)
					max_y = j;
			}
			r_buf1 ++;
		}
	}
    //将包含肿块的最小矩形向外扩充20个像素点，得到的新矩形区域为受限肿块区域；
	if(contour_pixels > 0)
	{
		cent_x /= contour_pixels;   // (cent_x,cent_y)为肿块的几何中心；
		cent_y /= contour_pixels;

		min_x -= 20;
		if(min_x < 0)
			min_x = 0;

		max_x += 20;
		if(max_x >= scol)
			max_x = scol - 1;

		min_y -= 20;
		if(min_y < 0)
			min_y = 0;

		max_y += 20;
		if(max_y >= srow)
			max_y = srow - 1;
	}

	// 将受限肿块区域之外的区域像素点灰度值设为4095；受限肿块区域灰度值设为0；

	for(j=0; j<min_y; j++)
		for(i=0; i<scol; i++)
			*(LocalGrowth_buf+j*scol+i) = 4095;

	for(j=max_y; j<srow; j++)
		for(i=0; i<scol; i++)
			*(LocalGrowth_buf+j*scol+i) = 4095;

	for(i=0; i<min_x; i++)
		for(j=0; j<srow; j++)
			*(LocalGrowth_buf+j*scol+i) = 4095;

	for(i=max_x; i<scol; i++)
		for(j=0; j<srow; j++)
			*(LocalGrowth_buf+j*scol+i) = 4095;
	for(k=0; k<contour_count; k++)
	{
		i = contour_x[k];
		j = contour_y[k];
		*(LocalGrowth_buf+j*scol+i) =0;
	}
					
	// save the growth suspicious region

	copy_svector(LocalGrowth_buf, growth_buf, srow, scol);  //将growth_buf赋值

	free_svector(LocalGrowth_buf);
	free_svector(LocalMorph_buf);

	return;
}
//Compute the map gradient
void gradient(USHORT*f ,	// image buffer of growth region defined by active contour
			    int rows, int cols)
{
	int i,j;
	PGdouble **fx, **fy;
    fx = pgDmatrix(0,rows-1,0,cols-1);
    fy = pgDmatrix(0,rows-1,0,cols-1);
	
    //处理corner
	fx[0][0]=fx[0][cols-1]=fx[rows-1][0]=fx[rows-1][cols-1]=0;
	fy[0][0]=fy[0][cols-1]=fy[rows-1][0]=fy[rows-1][cols-1]=0;
	//处理right和left列
	for (i=1; i<rows-1; i++)
	{
		fx[i][0]=f[i*cols+1]-f[i*cols];
		fx[i][cols-1]=f[i*cols+cols-1]-f[i*cols+cols-2];
        fy[i][0]=f[(i+1)*cols]-f[i*cols];
		fy[i][cols-1]=f[(i+1)*cols+cols-1]-f[i*cols+cols-1];
		
	}
	//处理bottom和top行
	for(j=1; j<cols-1; j++)
	{
		fx[0][j]=f[j+1]-f[j];
		fx[rows-1][j]=f[(rows-1)*cols+j+1]-f[(rows-1)*cols+j];
		fy[0][j]=f[cols+j]-f[j];
		fy[rows-1][j]=f[(rows-1)*cols+j]-f[(rows-2)*cols+j];
	}
	//处理内部梯度
	for(i=1; i<rows-1; i++)
		for(j=1; j<cols-1; j++)
		{
			fx[i][j]=(f[i*cols+j+1]-f[i*cols+j-1])/2;
			fy[i][j]=(f[(i+1)*cols+j]-f[(i-1)*cols+j])/2;
		}
    USHORT *pTemp;
	pTemp=f;
	for(i=1; i<rows-1; i++)
		for(j=1; j<cols-1; j++) 
		{
			pTemp = f+i*cols+j;
			*pTemp=sqrt(fx[i][j]*fx[i][j]+fy[i][j]*fy[i][j]);
		}

   pgFreeDmatrix(fx,0,rows-1,0,cols-1);
   pgFreeDmatrix(fy,0,rows-1,0,cols-1);
	  
}

//输入：肿块区域图； 
//输出和返回：种子点；
CvPoint  FindSeed(IplImage *flagImg)
{
	int i;
	int steps_PL = flagImg->widthStep;
	BYTE* pFlagData = (BYTE*)flagImg->imageData;
	CvPoint candi[4] = {0, 0, flagImg->width-1, 0, 0, flagImg->height-1, flagImg->width-1, flagImg->height-1};
	for(i=0; i<4; i++)
		if(pFlagData[candi[i].y*steps_PL+candi[i].x] != 255)
			return candi[i];
}

//根据边界点修改区域标志：0---背景，178---区域内像素，255---边界像素
void  Cal_RegionFlag(IplImage* flagImg)//flagImg: 8U
{
	using namespace std;
	IplImage* tempImg = cvCloneImage(flagImg);
	int i, j;
	int steps_PL = flagImg->widthStep;
	BYTE* pTempData = (BYTE*)tempImg->imageData;
	BYTE* pFlagData = (BYTE*)flagImg->imageData;

	CvPoint seed = FindSeed(flagImg);
	bool** visited = new bool*[flagImg->height];
	for(i=0; i<flagImg->height; i++)
		visited[i] = new bool[flagImg->width];
	for(i=0; i<flagImg->height; i++)
		for(j=0; j<flagImg->width; j++)
		{
			visited[i][j] = false;
		}

	//生长
	list<CvPoint> stack;
	stack.push_back(seed);
	pFlagData[seed.y* steps_PL + seed.x] = 0;
	visited[seed.y][seed.x] = true;
	int dir[4][2] = {{1, 0}, {0, 1}, {0, -1}, {-1, 0}};
	while(!stack.empty())
	{
		CvPoint curPoint = stack.back();//取出
		stack.pop_back();//弹出
		for(i=0; i<4; i++)
		{
			CvPoint nextPoint = cvPoint(curPoint.x + dir[i][0], curPoint.y + dir[i][1]);
			if(nextPoint.x < 0 || nextPoint.y <0
			|| nextPoint.x > flagImg->width-1 || nextPoint.y > flagImg->height-1)
				continue;
			if(visited[nextPoint.y][nextPoint.x])
				continue;
			if(pTempData[nextPoint.y* steps_PL + nextPoint.x] == 255)//遇到边界点
				continue;
			pFlagData[nextPoint.y* steps_PL + nextPoint.x] = 0;
			stack.push_back(nextPoint);
			visited[nextPoint.y][nextPoint.x] = true;
		}
	}
	for(i=0; i<flagImg->height; i++)
		delete[] visited[i];
	delete[] visited;
	cvReleaseImage(&tempImg);
}

//实现平滑；
void Gaussion(IplImage *img)
{

	int iTempW=3,iTempH=3;
	int iTempMY=1,iTempMX=1;
	float fCoef = 0.065;
	float fpArray[]={1,2,1,2,4,2,1,2,1};    //平滑模板；
	float fResult;
	BYTE * m_pData=(BYTE *)img->imageData;
	int steps_PL = img->widthStep * 8 / img->depth;
	int m_height=img->height;
	int m_width=img->width;
	BYTE **pTemp=new BYTE*[m_height];
	for (int i=0;i<m_height;i++)
	{
		pTemp[i]=new BYTE[m_width];
		for(int j=0;j<m_width;j++)
			pTemp[i][j]=m_pData[steps_PL*i+j];
	}		
	for( i=iTempMY;i<m_height+iTempMY+1-iTempH;i++)//
		
	{  
		for(int j=iTempMX;j<m_width-iTempW+iTempMX+1;j++)
		{
			//指向新dib第i行，第j个象素的指针			
			fResult=0;
			for(int k=0;k<iTempH;k++)
			{
				for(int l=0;l<iTempW;l++)
				{		

					fResult+=(pTemp[m_height-1-i+iTempMY-k][j-iTempMX+1])*fpArray[k*iTempW+l];
				}
			}
			fResult*=fCoef;
			fResult=(float)fabs(fResult);
			if(fResult>255)
			m_pData[steps_PL*(m_height-1-i)+j]=255;
			else	
			m_pData[steps_PL*(m_height-1-i)+j]=(BYTE)(fResult+0.5);
			
		}
	}

	for(i=0;i<m_height;i++)
		delete []pTemp[i];
	delete []pTemp;
	pTemp=NULL;
		 
}

/****输入：梯度图（边缘图），输出和返回：圆形初始轮廓的半径；
     算法：以肿块中心为中心，求0、45、90、135、180、225、270、315这八个方向射线上的
	 梯度值最大位置。******/
int  FindMassRadius(USHORT * EdgeMap)
{
	int  radius;
	int  maxg1,maxg2,maxg3,maxg4,maxg5,maxg6,maxg7,maxg8;    //定义八个方向梯度的最大值变量；
	maxg1=maxg2=maxg3=maxg4=maxg5=maxg6=maxg7=maxg8=0;   
	CPoint  pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8;    //定义八个方向梯度值为最大的位置变量；
	pos1.x=pos2.x=pos3.x=pos4.x=pos5.x=pos6.x=pos7.x=pos8.x=63;
	pos1.y=pos2.y=pos3.y=pos4.y=pos5.y=pos6.y=pos7.y=pos8.y=63;
	CPoint  centre(63,63);
	int i,j;
	for(i=0;i<9;i++)
	{
		int  sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8;
		sum1=sum2=sum3=sum4=sum5=sum6=sum7=sum8=0;
		int step;
		for(j=0;j<5;j++)    //求连续五个点的梯度值和；
		{
			step=5*i+j;
			sum1=sum1+*(EdgeMap+125*centre.x+(centre.y+step));
			sum2=sum2+*(EdgeMap+125*(centre.x-step)+(centre.y+step));
			sum3=sum3+*(EdgeMap+125*(centre.x-step)+centre.y);
			sum4=sum4+*(EdgeMap+125*(centre.x-step)+(centre.y-step));
			sum5=sum5+*(EdgeMap+125*centre.x+(centre.y-step));
			sum6=sum6+*(EdgeMap+125*(centre.x+step)+(centre.y-step));
			sum7=sum7+*(EdgeMap+125*(centre.x+step)+centre.y);
			sum8=sum8+*(EdgeMap+125*(centre.x+step)+(centre.y+step));
		}
		//取平均；
		sum1=sum1/5;
		sum2=sum2/5;
		sum3=sum3/5;
		sum4=sum4/5;
		sum5=sum5/5;
		sum6=sum6/5;
		sum7=sum7/5;
		sum8=sum8/5;
		//更新八个方向最大的梯度值以及其位置；
		if(sum1>maxg1)
		{
			maxg1=sum1;
			pos1.y=centre.y+(5*i+2);
		}
		if(sum2>maxg2)
		{
			maxg2=sum2;
			pos2.x=centre.x-(5*i+2);
			pos2.y=centre.y+(5*i+2);
		}
		if(sum3>maxg3)
		{
			maxg3=sum3;
			pos3.x=centre.x-(5*i+2);
		}
		if(sum4>maxg4)
		{
			maxg4=sum4;
			pos4.x=centre.x-(5*i+2);
			pos4.y=centre.y-(5*i+2);
		}
		if(sum5>maxg5)
		{
			maxg5=sum5;
			pos5.y=centre.y-(5*i+2);
		}
		if(sum6>maxg6)
		{
			maxg6=sum6;
			pos6.x=centre.x+(5*i+2);
			pos6.y=centre.y-(5*i+2);
		}
		if(sum7>maxg7)
		{
			maxg7=sum7;
			pos7.x=centre.x+(5*i+2);
		}
		if(sum8>maxg8)
		{
			maxg8=sum8;
			pos8.x=centre.x+(5*i+2);
			pos8.y=centre.y+(5*i+2);
		}
	}
	//取八个方向上的半径的平均值为最后求得的半径；
radius=0;
radius=radius+sqrt((pos1.x-centre.x)*(pos1.x-centre.x)+(pos1.y-centre.y)*(pos1.y-centre.y))
             +sqrt((pos2.x-centre.x)*(pos2.x-centre.x)+(pos2.y-centre.y)*(pos2.y-centre.y))
			 +sqrt((pos3.x-centre.x)*(pos3.x-centre.x)+(pos3.y-centre.y)*(pos3.y-centre.y))
			 +sqrt((pos4.x-centre.x)*(pos4.x-centre.x)+(pos4.y-centre.y)*(pos4.y-centre.y))
			 +sqrt((pos5.x-centre.x)*(pos5.x-centre.x)+(pos5.y-centre.y)*(pos5.y-centre.y))
			 +sqrt((pos6.x-centre.x)*(pos6.x-centre.x)+(pos6.y-centre.y)*(pos6.y-centre.y))
			 +sqrt((pos7.x-centre.x)*(pos7.x-centre.x)+(pos7.y-centre.y)*(pos7.y-centre.y))
			 +sqrt((pos8.x-centre.x)*(pos8.x-centre.x)+(pos8.y-centre.y)*(pos8.y-centre.y));
radius=radius/8+2;
assert(radius<55); 
return  radius;
}
//输入：轮廓半径；  输出：记录轮廓的缓冲区，轮廓点序列数组；
int  GetContourPoint(int radius , USHORT * Contour_buf,int * x_array,int * y_array)
{
	const int centre=63;
	CPoint lefttop,rightbottom;
    lefttop.x=lefttop.y=63-radius;
	rightbottom.x=rightbottom.y=63+radius;
	int i,j;
	int dist;
	int index=0;
	//将到中心的距离为radius的点记为轮廓点；
	for(i=lefttop.x;i<=rightbottom.x;i++)
		for(j=lefttop.y;j<=rightbottom.y;j++)
		{
			dist=sqrt((i-centre)*(i-centre)+(j-centre)*(j-centre));
			if(dist==radius)
			{
				*(Contour_buf+i*125+j)=4095;
				x_array[index]=i;
				y_array[index]=j;
				index++;
			}
		}
	return  index;
}
/****计算轮廓点数和肿块的面积，若轮廓点数和肿块的面积大于特定的阈值，就要重新获得初始轮廓，
     若重新定初始轮廓，返回true,否则返回false ***/
bool  IsGetRoundContour(IplImage * RoiRegion)
{
    const int cont_threshold=200;   //定义轮廓点阈值为200；
	const int area_threshold=1200;  //定义肿块面积阈值为1200；
	int steps_PL = RoiRegion->widthStep * 8 / RoiRegion->depth;
	BYTE * pData = (BYTE*) RoiRegion->imageData;
	int i,j;
	//求肿块的轮廓点数和肿块的面积；
	int contnum,areanum;
	contnum=areanum=0; 
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			if(pData[i*steps_PL+j]==255)
				contnum++;
			else if(pData[i*steps_PL+j]==178)
				areanum++;
			else 
				;
		}
	if(contnum>=cont_threshold&&areanum>=area_threshold)   //判断；
		return   true;
	else 
		return   false;
}
//输入：img
//输出：resultImg
//功能：最小二乘法背景趋势去除；
void  BackgrdCorrectMatlab(IplImage *img,IplImage * resultImg)
{
	cvSetZero(resultImg);
	//求最小二乘法拟和的平面参数
	int i,j;
	double a,b,c;//参数
	int steps_PL = img->widthStep * 8 / img->depth;
	BYTE * pData = (BYTE*) img->imageData;
	BYTE* pDest = (BYTE*) resultImg->imageData;
	int nRows,nCols;
	nRows = nCols = img->width;
	double sumII,sumJJ,sumIJ,sumI,sumJ,sum1;
	sumII=sumJJ=sumIJ=sumI=sumJ=sum1=0;
	double sumIPij,sumJPij,sumPij;
	sumIPij=sumJPij=sumPij=0;
	for(i=0;i<nRows;i++)
	{
		for(j=0;j<nCols;j++)
		{
			sumII += i*i;
			sumJJ += j*j;
			sumIJ += i*j;
			sumI  += i;
			sumJ  += j;
			sum1  += 1;

			sumIPij += i * pData[i*steps_PL+j];
			sumJPij += j * pData[i*steps_PL+j];
			sumPij  += pData[i*steps_PL+j];
		}
	}
	double value[9]={
		sumII,sumIJ,sumI,
		sumIJ,sumJJ,sumJ,
		sumI ,sumJ ,sum1   
	};
	CMatrix matrix(3,3,value);
	BOOL flag = matrix.InvertSsgj(); 
	if(flag)
	{
		double valueTmp[3]={
			sumIPij,sumJPij,sumPij
		};
		CMatrix matrixTmp(3,1,valueTmp);
		CMatrix para(3,1);
		para = matrix*matrixTmp;
		a=para.GetElement(0,0);
		b=para.GetElement(1,0);
		c=para.GetElement(2,0);
	}
	//进行校正 ,保存成原始数据
	double** rawResult = new double*[nRows];
	for(i=0; i<nRows; i++)
		rawResult[i] = new double[nCols];
	double min = 255; 
	double max = 0;
	for(i=0;i<nRows;i++)
	{
		for(j=0;j<nCols;j++)
		{
			double temp = a*i + b*j + c;
			rawResult[i][j] = pData[i*steps_PL+j] - temp;
			if(rawResult[i][j] > max)
				max = rawResult[i][j];
			if(rawResult[i][j] < min)
				min = rawResult[i][j];
		}
	}
	//进行标定，保存成结果图像
	for(i=0;i<nRows;i++)
	{
		for(j=0;j<nCols;j++)
		{
			if(max == min)
			{
				pDest[i*steps_PL + j] = 0;
				continue;
			}
			double temp = 255 * ((rawResult[i][j] - min)/ (max - min));
			pDest[i*steps_PL + j] = (BYTE) temp;
		}
	}
}
void CalHistogram8U(const IplImage* PImage, int* histogram, int totalgraylevel, IplImage* mask)
{
	assert(PImage != NULL);
	assert(PImage->depth == 8);
	assert(PImage->nChannels == 1);
	int i, j;
	int width = PImage->width;
	int height = PImage->height;
	int bytesPL = 8 * PImage->widthStep / PImage->depth;//4字节对齐
	BYTE* pData = (BYTE*)PImage->imageData;
	BYTE* pMask ;
	if(mask != NULL)
		pMask = (BYTE*)mask->imageData;
	else
		pMask = NULL;
	for(i=0; i<totalgraylevel; i++)
		histogram[i] = 0;
	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
			int value = (int) pData[i*bytesPL + j];
			if(pMask == NULL)
				histogram[value] ++;
			else
				if(pMask[i*bytesPL + j] > 0)
					histogram[value] ++;
		}
	}
}
//实现Gama校正
IplImage * GamaCorrectMatlab(IplImage *image)
{	
	int i,j;
	const double  threshold=0.005;
	BYTE* pData = (BYTE*) image->imageData;
	int steps_PL = image->widthStep * 8 / image->depth;
	int height,width;
	height = image->height; 
	width = image->width;
	double maxValue,minValue;
	maxValue=pData[0];
	minValue=pData[0];
	int count = (int)(threshold * image->width * image->height);
	int* histogram = new int[256];
	CalHistogram8U(image, histogram, 256, NULL);
	int sumLo = 0;
	int sumHi = 0;
	for(i=0; i<255; i++)
	{
		sumLo += histogram[i];
		if(sumLo > count)
		{
			minValue = i;
			break;
		}
	}
	for(i=255; i>=0; i--)
	{
		sumHi += histogram[i];
		if(sumHi > count)
		{
			maxValue = i;
			break;
		}
	}
	
	//进行标定，保存成结果图像
	for(i=0;i<height;i++)
	{
		for(j=0;j<width;j++)
		{
			double rate = (pData[i*steps_PL + j]-minValue)/(maxValue-minValue);
			if(rate > 1)
				rate = 1;
			if(rate < 0)
				rate = 0;
			double temp = 255 * rate * rate ;	
			pData[i*steps_PL + j] = (BYTE) temp;
		}
	}
	delete histogram;
	return image;
}

double  MaxInterstify(int m,double a[])
{
	double m_pix1=0;
	for(int i=0;i<=m;i++)
	{
		m_pix1+=a[i];
	}
	double m_pix2=1-m_pix1;
	double h=0;
	for(int j=45;j<=m;j++)
	{
		if(a[j] > 0.000001)
			h-=a[j]*log(a[j]);
	}
	h=h/m_pix1;
	double h2=0;
	for(int t=m+1;t<256;t++)
	{
		if(a[t] > 0.000001)
			h2-=a[t]*log(a[t]);
	}
	h2=h2/m_pix2;
	h+=h2;
	h+=logf(m_pix1);
	h+=logf(m_pix2);
	return h;
}

//最大熵分割法
void  MaxEntropy(IplImage *image)
{	
	int  m_height,m_width;
	m_height=image->height;
	m_width=image->width;
	BYTE * m_pData=(BYTE *)image->imageData;
	int steps_PL = image->widthStep * 8 / image->depth;
	double m_Intensify[256]={0},sum_pix=0;	 
	for(int i=0;i<m_height;i++)
		for(int j=0;j<m_width;j++)
		{
			
			m_Intensify[m_pData[i*steps_PL + j]]++;
			sum_pix++;
		}
		double p[256]={0};
		for(int j=0;j<256;j++)
		{
			p[j]=double(m_Intensify[j]/sum_pix);
		}		
		double entropy[256]={0};
		for(j=50;j<150;j++)
			entropy[j]=MaxInterstify(j,p);
		int max_inter=1;
		double maxen=entropy[1];
		for(i=50;i<150;i++)
			if(entropy[i]>maxen)
			{
				maxen=entropy[i];
				max_inter=i;
			}
			for(i=0;i<m_height;i++)
				for(j=0;j<m_width;j++)
				{
					
					if(m_pData[i*steps_PL + j]>max_inter)
						m_pData[i*steps_PL + j]=255;
					else
						m_pData[i*steps_PL + j]=0;
				}		
				
}
//输入：原始图像 img
//输出和返回：经过背景趋势去除和最大熵方法得到的肿块区域图，返回的图像是一个二值图像，
//            肿块区域灰度值为255，背景区域灰度值为0；

IplImage * BackgroundTrendRemoveAndMaxEntropy(IplImage * img)
{
	IplImage *  resultimg;
	resultimg=cvCreateImage(cvSize(125, 125), 8, 1);
	BackgrdCorrectMatlab(img,resultimg);
	GamaCorrectMatlab(resultimg);   //Gama校正；
	Gaussion(resultimg);    //高斯滤波；
	MaxEntropy(resultimg);  //最大熵方法分割肿块；
	return resultimg;      //返回结果图像；
}

//输入：img
//输出：resultimg
//功能：实现背景趋势去除；
 void  BackgroundTrendRemove(IplImage * img,IplImage *resultimg)
{
	BackgrdCorrectMatlab(img,resultimg);
	GamaCorrectMatlab(resultimg);   //Gama校正；
}
/*******基于四邻域取肿块区域轮廓，令轮廓点灰度值为255，背景区域灰度为0，
        肿块内部点灰度为178 ******/
void  GetRegionContour(IplImage * ImgRegion)  //ImgRegion为8位区域图，
                                              //肿块轮廓和背景的灰度为0，肿块区域内部灰度为178。
{
	int stepsPL= ImgRegion->widthStep;
	BYTE * pData=(BYTE*)ImgRegion->imageData;
	int i,j;
	BYTE top,bottom,left,right;
	for(i=0;i<ImgRegion->height;i++)
		for(j=0;j<ImgRegion->width;j++)
		{
			if((pData[stepsPL*i+j]==0)&&i>1&&j>1&&i<(ImgRegion->height-1)&&j<(ImgRegion->width-1))
			{   //四领域；
				top=pData[stepsPL*(i-1)+j];
				bottom=pData[stepsPL*(i+1)+j];
				left=pData[stepsPL*i+(j-1)];
				right=pData[stepsPL*i+(j+1)];
				if(top==178||bottom==178||left==178||right==178)
				{
					pData[stepsPL*i+j]=255;
				}
			}
		}
}

//找轮廓，并将轮廓灰度值值为255；
void FindContour(IplImage * roi_region)
{
	int i,j;
	int stepsPL =roi_region->widthStep;
	BYTE * pData = (BYTE*)roi_region->imageData;
	int lheight=roi_region->height;
	int lwidth=roi_region->width;
	CvMemStorage* storage = cvCreateMemStorage(0);	
	CvSeq* contours = NULL;
	for(i=0;i<lheight;i++)
		for(j=0;j<lwidth;j++)
		{
			if(pData[i*stepsPL+j]==255)
				pData[i*stepsPL+j]=178;
		}
    cvFindContours( roi_region, storage, &contours, sizeof(CvContour),    
                    CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, cvPoint(0,0) );
	if (contours != NULL)
	{
		CvPoint* pointArray = (CvPoint*)malloc(sizeof(CvPoint) * contours->total);
		cvCvtSeqToArray(contours, pointArray, CV_WHOLE_SEQ);
		int total = contours->total;
	    for(i=0; i<total; i++)
		{
		    int x = pointArray[i].x;
		    int y = pointArray[i].y;
		    pData[y*stepsPL+x] = 255;
		}
		delete pointArray;
	}
    cvReleaseMemStorage(&storage);

}
//区域生长
long RegionGrowth(CPoint &seed, bool **visited ,IplImage * roi_region)
{
	//循环变量
	LONG i,j;
	USHORT lpSrc;
	int stepsPL =roi_region->widthStep;
	BYTE * pData = (BYTE*)roi_region->imageData;
	visited[seed.x][seed.y]=true;   //入队时置1
	CList <CPoint, CPoint&> queue(10) ;
	//种子入队  
	queue.AddTail(seed);    	
	long countInRegion=1;
	while(!queue.IsEmpty())  //队列非空
	{	
		CPoint curPoint;
		//取队列头元素
		curPoint=(CPoint)queue.RemoveHead(); 
		//扩展种子点的8领域
		for(i=-1;i<=1;i++)
		{
		    for(j=-1;j<=1;j++)
			{
				if(curPoint.x+i<0  || curPoint.x+i>=125 ||
				   curPoint.y+j<0  || curPoint.y+j>=125 )   //越界处理
				   continue;
                if(i==0 && j==0 )
				   continue;
				lpSrc=pData[(curPoint.x+i)*stepsPL+(curPoint.y+j)];
			    if(!visited[curPoint.x+i][curPoint.y+j])
				{
				    visited[curPoint.x+i][curPoint.y+j]=true;
				    if(lpSrc>0) 
					{
						//入队
						CPoint newSeed(curPoint.x+i,curPoint.y+j);
						queue.AddTail(newSeed);
						countInRegion++;
					}

				}
			} 
		}  
	} //队列为空
	return countInRegion;
}
//找最大的连通区域作为肿块区域
void FindMostLargeRegion(IplImage * roi_region)
{
	 CPoint resultCenter(0,0);   //最大连通区域内部的一点；
	 LONG maxCountInRegion=0;     //最大连通区域的点数；
	 LONG tempCountInRegion;      //临时变量；
	 int stepsPL =roi_region->widthStep;
	 BYTE * pData = (BYTE*)roi_region->imageData;
	 int lHeight=roi_region->height;
	 int lWidth=roi_region->width;
     int i,j;
	 USHORT lpSrc;
	 bool** visited=new bool*[lHeight];   //访问标志
	 for(i=0;i<lHeight;i++)
		visited[i]=new bool[lWidth]; 
	 for(i=0;i<lHeight;i++)
		for(j=0;j<lWidth;j++)
		{
			visited[i][j]=false;
		}
	 for(i=0;i<lHeight;i++)
		for(j=0;j<lWidth;j++)
		{
			lpSrc=pData[i*stepsPL+j];
			if(lpSrc>0&&(visited[i][j]==false))  //Not the background and not visited
			{
				CPoint startPoint(i,j);
				tempCountInRegion=RegionGrowth(startPoint,visited,roi_region);     //区域增长；
				if (tempCountInRegion>maxCountInRegion)
				{
					maxCountInRegion=tempCountInRegion;
					resultCenter.x=i;
					resultCenter.y=j;
				} //if
			}   //if
			else
				visited[i][j]=true;
		}  //for
	 for(i=0;i<lHeight;i++)
		for(j=0;j<lWidth;j++)
		{
			visited[i][j]=false;
		}  //for
    //对该最大连通区域进行再次区域生长；以将背景区域都抹黑；
	RegionGrowth(resultCenter,visited,roi_region);
	for(i=0;i<lHeight;i++)
		for(j=0;j<lWidth;j++)
		{
			if(!visited[i][j])
			{
				pData[i*stepsPL+j]=0;
			}
		}
	for(i=0;i<lHeight;i++)
	{
		delete[] visited[i];
	}
	delete visited;
}

/*******输入：一个CROI对象,输出: bool值，如果初始肿块经过snake分割之后成为一条曲线，
        那么返回false,否则返回true;如果初始轮廓太大，那么重新取初始轮廓；********/
bool AddSnakeSegmentWithRound(CROI& roi)
{
	USHORT * roi_region;    //定义原始roi图像，在这里为背景趋势去除后的图像；
	USHORT * img_buf;
	USHORT * Contour_buf;    //记录肿块轮廓的缓冲区；
	int Contour_count;       //轮廓点数；
	USHORT *growth_buf;      //记录肿块区域的缓冲区；
	int i,j,k;
	int count=0;       //snake后轮廓内部点数。
	int *x_array;      //轮廓点的横坐标（列向） ；
	int *y_array;      //轮廓点的纵坐标（横向） ；
      //申请空间；
	x_array=ivector(1,125*125/3);
	y_array=ivector(1,125*125/3);
	Contour_buf=svector(125,125);
	growth_buf=svector(125,125);
	roi_region=svector(125,125);
	img_buf=svector(125,125);
	IplImage* enhanceImg = cvCreateImage(cvSize(125, 125), 8, 1);
	IplImage* temp = cvCreateImage(cvSize(125, 125), 8, 1);
	//背景趋势去除		
	ImageProcess::Convert32Fto8U(roi.m_roiImg,temp);
	BackgroundTrendRemove(temp,enhanceImg);

	int stepsPLRoi = 8 * enhanceImg->widthStep / enhanceImg->depth;
	BYTE * pRoiData = (BYTE*)enhanceImg->imageData;

	int stepsPLRegion= roi.m_region->widthStep;
	BYTE * pRegionData=(BYTE*)roi.m_region->imageData;
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			*(img_buf+125*i+j)=(USHORT)(pRoiData[i*stepsPLRoi+j]*16);  //将enhanceImg拷贝到img_buf中；
			*(Contour_buf+125*i+j)=0;
		}
    if(IsGetRoundContour(roi.m_region))  //判断是否重新取初始轮廓；
	{
	    USHORT  * Para_a, * Para_b, * Para_c, * OrigImg,*EdgeMap_buf;
		int rows=125;
		int cols=125;
	    Para_a = svector(rows, cols);
	    Para_b = svector(rows, cols);
	    Para_c = svector(rows, cols);
	    OrigImg = svector(rows, cols);
		EdgeMap_buf= svector(rows, cols);
		int radius;
		copy_svector(img_buf, OrigImg, rows, cols);	
	    MinQuare(Para_a,Para_b,Para_c,EdgeMap_buf,OrigImg);  //最小二乘法就梯度图；
        radius=FindMassRadius(EdgeMap_buf);     //求新的圆形初始轮廓的半径；
		GetContourPoint(radius,Contour_buf,x_array,y_array);    //获得新的轮廓点数；
		//释放空间；
		free_svector(Para_a);
		free_svector(Para_b);
		free_svector(Para_c);
		free_svector(OrigImg);
		free_svector(EdgeMap_buf);
	}
	else
	{
		int index;
        index=0;
	    for(i=0;i<125;i++)
		   for(j=0;j<125;j++)
		   {
			   if(pRegionData[i*stepsPLRegion+j]==255)    //轮廓点的灰度值为最大值255；
			   {
				   *(Contour_buf+125*i+j)=4095;    
				   x_array[index]=i;             //初始化轮廓点数组；
				   y_array[index]=j;
				   index++;
			   }
			   else
			   {
				   *(Contour_buf+125*i+j)=0;     //非轮廓点的灰度值令其为0；
			   }
		   }
	}

    Contour_count=record_contour_array(Contour_buf,125,125,x_array,y_array); //将曲线上的点按行走路线有序存储；
	ActiveContourRegion(growth_buf,125,125,Contour_count,x_array,y_array);   //求包含肿块区域的矩形区域，并记录在缓冲区growth_buf中；
	Contour_count=Computing_ActiveContour(img_buf,growth_buf,125,125,Contour_count,x_array,y_array);  //计算活动轮廓；
	//令轮廓灰度值为255，其他点的灰度值为178；
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			pRegionData[i*stepsPLRegion+j]=178;
		}
	for(k=0;k<Contour_count;k++)
	{
		i=x_array[k];
		j=y_array[k];
		pRegionData[j*stepsPLRegion+i]=255;
	}
	//根据边界点修改区域标志：0---背景，178---区域内像素，255---边界像素					
	Cal_RegionFlag(roi.m_region);
	//计算肿块内部点数，若为0，就释放空间，返回false；
	count=0;
    for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			if(pRegionData[i*stepsPLRegion+j]==178)
				count++;
		}
	if(count==0)
	{
		free_ivector(x_array,1,125*125/3);
	    free_ivector(y_array,1,125*125/3);
	    free_svector(Contour_buf);
	    free_svector(growth_buf);
	    free_svector(roi_region);
	    free_svector(img_buf);
	    cvReleaseImage(&enhanceImg);
	    cvReleaseImage(&temp);
		return false;
	}
/**************************后处理***************************/
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			if(pRegionData[i*stepsPLRegion+j]==255)
				pRegionData[i*stepsPLRegion+j]=0;
		}
	GetRegionContour(roi.m_region);  //基于四邻域取肿块区域轮廓；
	
//找最大的连通区域作为肿块区域。

	FindMostLargeRegion(roi.m_region);
	FindContour(roi.m_region);   //找轮廓，并将轮廓灰度值值为255；

/**************************释放空间***************************/						
	free_ivector(x_array,1,125*125/3);
	free_ivector(y_array,1,125*125/3);
	free_svector(Contour_buf);
	free_svector(growth_buf);
	free_svector(roi_region);
	free_svector(img_buf);
	cvReleaseImage(&enhanceImg);
	cvReleaseImage(&temp);
	return true;
}
/*******输入：一个CROI对象,输出: bool值，如果初始肿块经过snake分割之后成为一条曲线，
        那么返回false,否则返回true。********/
bool AddSnakeSegment(CROI& roi)
{
	USHORT * roi_region;    //定义原始roi图像，在这里为背景趋势去除后的图像；
	USHORT * img_buf;
	USHORT * Contour_buf;    //记录肿块轮廓的缓冲区；
	int Contour_count;       //轮廓点数；
	USHORT *growth_buf;      //记录肿块区域的缓冲区；
	int i,j,k;
	int count=0;       //snake后轮廓内部点数。
	int *x_array;      //轮廓点的横坐标（列向） ；
	int *y_array;      //轮廓点的纵坐标（横向） ；
	//申请空间；
	x_array=ivector(1,125*125/3);
	y_array=ivector(1,125*125/3);
	Contour_buf=svector(125,125);
	growth_buf=svector(125,125);
	roi_region=svector(125,125);
	img_buf=svector(125,125);
	IplImage* enhanceImg = cvCreateImage(cvSize(125, 125), 8, 1);
	IplImage* temp = cvCreateImage(cvSize(125, 125), 8, 1);
	//背景趋势去除		
	ImageProcess::Convert32Fto8U(roi.m_roiImg,temp);  //将32位浮点数转换为8位的；
/*
	cvNamedWindow("jyx",1);
    cvShowImage("jyx",temp);
    AfxMessageBox("nihao",MB_OK);*/

	BackgroundTrendRemove(temp,enhanceImg); 
	

		
	int stepsPLRoi = 8 * enhanceImg->widthStep / enhanceImg->depth;
	BYTE * pRoiData = (BYTE*)enhanceImg->imageData;

	int stepsPLRegion= roi.m_region->widthStep;
	BYTE * pRegionData=(BYTE*)roi.m_region->imageData;
	int index;
    index=0;
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			*(img_buf+125*i+j)=(USHORT)(pRoiData[i*stepsPLRoi+j]*16);   //将enhanceImg拷贝到img_buf中；
			if(pRegionData[i*stepsPLRegion+j]==255)   //轮廓点的灰度值为最大值255；
			{
				*(Contour_buf+125*i+j)=4095;
				x_array[index]=i;
				y_array[index]=j;        //初始化轮廓点数组；
				index++;
			}
			else
			{
				*(Contour_buf+125*i+j)=0;    //非轮廓点的灰度值令其为0；        
			}
		}
    Contour_count=record_contour_array(Contour_buf,125,125,x_array,y_array);   //将曲线上的点按行走路线有序存储；
	ActiveContourRegion(growth_buf,125,125,Contour_count,x_array,y_array);     //求包含肿块区域的矩形区域，并记录在缓冲区growth_buf中；
	Contour_count=Computing_ActiveContour(img_buf,growth_buf,125,125,Contour_count,x_array,y_array);       //计算活动轮廓；
	//令轮廓灰度值为255，其他点的灰度值为178；
	for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			pRegionData[i*stepsPLRegion+j]=178;
		}
	for(k=0;k<Contour_count;k++)
	{
		i=x_array[k];
		j=y_array[k];
		pRegionData[j*stepsPLRegion+i]=255;
	}
	//根据边界点修改区域标志：0---背景，178---区域内像素，255---边界像素				
	Cal_RegionFlag(roi.m_region);  
	//计算肿块内部点数，若为0，就释放空间，返回false；
	count=0;
    for(i=0;i<125;i++)
		for(j=0;j<125;j++)
		{
			if(pRegionData[i*stepsPLRegion+j]==178)
				count++;
		}
	if(count==0)
	{
		free_ivector(x_array,1,125*125/3);
	    free_ivector(y_array,1,125*125/3);
	    free_svector(Contour_buf);
	    free_svector(growth_buf);
	    free_svector(roi_region);
	    free_svector(img_buf);
	    cvReleaseImage(&enhanceImg);
	    cvReleaseImage(&temp);
		return false;
	}
  //测试；
/*
  cvNamedWindow("jyx",CV_WINDOW_AUTOSIZE);
    cvShowImage("jyx",roi.m_region);
    AfxMessageBox("nihao!");*/
  
/**************************后处理***************************/
for(i=0;i<125;i++)
	for(j=0;j<125;j++)
	{
		if(pRegionData[i*stepsPLRegion+j]==255)
			pRegionData[i*stepsPLRegion+j]=0;
	}
GetRegionContour(roi.m_region);   //基于四邻域取肿块区域轮廓；
				
//找最大的连通区域作为肿块区域。
			
FindMostLargeRegion(roi.m_region);			
FindContour(roi.m_region);    //找轮廓，并将轮廓灰度值值为255；

/**************************释放空间***************************/						
	free_ivector(x_array,1,125*125/3);
	free_ivector(y_array,1,125*125/3);
	free_svector(Contour_buf);
	free_svector(growth_buf);
	free_svector(roi_region);
	free_svector(img_buf);
	cvReleaseImage(&enhanceImg);
	cvReleaseImage(&temp);
	return true;
}