//
//    RRT - A numerical kernel for Relativistic Ray Tracing.
//    Copyright (C) 2011  Renaud Sirdey (renaud.sirdey@gmail.com)
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#define SCREEN_Z 100.0
#define SCREEN_BASE 1.0 //0.0625 //0.125//0.5
#define SCREEN_LEFT_X -1.0 //-SCREEN_BASE
#define SCREEN_TOP_Y 1.0 //0.75 //1.0 //SCREEN_BASE
#define SCREEN_RIGHT_X 1.0 //SCREEN_BASE
#define SCREEN_BOTTOM_Y -1.0 //-0.75 //-1.0 //-SCREEN_BASE

#define OBS_DELTA_Z 1.0

#define BLACK_HOLE_MASS 4 //6 //8 //4 //16 //4//2//0.5//0.01

#define WALL_DIST 100.0
#define WALL_PAPER_LENGTH 10.0

#define IMAG_WIDTH 512//2000//2000//1024//3648 //256 //7000
#define IMAG_HEIGHT 512//2000//2000//1024//3648 //256 //7000

#define PREC 0.001

#define MIN_DIST (WALL_DIST-10.0*PREC)

#define SQR(x) ((x)*(x))

#define BLACK 0//0
#define RED 1
#define ORANGE 2
#define YELLOW 3
#define GREEN 4
#define BLUE 5
#define PURPLE 6
#define WHITE 7

#define MAX(a,b) ((a)>(b)?(a):(b))

unsigned char imag[IMAG_WIDTH*IMAG_HEIGHT];

#ifdef __TEXTURE__

/*#define TEXTURE_FILE "Scientifique-relativite-carre.bmp"
#define TEXTURE_WIDTH 1000
#define TEXTURE_HEIGHT 1000*/

/*#define TEXTURE_FILE "Cray_200x200.bmp"
#define TEXTURE_WIDTH 200
#define TEXTURE_HEIGHT 200*/

// Escher/
#define TEXTURE_FILE "LW306_512x512.bmp"
#define TEXTURE_WIDTH 512
#define TEXTURE_HEIGHT 512

/*#define TEXTURE_FILE "Flash_dragon_614x614.bmp"
#define TEXTURE_WIDTH 614
#define TEXTURE_HEIGHT 614
*/
/*#define TEXTURE_FILE "Buzz-4.bmp"
#define TEXTURE_WIDTH 1600
#define TEXTURE_HEIGHT 1600*/
/*#define TEXTURE_FILE "Buzz-3.bmp"
#define TEXTURE_WIDTH 512
#define TEXTURE_HEIGHT 512*/
/*#define TEXTURE_FILE "Buzz-2.bmp"
#define TEXTURE_WIDTH 512
#define TEXTURE_HEIGHT 512*/
/*#define TEXTURE_FILE "Buzz.bmp"
#define TEXTURE_WIDTH 1024
#define TEXTURE_HEIGHT 1024*/
/*#define TEXTURE_FILE "Hubble_deep_field.bmp"
#define TEXTURE_WIDTH 3100
#define TEXTURE_HEIGHT 3100*/
/*
#define TEXTURE_FILE "Flash.bmp"
#define TEXTURE_WIDTH 700
#define TEXTURE_HEIGHT 700 */
unsigned char texture[TEXTURE_WIDTH*TEXTURE_HEIGHT];
#endif

double schwarzschildPotential(const double r,const double M)
{
	assert(r>=2.0*M);
	return((1-2*M/r)/SQR(r));
}

double invSchwarzschildPotential(const double b,const double M)
{
	const double err=1e-6;
	double lhs,a1,a2;

	assert(b>=sqrt(27)*M);

	lhs=1/SQR(b);

    a1=3.0*M;
    a2=3.0*M;

	// Left widening.
	while(schwarzschildPotential(a1,M)<=lhs) a1=a1/2;

	// Right widening.
	while(schwarzschildPotential(a2,M)>=lhs) a2=a2*2;

	// Dichotomy.
	while(a2-a1>=err)
	{
		double a;
		a=(a1+a2)/2.0;
		if(schwarzschildPotential(a,M)>lhs)
			a1=a;
		else
			a2=a;
	}

	return((a1+a2)/2.0);
}

int checkIntersection(const int i,const int j,const double x,const double y,const double z)
{
	if(x>=WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(y+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(z+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		//imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];
		
		imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];//255;//texture[k*TEXTURE_WIDTH+l];
		
#else
#ifdef __CIRCLES__
		if(SQR(y)+SQR(z)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=RED;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(y+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(z+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l)%2==0)
			imag[i+j*IMAG_WIDTH]=RED;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}
	if(x<=-WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(y+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(z+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		//imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+(TEXTURE_WIDTH-l-1)];
		
		imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];//255;//texture[k*TEXTURE_WIDTH+l];
		
#else
#ifdef __CIRCLES__
		if(SQR(y)+SQR(z)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=ORANGE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(y+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(z+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l+1)%2==0)
			imag[i+j*IMAG_WIDTH]=ORANGE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}

	if(y>=WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(x+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(z+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		//imag[i+j*IMAG_WIDTH]=texture[k+l*TEXTURE_WIDTH];
		
		imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];//255;//texture[k*TEXTURE_WIDTH+l];
		
#else
#ifdef __CIRCLES__
		if(SQR(x)+SQR(z)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=YELLOW;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(x+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(z+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l)%2==0)
			imag[i+j*IMAG_WIDTH]=YELLOW;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}
	if(y<=-WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(x+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(z+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		//imag[i+j*IMAG_WIDTH]=texture[k+(TEXTURE_HEIGHT-l-1)*TEXTURE_WIDTH];
		
		imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];//255;//texture[k*TEXTURE_WIDTH+l];
		
#else
#ifdef __CIRCLES__
		if(SQR(x)+SQR(z)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=GREEN;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(x+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(z+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l+1)%2==0)
			imag[i+j*IMAG_WIDTH]=GREEN;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}

	if(z>=WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(x+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(y+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		//imag[i+j*IMAG_WIDTH]=texture[(TEXTURE_WIDTH-k-1)+l*TEXTURE_WIDTH];
		
		imag[i+j*IMAG_WIDTH]=texture[k*TEXTURE_WIDTH+l];//255;//texture[k*TEXTURE_WIDTH+l];
		
#else
#ifdef __CIRCLES__
		if(SQR(x)+SQR(y)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=BLUE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(x+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(y+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l+1)%2==0)
			imag[i+j*IMAG_WIDTH]=BLUE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}

	if(z<=-WALL_DIST)
	{
		int k,l;
#ifdef __TEXTURE__
		k=(int)MAX(0.0,(x+WALL_DIST)*(TEXTURE_WIDTH-1)/(2.0*WALL_DIST));
		l=(int)MAX(0.0,(y+WALL_DIST)*(TEXTURE_HEIGHT-1)/(2.0*WALL_DIST));
		assert(k>=0 && l>=0);
		assert(k<TEXTURE_WIDTH && l<TEXTURE_HEIGHT);
		imag[i+j*IMAG_WIDTH]=texture[k+l*TEXTURE_WIDTH];
#else
#ifdef __CIRCLES__
		if(SQR(x)+SQR(y)<SQR(WALL_DIST))
			imag[i+j*IMAG_WIDTH]=PURPLE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#else			
		k=(int)floor(MAX(0.0,(x+WALL_DIST)/WALL_PAPER_LENGTH));
		l=(int)floor(MAX(0.0,(y+WALL_DIST)/WALL_PAPER_LENGTH));
		assert(k>=0 && l>=0);
		if((k+l)%2==0)
			imag[i+j*IMAG_WIDTH]=PURPLE;
		else
			imag[i+j*IMAG_WIDTH]=WHITE;
#endif
#endif
		return(1);
	}

	return(0);
}

void processOneRay(const int i,const int j,const double M)
{
	double x,y,r0,theta0,b;

	x=SCREEN_LEFT_X+((double)i)*(SCREEN_RIGHT_X-SCREEN_LEFT_X)/IMAG_WIDTH;
	y=SCREEN_TOP_Y+((double)j)*(SCREEN_BOTTOM_Y-SCREEN_TOP_Y)/IMAG_HEIGHT;

	r0=SCREEN_Z+OBS_DELTA_Z;

	{
		double prod,norm,cos_theta0;
		prod=-OBS_DELTA_Z;
		norm=sqrt(SQR(x)+SQR(y)+SQR(-OBS_DELTA_Z));
		cos_theta0=prod/norm;
		theta0=acos(cos_theta0);
	}

	assert(theta0>=M_PI/2.0 && theta0<=M_PI);

	b=cos(theta0-M_PI/2)*r0;

	if(b<=M*sqrt(27)) /* Rayon absorbé */
	{
		imag[i+j*IMAG_WIDTH]=BLACK;
		imag[(IMAG_WIDTH-i-1)+j*IMAG_WIDTH]=BLACK;
		imag[i+(IMAG_HEIGHT-j-1)*IMAG_WIDTH]=BLACK;
		imag[(IMAG_WIDTH-i-1)+(IMAG_HEIGHT-j-1)*IMAG_WIDTH]=BLACK;
	}
	else /* Rayon dévié */
	{
		const double err=1e-4,dt=PREC;
		double d,phi,r,sr,cos_alpha,sin_alpha;

		d=invSchwarzschildPotential(b,M);

		cos_alpha=x/sqrt(pow(x,2.0)+pow(y,2.0));
		sin_alpha=y/sqrt(pow(x,2.0)+pow(y,2.0));

		phi=0.0;
		r=r0;
		sr=-1.0;

		while(1)
		{
			double phi1,r1,f,f2,g;
			f=1.0-2.0*M/r;
			f2=f*f;
			g=f2-f2*f*(b*b)/(r*r);
			//g=pow(f,2.0)-pow(f,3.0)*pow(b/r,2.0);
			assert(g>=0.0);
			r1=r + sr * sqrt(g) * dt;
			phi1=phi + dt*b*(f)/(r*r);
			r=r1;
			phi=phi1;

			if(r<d+err && sr==-1.0)
			{
				sr=1.0;
			}

			if(sr==1.0 && r>=MIN_DIST)
			{
				double x0,z0,x,y,z;

				z0=r*cos(phi);
				x0=r*sin(phi);

				x=cos_alpha*x0;
				y=-sin_alpha*x0;
				z=z0;

				if(checkIntersection(i,j,x,y,z))
				{
					int k=0;
					x=-cos_alpha*x0;
					y=-sin_alpha*x0;
					k+=checkIntersection(IMAG_WIDTH-i-1,j,x,y,z);
					x=cos_alpha*x0;
					y=sin_alpha*x0;
					k+=checkIntersection(i,IMAG_HEIGHT-j-1,x,y,z);
					x=-cos_alpha*x0;
					y=sin_alpha*x0;
					k+=checkIntersection(IMAG_WIDTH-i-1,IMAG_HEIGHT-j-1,x,y,z);
					assert(k==3);

					break;
				}
			}
		}
	}
}

typedef double v2df __attribute__ ((mode(V2DF)));

typedef union
{
	v2df v;
	double d[2];
}double_vec2;

void processTwoRays(const int i0,const int j0,const int i1,const int j1,const double M)
{
	const int i[2]={i0,i1},j[2]={j0,j1};
	double_vec2 x,y,theta0,b;
	double r0;
	int k;

	for(k=0;k<2;k++)
	{
		x.d[k]=SCREEN_LEFT_X+((double)i[k])*(SCREEN_RIGHT_X-SCREEN_LEFT_X)/IMAG_WIDTH;
		y.d[k]=SCREEN_TOP_Y+((double)j[k])*(SCREEN_BOTTOM_Y-SCREEN_TOP_Y)/IMAG_HEIGHT;

		r0=SCREEN_Z+OBS_DELTA_Z;

		{
			double prod,norm,cos_theta0;
			prod=-OBS_DELTA_Z;
			norm=sqrt(SQR(x.d[k])+SQR(y.d[k])+SQR(-OBS_DELTA_Z));
			cos_theta0=prod/norm;
			theta0.d[k]=acos(cos_theta0);
		}

		assert(theta0.d[k]>=M_PI/2.0 && theta0.d[k]<=M_PI);

		b.d[k]=cos(theta0.d[k]-M_PI/2)*r0;

		if(b.d[k]<=M*sqrt(27)) /* Rayon 0 absorbé */
		{
			imag[i[k]+j[k]*IMAG_WIDTH]=BLACK;
			imag[(IMAG_WIDTH-i[k]-1)+j[k]*IMAG_WIDTH]=BLACK;
			imag[i[k]+(IMAG_HEIGHT-j[k]-1)*IMAG_WIDTH]=BLACK;
			imag[(IMAG_WIDTH-i[k]-1)+(IMAG_HEIGHT-j[k]-1)*IMAG_WIDTH]=BLACK;
			/* Traitement de l'autre rayon individuellement. */
			processOneRay(i[(k+1)&1],j[(k+1)&1],M);
			return;
		}
	}

	/* Les deux rayons sont non absorbés. */
	{
		const double err=1e-4;
		const double_vec2 dt={.d={PREC,PREC}}, one={.d={1.0,1.0}},twoM={.d={2.0*M,2.0*M}};
		double_vec2 d,phi,r,sr,cos_alpha,sin_alpha;

		int inter[2]={0,0},n_inter=0; // Nombre de rayons à avoir intersectés. */

		for(k=0;k<2;k++)
		{
			imag[i[k]+j[k]*IMAG_WIDTH]=BLACK;

			d.d[k]=invSchwarzschildPotential(b.d[k],M);

			{
				double d;
				d=sqrt(x.d[k]*x.d[k]+y.d[k]*y.d[k]);
				cos_alpha.d[k]=x.d[k]/d;
				sin_alpha.d[k]=y.d[k]/d;
			}

			phi.d[k]=0.0;
			r.d[k]=r0;
			sr.d[k]=-1.0;
		}

		while(n_inter<2)
		{
			double_vec2 phi1,r1,f,f2,g,sqrtg;

#ifndef __NO_GCC_VEC_EXT__
			double_vec2 r2;
			r2.v=r.v*r.v;
			f.v=one.v-twoM.v/r.v;
			f2.v=f.v*f.v;
			g.v=f2.v-f2.v*f.v*(b.v*b.v)/r2.v;
			sqrtg.v=__builtin_ia32_sqrtpd(g.v);
			assert(g.d[0]>=0.0 && g.d[1]>=0.0);
			r1.v=r.v + sr.v * sqrtg.v * dt.v;
			phi1.v=phi.v + dt.v*b.v*f.v/r2.v;
			r=r1;
			phi=phi1;
#else
			f.d[0]=one.d[0]-twoM.d[1]/r.d[0];
			f.d[1]=one.d[1]-twoM.d[1]/r.d[1];
			f2.d[0]=f.d[0]*f.d[0];
			f2.d[1]=f.d[1]*f.d[1];
			g.d[0]=f2.d[0]-f2.d[0]*f.d[0]*(b.d[0]*b.d[0])/(r.d[0]*r.d[0]);
			g.d[1]=f2.d[1]-f2.d[1]*f.d[1]*(b.d[1]*b.d[1])/(r.d[1]*r.d[1]);
			assert(g.d[0]>=0.0 && g.d[1]>=0.0);
			r1.d[0]=r.d[0] + sr.d[0] * sqrt(g.d[0]) * dt.d[0];
			r1.d[1]=r.d[1] + sr.d[1] * sqrt(g.d[1]) * dt.d[1];
			phi1.d[0]=phi.d[0] + dt.d[0]*b.d[0]*f.d[0]/(r.d[0]*r.d[0]);
			phi1.d[1]=phi.d[1] + dt.d[1]*b.d[1]*f.d[1]/(r.d[1]*r.d[1]);
			r.d[0]=r1.d[0];
			r.d[1]=r1.d[1];
			phi.d[0]=phi1.d[0];
			phi.d[1]=phi1.d[1];
#endif

			for(k=0;k<2;k++)
			{
				if(r.d[k]<d.d[k]+err && sr.d[k]==-1.0)
				{
					sr.d[k]=1.0;
				}

				if(inter[k]==0 && sr.d[k]==1.0 && r.d[k]>=MIN_DIST)
				{
					double x0,z0,x,y,z;

					z0=r.d[k]*cos(phi.d[k]);
					x0=r.d[k]*sin(phi.d[k]);

					x=cos_alpha.d[k]*x0;
					y=-sin_alpha.d[k]*x0;
					z=z0;

					if(checkIntersection(i[k],j[k],x,y,z))
					{
						int ctr=0;
						x=-cos_alpha.d[k]*x0;
						y=-sin_alpha.d[k]*x0;
						ctr+=checkIntersection(IMAG_WIDTH-i[k]-1,j[k],x,y,z);
						x=cos_alpha.d[k]*x0;
						y=sin_alpha.d[k]*x0;
						ctr+=checkIntersection(i[k],IMAG_HEIGHT-j[k]-1,x,y,z);
						x=-cos_alpha.d[k]*x0;
						y=sin_alpha.d[k]*x0;
						ctr+=checkIntersection(IMAG_WIDTH-i[k]-1,IMAG_HEIGHT-j[k]-1,x,y,z);
						assert(ctr==3);

						inter[k]=1;
						n_inter++;
					}
				}
			}
		}
	}
}

#include <pthread.h>

#define NUM_OF_THREADS 4

void *doWork(const int * const firstLine)
{
	int i,j,k=0;

	assert(*firstLine>=0);

	for(j=*firstLine;j<IMAG_HEIGHT/2;j+=NUM_OF_THREADS)
#ifdef __NO_SIMD__
		for(i=0;i<IMAG_WIDTH/2;i++)
		{
			processOneRay(i,j,BLACK_HOLE_MASS);
			if(*firstLine==0)
			{
				printf("[%3d%%]\r",(int)(100.0*4.0*NUM_OF_THREADS*((double)k)/((double)(IMAG_WIDTH*IMAG_HEIGHT))));
				fflush(stdout);
			}
			k++;
		}
#else
		for(i=0;i<IMAG_WIDTH/2;i+=2)
		{
			processTwoRays(i,j,i+1,j,BLACK_HOLE_MASS);
			if(*firstLine==0)
			{
				printf("[%3d%%]\r",(int)(100.0*4.0*NUM_OF_THREADS*((double)k)/((double)(IMAG_WIDTH*IMAG_HEIGHT))));
				fflush(stdout);
			}
			k+=2;
		}
#endif
	if(*firstLine==0)
		printf("[100%%]\n");

	return(NULL);
}

#ifdef __TEXTURE__
void readBmp(const char * const,unsigned char * const);
#endif
void writeBmp(const char * const,const unsigned char * const);

int main(void)
{

	int k;
	pthread_t tid[NUM_OF_THREADS];
	int firstLines[NUM_OF_THREADS];

	printf("Resolution : %dx%d\n",IMAG_WIDTH,IMAG_HEIGHT);

#ifdef __TEXTURE__
	printf("Loading %s\n",TEXTURE_FILE);
	readBmp(TEXTURE_FILE,texture);
#endif

#ifdef __NO_THREADS__
	assert(NUM_OF_THREADS==1);
	firstLines[0]=0;
	doWork(&firstLines[0]);
#else
	printf("Spawning %d threads\n",NUM_OF_THREADS);

	for(k=0;k<NUM_OF_THREADS;k++)
	{
		firstLines[k]=k;
		if(pthread_create(&tid[k],NULL,(void*(*)(void*))doWork,&firstLines[k])!=0)
		{
			printf("Thread creation error\n");
			exit(1);
		}
	}

	for(k=0;k<NUM_OF_THREADS;k++)
	{
		if(pthread_join(tid[k],NULL)!=0)
			printf("Suspicious termination of thread %d\n",k);
		else
			printf("Thread %d successfully terminated\n",k);
	}
#endif

	writeBmp("out-sse.bmp",imag);
}



/* Unimportant : helpers for reading/writing BMP files. */

typedef struct
{
	unsigned short bfType;
	unsigned int bfSize;
	unsigned short bfReserved1;
	unsigned short bfReserved2;
	unsigned int bfOffBits;
}BmpFileHeader;

typedef struct
{
	unsigned int biSize;
	unsigned int biWidth;
	unsigned int biHeight;
	unsigned short biPlanes;
	unsigned short biBitCount;
	unsigned int biCompression;
	unsigned int biSizeImage;
	unsigned int biXPelsPerMeter;
	unsigned int biYPelsPerMeter;
	unsigned int biClrUsed;
	unsigned int biClrImportant;
}BmpInfoHeader;

typedef struct
{
	unsigned char rgbBlue;
	unsigned char rgbGreen;
	unsigned char rgbRed;
	unsigned char rgbReserved;
}BmpCmapEntry;

#ifdef __TEXTURE__

BmpCmapEntry cmap[256];

void readBmp(const char * const file,unsigned char * const data)
{
	BmpFileHeader fileHeader;
	BmpInfoHeader infoHeader;
	FILE *desc=fopen(file,"rb");
	assert(desc!=NULL);

	assert(sizeof(unsigned int)==4);
	assert(sizeof(unsigned short)==2);

	/* Reading one by one to avoid padding issues. */
	fread(&fileHeader.bfType,sizeof(unsigned short),1,desc);
	fread(&fileHeader.bfSize,sizeof(unsigned int),1,desc);
	fread(&fileHeader.bfReserved1,sizeof(unsigned short),1,desc);
	fread(&fileHeader.bfReserved2,sizeof(unsigned short),1,desc);
	fread(&fileHeader.bfOffBits,sizeof(unsigned int),1,desc);

	fread(&infoHeader,sizeof(BmpInfoHeader),1,desc);

	assert(TEXTURE_WIDTH==infoHeader.biWidth && TEXTURE_HEIGHT==infoHeader.biHeight);
	assert(infoHeader.biBitCount==8);

	// Les mystères du format BMP...
	//fseek(desc,20*4,SEEK_CUR);

	fread(cmap,sizeof(BmpCmapEntry),256,desc);

	fseek(desc,fileHeader.bfOffBits,SEEK_SET);

	fread(data,sizeof(unsigned char),infoHeader.biWidth*infoHeader.biHeight,desc);

	fclose(desc);
}
#endif

void writeBmp(const char * const file,const unsigned char * const data)
{
	BmpFileHeader fileHeader;
	BmpInfoHeader infoHeader;
#ifndef __TEXTURE__
	BmpCmapEntry cmap[256];
#endif
	FILE *desc=fopen(file,"wb");
	assert(desc!=NULL);

	fileHeader.bfType=19778;
	fileHeader.bfSize=1078+IMAG_WIDTH*IMAG_HEIGHT;
	fileHeader.bfReserved1=0;
	fileHeader.bfReserved2=0;
	fileHeader.bfOffBits=1078;

	infoHeader.biSize=40;
	infoHeader.biWidth=IMAG_WIDTH;
	infoHeader.biHeight=IMAG_HEIGHT;
	infoHeader.biPlanes=1;
	infoHeader.biBitCount=8;
	infoHeader.biCompression=0;
	infoHeader.biSizeImage=0;
	infoHeader.biXPelsPerMeter=0;
	infoHeader.biYPelsPerMeter=0;
	infoHeader.biClrUsed=0;
	infoHeader.biClrImportant=0;

#ifndef __TEXTURE__
	{
		int i;
		for(i=0;i<256;i++)
		{
			cmap[i].rgbBlue=(unsigned char)i;
			cmap[i].rgbGreen=(unsigned char)i;
			cmap[i].rgbRed=(unsigned char)i;
			cmap[i].rgbReserved=0;
		}
		cmap[BLACK].rgbRed=   0; cmap[BLACK].rgbGreen=   0; cmap[BLACK].rgbBlue=   0; //
		cmap[RED].rgbRed=   255; cmap[RED].rgbGreen=     0; cmap[RED].rgbBlue=     0; //
		cmap[ORANGE].rgbRed=255; cmap[ORANGE].rgbGreen=165; cmap[ORANGE].rgbBlue=  0; //
		cmap[YELLOW].rgbRed=255; cmap[YELLOW].rgbGreen=255; cmap[YELLOW].rgbBlue=  0; //
		cmap[GREEN].rgbRed=   0; cmap[GREEN].rgbGreen= 255; cmap[GREEN].rgbBlue=   0; //
		cmap[BLUE].rgbRed=    0; cmap[BLUE].rgbGreen=    0; cmap[BLUE].rgbBlue=  255; //
		cmap[PURPLE].rgbRed=102; cmap[PURPLE].rgbGreen=  0; cmap[PURPLE].rgbBlue=153; //
		cmap[WHITE].rgbRed= 255; cmap[WHITE].rgbGreen= 255; cmap[WHITE].rgbBlue= 255; //
	}
#endif

	assert(sizeof(unsigned int)==4);
	assert(sizeof(unsigned short)==2);

	/* Writing one by one to avoid padding issues. */
	fwrite(&fileHeader.bfType,sizeof(unsigned short),1,desc);
	fwrite(&fileHeader.bfSize,sizeof(unsigned int),1,desc);
	fwrite(&fileHeader.bfReserved1,sizeof(unsigned short),1,desc);
	fwrite(&fileHeader.bfReserved2,sizeof(unsigned short),1,desc);
	fwrite(&fileHeader.bfOffBits,sizeof(unsigned int),1,desc);

	fwrite(&infoHeader,sizeof(BmpInfoHeader),1,desc);

	fwrite(cmap,sizeof(BmpCmapEntry),256,desc);

	fwrite(data,sizeof(unsigned char),infoHeader.biWidth*infoHeader.biHeight,desc);

	fclose(desc);

}



/*

dim. 16 nov. 2014 15:56:19
Resolution : 2000x2000
Loading Scientifique-relativite-carre.bmp
Spawning 2 threads
[100%]
Thread 0 successfully terminated
Thread 1 successfully terminated

real    88m20.855s
user    174m8.307s
sys     0m1.466s
dim. 16 nov. 2014 17:24:40

*/
