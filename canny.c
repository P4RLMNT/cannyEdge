/*********************
*  Author: Andrew Gates
*  Date: 10/12/2017
*  git repo:
*  For: CAP 4453 - Robot Vision Assignment One
*  
*/

/* DEVELOPMENT NOTES
do 8 cases for direction of gradient to trace possible edges to find real edge to assert possible edge in question.

change thesh array to have values of -1,0,1. and if an Edge gets asserted give it a 1 value;

also see if you can automatically come up with threshold values.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define  MAXMASK 100
#define  PICSIZE 256
#define  SIGONE 1
#define  SIGTWO 2

int getThreshold(double pixel);
void hysteresis(void);
void imageToFile(char*, int[256][256]);
int inBounds(int i, int j);
void magToFile(char *file, double output[256][256]);
void nonMaxSup(void);
void normalize(double maxVal);
void tracer(int i, int j);


double ival[PICSIZE][PICSIZE];
int final[PICSIZE][PICSIZE];
double gradDir[PICSIZE][PICSIZE];
int    pic[PICSIZE][PICSIZE];
int outpicx[PICSIZE][PICSIZE];
int outpicy[PICSIZE][PICSIZE];
int    peaks[PICSIZE][PICSIZE]; // non-maximum suppression
double xmask[MAXMASK][MAXMASK], ymask[MAXMASK][MAXMASK];
int thresh[PICSIZE][PICSIZE];
double hThresh = 110.0, lThresh = 30.0;
int visited[PICSIZE][PICSIZE];


int main(int argc, char** argv){
	int     i,j,p,q,mr,centx,centy;
	double  xmaskval,ymaskval,Gx,Gy,sig,maxival;
	FILE    *fo1,*fp1;
	char    *foobar, *pgmType;
	int picLength, picWidth, brightness;

	argc--; argv++;
	foobar = *argv;
	fp1=fopen(foobar,"rb");

	argc--; argv++;
	foobar = *argv;
	sig = atof(foobar);

	mr = (int)(sig * 3);
	centx = (MAXMASK / 2);
	centy = (MAXMASK / 2);

	fscanf(fp1, "%s", pgmType);
	fscanf(fp1, "%d %d", &picLength, &picWidth);
	fscanf(fp1, "%d", &brightness);
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
					{
						pic[i][j]  =  getc (fp1);
						final[i][j] = 0; // initializing final image to ZERO
						thresh[i][j] = 0;
						visited[i][j] = 0;
					}
	}

// convolution to create the x and y mask
	for (p=-mr;p<=mr;p++)
	{  for (q=-mr;q<=mr;q++)
		 {
				xmaskval = (q * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				ymaskval = (p * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				xmask[p+centy][q+centx] = xmaskval;
				ymask[p+centy][q+centx] = ymaskval;
		 }
	}
	for (i=mr;i<=255-mr;i++)
	{ for (j=mr;j<=255-mr;j++)
		{
			 Gx = 0.0;
			 Gy = 0.0;
			 for (p=-mr;p<=mr;p++)
			 {
					for (q=-mr;q<=mr;q++)
					{
						 Gx += pic[i+p][j+q] * xmask[p+centy][q+centx];
						 Gy += pic[i+p][j+q] * ymask[p+centy][q+centx];
					}
			 }
			 outpicx[i][j] = (int)Gx;
			 outpicy[i][j] = (int)Gy;
			 gradDir[i][j] = (atan2(Gx,Gy)/M_PI) * 180.0; // calculating gradient direction
		}
	}
// sobel code
maxival = 0.0; // value for scaling
// getting magnitude for each pixel
for (i=mr;i<256-mr;i++)
{ for (j=mr;j<256-mr;j++)
	{
		 ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
															(outpicy[i][j]*outpicy[i][j])));
		 if (ival[i][j] > maxival)
				maxival = ival[i][j];

	 }
}

normalize(maxival);
nonMaxSup();
hysteresis();

magToFile("cannyMag.pgm", ival);
imageToFile("cannyXGrad.pgm", outpicx);
imageToFile("cannyYGrad.pgm", outpicy);
imageToFile("cannyPeaks.pgm", peaks);
imageToFile("cannyThres.pgm", thresh);
imageToFile("cannyFinal.pgm", final);

 return 0;
} // end main

void hysteresis()
{
	int i, j;
	for(i = 1; i < PICSIZE-1; i++)
	{
		for(j = 1; j < PICSIZE-1; j++)
		{
			if(peaks[i][j] == 0)
			{
				final[i][j] = 0;
				visited[i][j] = 1;
			}
			else if(ival[i][j] > hThresh)
				{
					final[i][j] = 1;
					visited[i][j] = 1;
					// do trace on 8 neighbors
					tracer(i, j-1); // W
					tracer(i, j+1); // E
					tracer(i-1,j); // N
					tracer(i+1, j); // S
					tracer(i-1, j-1); // NW
					tracer(i-1, j+1); // NE
					tracer(i+1, j-1); // SW
					tracer(i+1, j+1); // SE
				}
		} // end j for
	} // end i for
} // end hysteresis

void tracer(int i, int j){
	if(inBounds(i,j) < 0 || visited[i][j] == 1) return;
	if(peaks[i][j] == 0)
	{
		final[i][j] = 0;
		visited[i][j] = 1;
		return;
	}
	else if(ival[i][j] > lThresh)
		{
			final[i][j] = 1;
			visited[i][j] = 1;
			tracer(i, j-1); // W
			tracer(i, j+1); // E
			tracer(i-1,j); // N
			tracer(i+1, j); // S
			tracer(i-1, j-1); // NW
			tracer(i-1, j+1); // NE
			tracer(i+1, j-1); // SW
			tracer(i+1, j+1); // SE
		}
} // end tracer

int getThreshold(double pixel)
{
	if(pixel > hThresh) return 1;
	else if(pixel > lThresh) return 0;
	else return -1;
} // end getThreshold

// peaks thresh gradDir
void nonMaxSup()
{
	int i, j;
	double theta;
	for(i = 0; i < PICSIZE; i++) // zeroing edges of peaks image
	{
		peaks[i][0] = 0; // left edge
		peaks[i][PICSIZE-1] = 0; // right edge

		thresh[i][0] = 0; // left edge
		thresh[i][PICSIZE-1] = 0; // right edge
	}
	for(j = 0; j < PICSIZE; j++)
	{
		peaks[0][j] = 0; // top edge
		peaks[PICSIZE-1][j] = 0; // bottom edge

		thresh[0][j] = 0; // top edge
		thresh[PICSIZE-1][j] = 0; // bottom edge
	}
	for(i = 1; i < PICSIZE-1; i++)
	{
		for(j = 1; j < PICSIZE-1; j++)
		{
			theta = gradDir[i][j];
			if(( (theta < 22.5) && (theta > -22.5) ) || (theta > 157.5) || (theta < -157.5)) // EW
			{
				if(ival[i-1][j] > ival[i][j] || ival[i+1][j] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
				gradDir[i][j] = 0.0;
			}
			if(( (theta > 22.5) && (theta < 67.5) ) || ( (theta < -112.5) && (theta > -157.5) )) // NE/SW
			{
				if(ival[i-1][j-1] > ival[i][j] || ival[i+1][j+1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
				gradDir[i][j] = 45.0;
			}
			if(( (theta > 67.5) && (theta < 112.5) ) || ( (theta < -67.5) && (theta > -112.5) )) // NS
			{
				if(ival[i][j-1] > ival[i][j] || ival[i][j+1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
				gradDir[i][j] = 90.0;
			}

			if(( (theta > 112.5) && (theta < 157.5) ) || ( (theta < -22.5) && (theta > -67.5) )) // SE/NW
			{
				if(ival[i-1][j+1] > ival[i][j] || ival[i+1][j-1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
				gradDir[i][j] = 135.0;
			}
			// making threshold image
			if(peaks[i][j] == 1 && getThreshold(ival[i][j]) >= 0)
				thresh[i][j] = 1;
			else thresh[i][j] = 0;

		}// end j for
	}// end i for
} // end candidates

void normalize(double maxVal)
{
	int i, j;
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
			{
				ival[i][j] = (ival[i][j] / maxVal) * 255.0;
			}
	}
	return;
} // end normalize

int inBounds(int i, int j){
	if((i > PICSIZE - 1 || i < 0) && (j > PICSIZE - 1 || j < 0)) return -3;
	if(i > PICSIZE-1 || i < 0) return -2;
	if(j > PICSIZE - 1 || j < 0) return -1;
	return 1;
}

void imageToFile(char *file, int output[256][256]){
	int i, j, flag = 0;
	if(strcmp(file, "cannyXGrad.pgm") == 0 || strcmp(file, "cannyYGrad.pgm") == 0) flag = 1;
	FILE *fOut = fopen(file, "wb");
	fprintf(fOut,"P5\n256 256\n255\n");  /* Print out the .pgm header */
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
			{
				if(flag) fprintf(fOut,"%c", ((char)output[i][j]));
				else fprintf(fOut,"%c", ((char)(output[i][j] * 255)));
			}
	}
	fclose(fOut);
} // end imageToFile

void magToFile(char *file, double output[256][256])
{
	int i, j;

	FILE *fOut = fopen(file, "wb");
	fprintf(fOut,"P5\n256 256\n255\n");  /* Print out the .pgm header */
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
			{
				fprintf(fOut,"%c", ((char)((int)output[i][j])));
			}
	}
} // end imageToFile
