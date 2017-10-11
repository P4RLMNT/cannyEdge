#include <stdio.h>                          /* Sobel.c */
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define  PICSIZE 256
#define  MAXMASK 100
#define  SIGONE 1
#define  SIGTWO 2

void imageToFile(char*, int[256][256]);
void isCand(double image[][PICSIZE]);
void nonMaxSup();
void normalize(double maxVal);



int    pic[PICSIZE][PICSIZE];
int outpicx[PICSIZE][PICSIZE];
int outpicy[PICSIZE][PICSIZE];
int    peaks[PICSIZE][PICSIZE]; // non-maximum suppression
double xmask[MAXMASK][MAXMASK], ymask[MAXMASK][MAXMASK];
double conv[PICSIZE][PICSIZE];
double ival[PICSIZE][PICSIZE];
double gradDir[PICSIZE][PICSIZE];
int lowOutPic[PICSIZE][PICSIZE];
int highOutPic[PICSIZE][PICSIZE];
int final[PICSIZE][PICSIZE];
double hThresh = 110.0, lThresh = 40.0;


int main(int argc, char** argv){
	int     i,j,p,q,s,t,mr,centx,centy;
	double  xmaskval,ymaskval,sum1,sum2,sig,maxival,minival,maxval,ZEROTOL;
	FILE    *fo1, *fo2,*fp1, *fopen();
	char    *foobar, *pgmType;
	int picLength, picWidth;

	argc--; argv++;
	foobar = *argv;
	fp1=fopen(foobar,"rb");


	fo1=fopen("cannyMag.pgm","wb");


	fo2=fopen("canny2output.pgm","wb");

	argc--; argv++;
	foobar = *argv;
	sig = atof(foobar);

	mr = (int)(sig * 3);
	centx = (MAXMASK / 2);
	centy = (MAXMASK / 2);

	fscanf(fp1, "%s", pgmType);
	fscanf(fp1, "%d", &picLength);
	fscanf(fp1, "%d", &picWidth);
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
					{
						pic[i][j]  =  getc (fp1);
					}
	}
// convolution to create the x and y mask
	for (p=-mr;p<=mr;p++)
	{  for (q=-mr;q<=mr;q++)
		 {
				xmaskval = (p * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				ymaskval = (q * (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
				(xmask[p+centy][q+centx]) = xmaskval;
				(ymask[p+centy][q+centx]) = ymaskval;
		 }
	}
	for (i=mr;i<=255-mr;i++)
	{ for (j=mr;j<=255-mr;j++)
		{
			 sum1 = 0;
			 sum2 = 0;
			 for (p=-mr;p<=mr;p++)
			 {
					for (q=-mr;q<=mr;q++)
					{
						 sum1 += pic[i+p][j+q] * ymask[p+centy][q+centx];
						 sum2 += pic[i+p][j+q] * xmask[p+centy][q+centx];
					}
			 }
			 outpicx[i][j] = (int)sum1;
			 outpicy[i][j] = (int)sum2;
			 gradDir[i][j] = atan(sum2/sum1);
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
// printing magnitude imageToFile

fprintf(fo1,"P5\n256 256\n255\n");  /* Print out the .pgm header */
for (i=0;i<256;i++)
{ for (j=0;j<256;j++)
		{
		 // applying high threshold
		 if(ival[i][j] > hThresh) highOutPic[i][j] = 255;
		 else highOutPic[i][j] = 0;
		 // applying low threshold
		 if(ival[i][j] > lThresh) lowOutPic[i][j] = 255;
		 else lowOutPic[i][j] = 0;
		 fprintf(fo1,"%c",(char)((int)(ival[i][j])));

		}
}

imageToFile("cannyXGrad.pgm", outpicx);
imageToFile("cannyYGrad.pgm", outpicy);
imageToFile("cannyPeaks.pgm", peaks);

 return 0;
} // end main

void isCand(double image[][PICSIZE])
{

}

void nonMaxSup()
{
	int i, j;
	double theta;
	for(i = 0; i < PICSIZE; i++) // zeroing edges of peaks image
	{
		peaks[i][0] = 0; // left edge
		peaks[i][PICSIZE-1] = 0; // right edge
	}
	for(j = 0; j < PICSIZE; j++)
	{
		peaks[0][j] = 0; // top edge
		peaks[PICSIZE-1][j] = 0; // bottom edge
	}
	for(i = 1; i < PICSIZE-1; i++)
	{
		for(j = 1; j < PICSIZE-1; j++)
		{
			theta = gradDir[i][j];
			if(theta > 3*M_PI/8 && theta <= -3*M_PI/8) // NS
			{
				if(ival[i][j-1] > ival[i][j] || ival[i][j+1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
			}
			else if(theta > M_PI/8 && theta <= 3*M_PI/8) // NE/SW
			{
				if(ival[i-1][j-1] > ival[i][j] || ival[i+1][j+1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
			}
			else if(theta > -1*M_PI/8 && theta <= M_PI/8) // EW
			{
				if(ival[i-1][j] > ival[i][j] || ival[i+1][j] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
			}
			else if(theta > -1*3*M_PI/8 && theta <= -1*M_PI/8) // SE/NW
			{
				if(ival[i-1][j+1] > ival[i][j] || ival[i+1][j-1] > ival[i][j]) peaks[i][j] = 0;
				else peaks[i][j] = 1;
			}

		}// end j for
	}// end i for
} // end candidates

int getDirection(){
	int result = 0;


	return result;
}

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
}

void imageToFile(char *file, int output[256][256]){
	int i, j, isPeaks = 0;
	if(strcmp(file, "cannyPeaks.pgm") == 0) isPeaks = 1;
	FILE *fOut = fopen(file, "wb");
	fprintf(fOut,"P5\n256 256\n255\n");  /* Print out the .pgm header */
	for (i=0;i<256;i++)
	{ for (j=0;j<256;j++)
			{
				if(isPeaks) fprintf(fOut,"%c", ((char)(output[i][j] * 255)));
				else fprintf(fOut,"%c", ((char)output[i][j]));
			}
	}
} // end imageToFile
